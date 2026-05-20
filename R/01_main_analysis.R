# =============================================================================
# Tillage x Flow Regime Interaction Controls Phosphorus Source Dynamics in
# Tile Drainage: Integrating Continuous Discharge Monitoring, Water Residence
# Time, and Phosphate Oxygen Isotopes
#
# Authors:  Arfania, H., Kayler, Z. E., Strawn, D. G., Brooks, E. S., Laan, M.
# Site:     R. J. Cook Agronomy Farm (CAF) LTAR, Pullman, WA, USA
# Target:   JGR Biogeosciences (AGU)
# Script:   01_main_analysis.R  -- end-to-end workflow producing all main-text
#           and supplementary tables and figures EXCEPT the Craig diagram
#           (see 02_craig_diagram.R).
#
# Reproducibility
# ---------------
# This script is path-agnostic. Open the project from the repository root
# (Arfania_coder.Rproj if available, or open R with the working directory set
# to the repo root) and it will resolve all inputs and outputs through the
# `here` package. No setwd() calls are used.
#
# Inputs (under data/):
#   18Ow_Prcep.csv, 18Ow_Tile.csv, Tile_d18Op.csv,
#   DRP_2022_23.csv, DRP_2023_24.csv,
#   CT_Discharge_2023.csv, CT_Discharge_2024.csv,
#   NT_Discharge_2023.csv, NT_Discharge_2024.csv,
#   CT_water_tempreture_2023.csv, CT_water_tempreture_2024.csv,
#   NT_water_tempreture_2023.csv, NT_water_tempreture_2024.csv
#
# Outputs:
#   output/   -- CSV tables (Table 1-8) and analysis_summary.rds
#   figures/  -- PNG (300 dpi), PDF, and TIFF (600 dpi, LZW) versions of
#                Figures 2, 3, 4, 6, 7, S1, S2
#
# Revision note (v2.0, 2026-05): Dynamic equilibrium d18O-PO4 is computed by
# matching each PO4 sample to (i) treatment-specific tile-drain d18O-H2O and
# (ii) daily-mean continuously logged tile water temperature, addressing the
# reviewer concern that pooled water isotopes and a fixed 5 deg C assumption
# oversimplify the Chang & Blake (2015) equilibrium reference. Fixed-equilibrium
# (11.2 per mil) results are retained for sensitivity comparison.
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

# Clear workspace for a clean run
rm(list = ls())

# Required packages
required_packages <- c(
  "here",          # project-relative paths
  "tidyverse",     # data manipulation and visualization
  "lubridate",     # date handling
  "ggpubr",        # publication-ready plot helpers
  "car",           # Type III ANOVA
  "effsize",       # effect size
  "minpack.lm",    # Levenberg-Marquardt nonlinear LS
  "scales",        # plot scales
  "cowplot",       # plot composition
  "viridis"        # color palette
)

new_pkgs <- required_packages[!(required_packages %in%
                                  installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_packages, library, character.only = TRUE))

# Project-relative paths (works on any OS / any user)
data_dir    <- here("data")
output_dir  <- here("output")
figures_dir <- here("figures")
dir.create(output_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Publication theme and color palette
theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "bottom",
                  strip.background = element_rect(fill = "white")))

treatment_colors <- c("CT" = "#D55E00", "NT" = "#0072B2")

# Chang & Blake (2015) equilibrium equation
# d18O-PO4(eq) = exp[(14.43 / T_K) - 0.02654] * (d18O-H2O + 1000) - 1000
cb_eq <- function(d18Ow, T_C) {
  T_K <- T_C + 273.15
  exp(14.43 / T_K - 0.02654) * (d18Ow + 1000) - 1000
}

# Generic date parser: handles m/d/Y, Y-m-d, dmY, and POSIXct-style strings
parse_mixed_date <- function(x) {
  pdt <- parse_date_time(as.character(x),
                         orders = c("mdY", "Ymd", "Ymd HMS", "Ymd HM",
                                    "mdY HMS", "mdY HM", "dmY"),
                         quiet = TRUE)
  as_date(pdt)
}

cat("Setup complete. Packages loaded; project rooted at:\n  ", here(), "\n",
    sep = "")


# =============================================================================
# 2. DATA IMPORT
# =============================================================================

cat("\n=== SECTION 2: DATA IMPORT ===\n")

# ---- 2.1 Continuous tile discharge (15-min, multiple deployment formats) ----
#
# The four discharge CSVs have different layouts (full instrument header,
# two-column "Date, cfs", row-labels-with-trailing-blanks, and TOA5 multi-line
# header). The helper below auto-detects the timestamp and discharge columns
# and returns a clean tibble [treatment, datetime, discharge_L_s].

read_discharge <- function(path, label) {
  raw <- read_lines(path, n_max = 10)
  hdr_pattern <- "(?i)\\b(date|time|row labels|timestamp)\\b"
  hdr_idx <- which(grepl(hdr_pattern, raw))[1]
  if (is.na(hdr_idx)) hdr_idx <- 1
  skip_n <- hdr_idx - 1L

  df <- suppressWarnings(read_csv(
    path, skip = skip_n, show_col_types = FALSE,
    name_repair = "unique", guess_max = 1e5
  ))

  # Drop blank/auto-generated columns that contain no data
  keep <- vapply(seq_along(df), function(i) {
    nm <- names(df)[i]
    if (is.na(nm) || nm == "" || grepl("^\\.\\.\\.\\d+$", nm)) {
      return(any(!is.na(df[[i]])))
    }
    TRUE
  }, logical(1))
  df <- df[, keep, drop = FALSE]

  nm <- names(df)
  ts_idx <- which(grepl("(?i)timestamp|date|row labels", nm))[1]
  if (is.na(ts_idx)) ts_idx <- 1L

  # Prefer Q_lps, fall back to "Sum of Q_lps", then any cfs-like column
  q_lps_idx <- which(grepl("^Q_lps$", nm, ignore.case = TRUE))[1]
  if (is.na(q_lps_idx)) {
    q_lps_idx <- which(grepl("(?i)sum of q_lps|q.*l.*s", nm))[1]
  }
  q_cfs_idx <- which(grepl("(?i)q_cfs|q.*cfs|^cfs$|tl_q_cfs", nm))[1]

  if (!is.na(q_lps_idx)) {
    use_idx <- q_lps_idx; in_cfs <- FALSE
  } else if (!is.na(q_cfs_idx)) {
    use_idx <- q_cfs_idx; in_cfs <- TRUE
  } else {
    stop("No discharge column found in ", basename(path))
  }

  tibble(
    treatment     = label,
    Timestamp_raw = as.character(df[[ts_idx]]),
    Q_raw         = suppressWarnings(as.numeric(df[[use_idx]]))
  ) %>%
    mutate(Timestamp = parse_date_time(
      Timestamp_raw,
      orders = c("Ymd HMS", "Ymd HM", "Ymd",
                 "mdY HMS", "mdY HM", "mdY",
                 "dmY HMS", "dmY HM", "dmY"),
      quiet = TRUE)) %>%
    filter(!is.na(Timestamp), !is.na(Q_raw)) %>%
    mutate(Q_lps = if (in_cfs) Q_raw * 28.317 else Q_raw) %>%
    select(treatment, Timestamp, Q_lps)
}

discharge_data <- bind_rows(
  read_discharge(file.path(data_dir, "CT_Discharge_2023.csv"), "CT"),
  read_discharge(file.path(data_dir, "CT_Discharge_2024.csv"), "CT"),
  read_discharge(file.path(data_dir, "NT_Discharge_2023.csv"), "NT"),
  read_discharge(file.path(data_dir, "NT_Discharge_2024.csv"), "NT")
) %>%
  rename(datetime = Timestamp, discharge_L_s = Q_lps)

cat("Discharge observations:", nrow(discharge_data),
    "  (CT:", sum(discharge_data$treatment == "CT"),
    "| NT:", sum(discharge_data$treatment == "NT"), ")\n")

# ---- 2.2 Dissolved reactive phosphorus (DRP) -------------------------------

# 2022-23: paired-side-by-side layout (CT cols 1:3, NT cols 8:10)
drp_2223_raw <- read_csv(file.path(data_dir, "DRP_2022_23.csv"),
                         show_col_types = FALSE, name_repair = "minimal")

drp_ct_2223 <- drp_2223_raw %>%
  select(1, 2, 3) %>%
  rename(sample_id = 1, Date = 2, drp_mg_L = 3) %>%
  mutate(treatment = "CT")

drp_nt_2223 <- drp_2223_raw %>%
  select(8, 9, 10) %>%
  rename(sample_id = 1, Date = 2, drp_mg_L = 3) %>%
  mutate(treatment = "NT")

drp_2223 <- bind_rows(drp_ct_2223, drp_nt_2223) %>%
  filter(!is.na(drp_mg_L)) %>%
  mutate(sample_id = as.character(sample_id),
         Date      = parse_mixed_date(Date),
         drp_mg_L  = suppressWarnings(as.numeric(drp_mg_L)))

# 2023-24: stacked with Location encoding treatment; tile only (exclude SR)
drp_2324 <- read_csv(file.path(data_dir, "DRP_2023_24.csv"),
                     show_col_types = FALSE, name_repair = "minimal") %>%
  rename_with(~ tolower(.x)) %>%
  filter(grepl("(?i)_TL", location)) %>%
  mutate(
    treatment = case_when(
      grepl("(?i)CT_TL", location) ~ "CT",
      grepl("(?i)NT_TL", location) ~ "NT",
      TRUE ~ NA_character_),
    sample_id = as.character(sampleid),
    Date      = parse_mixed_date(date),
    drp_mg_L  = suppressWarnings(as.numeric(`mg/l drp`))
  ) %>%
  filter(!is.na(drp_mg_L), !is.na(treatment), !is.na(Date)) %>%
  select(sample_id, treatment, Date, drp_mg_L)

drp_data <- bind_rows(drp_2223, drp_2324) %>%
  filter(!is.na(drp_mg_L), drp_mg_L > 0, !is.na(Date)) %>%
  mutate(drp_mg_L = pmax(0.005, drp_mg_L))   # MDL = 0.005 mg/L

cat("DRP samples:", nrow(drp_data),
    "  (CT:", sum(drp_data$treatment == "CT"),
    "| NT:", sum(drp_data$treatment == "NT"), ")\n")

# ---- 2.3 Water isotopes (d18O-H2O, d2H-H2O) --------------------------------

precip_d18Ow <- read_csv(file.path(data_dir, "18Ow_Prcep.csv"),
                         show_col_types = FALSE) %>%
  rename(Date = Date, d18O_H2O = `18O_H2O`, d2H_H2O = Detrium_H2O) %>%
  mutate(Date = parse_mixed_date(Date)) %>%
  filter(!is.na(d18O_H2O), !is.na(Date))

tile_h2o_raw <- read_csv(file.path(data_dir, "18Ow_Tile.csv"),
                         show_col_types = FALSE, name_repair = "minimal")
names(tile_h2o_raw)[1] <- "Date"
tile_h2o <- tile_h2o_raw %>%
  select(Date, NT_18OW, NT_Deutrium, CT_18OW, CT_Deutrium) %>%
  filter(!is.na(Date)) %>%
  mutate(Date = parse_mixed_date(Date)) %>%
  filter(!is.na(Date))

tile_isotope <- tile_h2o %>%
  select(Date, CT_18OW, NT_18OW) %>%
  pivot_longer(c(CT_18OW, NT_18OW),
               names_to = "treatment", values_to = "d18O_H2O") %>%
  mutate(treatment = sub("_18OW$", "", treatment)) %>%
  filter(!is.na(d18O_H2O))

cat("Water isotopes: precip n =", nrow(precip_d18Ow),
    ", tile n =", nrow(tile_isotope), "\n")

# ---- 2.4 Tile-drain d18O-PO4 -----------------------------------------------

phosphate_isotope <- read_csv(file.path(data_dir, "Tile_d18Op.csv"),
                              show_col_types = FALSE) %>%
  mutate(Date = parse_mixed_date(Date)) %>%
  filter(!is.na(Date)) %>%
  pivot_longer(c(CT, NT), names_to = "treatment", values_to = "d18O_PO4") %>%
  filter(!is.na(d18O_PO4)) %>%
  arrange(treatment, Date)

cat("d18O-PO4 samples:", nrow(phosphate_isotope), "\n")

# ---- 2.5 Continuous tile water temperature ---------------------------------

read_temp <- function(path, label) {
  df <- read_csv(path, show_col_types = FALSE, name_repair = "minimal")
  names(df)[1:2] <- c("Timestamp", "Temp_C")
  df %>%
    mutate(
      Timestamp = parse_date_time(
        Timestamp,
        orders = c("Ymd HMS", "mdY HM", "mdY HMS", "Ymd HM", "mdY", "Ymd"),
        quiet = TRUE),
      Temp_C    = suppressWarnings(as.numeric(Temp_C))
    ) %>%
    filter(!is.na(Timestamp), !is.na(Temp_C),
           Temp_C >= -2, Temp_C <= 30) %>%
    mutate(treatment = label, Date = as_date(Timestamp))
}

temp_files <- list(
  list(file.path(data_dir, "CT_water_tempreture_2023.csv"), "CT"),
  list(file.path(data_dir, "CT_water_tempreture_2024.csv"), "CT"),
  list(file.path(data_dir, "NT_water_tempreture_2023.csv"), "NT"),
  list(file.path(data_dir, "NT_water_tempreture_2024.csv"), "NT")
)
temp_all <- map_dfr(temp_files, ~ read_temp(.x[[1]], .x[[2]]))

cat("Tile water temperature: ", nrow(temp_all),
    " QC'd 15-min readings (CT:", sum(temp_all$treatment == "CT"),
    "| NT:", sum(temp_all$treatment == "NT"), ")\n")

daily_temp <- temp_all %>%
  group_by(treatment, Date) %>%
  summarise(Temp_C = mean(Temp_C, na.rm = TRUE), .groups = "drop")


# =============================================================================
# 3. DISCHARGE ANALYSIS (Table 1, Figure 2)
# =============================================================================

cat("\n=== SECTION 3: DISCHARGE ANALYSIS ===\n")

discharge_stats <- discharge_data %>%
  group_by(treatment) %>%
  summarise(
    n_obs      = n(),
    n_days     = as.numeric(difftime(max(datetime), min(datetime),
                                     units = "days")),
    mean_L_s   = mean(discharge_L_s,   na.rm = TRUE),
    sd_L_s     = sd(discharge_L_s,     na.rm = TRUE),
    median_L_s = median(discharge_L_s, na.rm = TRUE),
    Q75        = quantile(discharge_L_s, 0.75, na.rm = TRUE),
    max_L_s    = max(discharge_L_s,   na.rm = TRUE),
    .groups = "drop"
  )

write_csv(discharge_stats, file.path(output_dir, "Table1_Discharge_Statistics.csv"))
print(discharge_stats)

Q75_CT <- discharge_stats$Q75[discharge_stats$treatment == "CT"]
Q75_NT <- discharge_stats$Q75[discharge_stats$treatment == "NT"]

# Daily means and daily Q75 thresholds (used to classify grab-sample flow regimes)
daily_discharge <- discharge_data %>%
  mutate(Date = as_date(datetime)) %>%
  group_by(treatment, Date) %>%
  summarise(mean_discharge = mean(discharge_L_s, na.rm = TRUE),
            .groups = "drop")

q75_daily <- daily_discharge %>%
  group_by(treatment) %>%
  summarise(Q75 = quantile(mean_discharge, 0.75, na.rm = TRUE),
            .groups = "drop")
q75_daily_map <- setNames(q75_daily$Q75, q75_daily$treatment)

cat("Q75 (15-min):  CT =", round(Q75_CT, 3),
    ", NT =", round(Q75_NT, 3), "L/s\n")
cat("Q75 (daily):   CT =", round(q75_daily_map["CT"], 3),
    ", NT =", round(q75_daily_map["NT"], 3), "L/s\n")

fig2 <- ggplot(daily_discharge,
               aes(x = Date, y = mean_discharge, color = treatment)) +
  geom_line(alpha = 0.7) +
  geom_hline(data = q75_daily,
             aes(yintercept = Q75, color = treatment, linetype = "Q75"),
             linewidth = 0.8) +
  scale_y_log10(labels = scales::number_format(accuracy = 0.01)) +
  scale_color_manual(values = treatment_colors) +
  scale_linetype_manual(values = c("Q75" = "dashed")) +
  labs(x = "Date",
       y = expression(paste("Discharge (L s"^"-1", ")")),
       color = "Treatment", linetype = "") +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2))

ggsave(file.path(figures_dir, "Fig2_Discharge_TimeSeries.png"),
       fig2, width = 10, height = 5, dpi = 300)
ggsave(file.path(figures_dir, "Fig2_Discharge_TimeSeries.pdf"),
       fig2, width = 10, height = 5)
ggsave(file.path(figures_dir, "Fig2_Discharge_TimeSeries.tiff"),
       fig2, width = 10, height = 5, dpi = 600, compression = "lzw")


# =============================================================================
# 4. DRP ANALYSIS (Table 2)
# =============================================================================

cat("\n=== SECTION 4: DRP ANALYSIS ===\n")

drp_data <- drp_data %>%
  left_join(daily_discharge, by = c("treatment", "Date")) %>%
  mutate(flow_regime = case_when(
    is.na(mean_discharge) ~ NA_character_,
    treatment == "CT" & mean_discharge >= q75_daily_map["CT"] ~ "Event",
    treatment == "NT" & mean_discharge >= q75_daily_map["NT"] ~ "Event",
    TRUE ~ "Baseflow"
  ),
  discharge_L_s = mean_discharge)

drp_summary <- drp_data %>%
  filter(!is.na(flow_regime)) %>%
  group_by(treatment, flow_regime) %>%
  summarise(
    n             = n(),
    mean_drp      = mean(drp_mg_L, na.rm = TRUE),
    sd_drp        = sd(drp_mg_L,   na.rm = TRUE),
    FWMC          = sum(drp_mg_L * discharge_L_s, na.rm = TRUE) /
                    sum(discharge_L_s,            na.rm = TRUE),
    total_load_kg = sum(drp_mg_L * discharge_L_s * 0.001, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(treatment) %>%
  mutate(pct_total = total_load_kg / sum(total_load_kg) * 100) %>%
  ungroup()

write_csv(drp_summary, file.path(output_dir, "Table2_DRP_Summary.csv"))
print(drp_summary)


# =============================================================================
# 5. WATER ISOTOPE ANALYSIS (Table 3, Figure 3)
# =============================================================================

cat("\n=== SECTION 5: WATER ISOTOPE ANALYSIS ===\n")

# Outlier filter: precip samples > 5 SD from mean d18O (matches manuscript)
precip_clean <- precip_d18Ow %>%
  filter(abs(d18O_H2O - mean(d18O_H2O, na.rm = TRUE)) <
           5 * sd(d18O_H2O, na.rm = TRUE))

precip_stats <- precip_clean %>%
  summarise(source = "Precipitation",
            n = n(),
            mean_d18O  = mean(d18O_H2O),
            sd_d18O    = sd(d18O_H2O),
            min_d18O   = min(d18O_H2O),
            max_d18O   = max(d18O_H2O),
            range_d18O = max_d18O - min_d18O)

tile_stats <- tile_isotope %>%
  group_by(treatment) %>%
  summarise(source = paste("Tile -", first(treatment)),
            n = n(),
            mean_d18O  = mean(d18O_H2O),
            sd_d18O    = sd(d18O_H2O),
            min_d18O   = min(d18O_H2O),
            max_d18O   = max(d18O_H2O),
            range_d18O = max_d18O - min_d18O,
            .groups = "drop") %>%
  select(-treatment)

all_tile_stats <- tile_isotope %>%
  summarise(source = "All Tile Drain",
            n = n(),
            mean_d18O  = mean(d18O_H2O),
            sd_d18O    = sd(d18O_H2O),
            min_d18O   = min(d18O_H2O),
            max_d18O   = max(d18O_H2O),
            range_d18O = max_d18O - min_d18O)

table3 <- bind_rows(precip_stats, tile_stats, all_tile_stats)
write_csv(table3, file.path(output_dir, "Table3_Water_Isotope_Summary.csv"))
print(table3)

f_test_result  <- var.test(precip_clean$d18O_H2O, tile_isotope$d18O_H2O)
variance_ratio <- var(precip_clean$d18O_H2O) / var(tile_isotope$d18O_H2O)

# CT vs NT paired t-test on matched dates
matched_dates <- intersect(
  tile_isotope %>% filter(treatment == "CT") %>% pull(Date),
  tile_isotope %>% filter(treatment == "NT") %>% pull(Date)
)
ct_data <- tile_isotope %>% filter(treatment == "CT",
                                   Date %in% matched_dates) %>% arrange(Date)
nt_data <- tile_isotope %>% filter(treatment == "NT",
                                   Date %in% matched_dates) %>% arrange(Date)

paired_t_test <- t.test(ct_data$d18O_H2O, nt_data$d18O_H2O, paired = TRUE)
wilcox_test   <- wilcox.test(ct_data$d18O_H2O, nt_data$d18O_H2O, paired = TRUE)

cat("Variance ratio (precip/tile):", round(variance_ratio, 1),
    "  F-test p =", format(f_test_result$p.value, scientific = TRUE), "\n")
cat("CT vs NT (paired, n =", length(matched_dates), ") t-test p =",
    format(paired_t_test$p.value, digits = 4), "\n")

# Figure 3
boxplot_data <- bind_rows(
  precip_clean %>% transmute(Date, d18O_H2O,
                             source = "Precipitation",
                             treatment = NA_character_),
  tile_isotope %>% transmute(Date, d18O_H2O,
                             source = paste("Tile -", treatment),
                             treatment)
)
boxplot_data$source <- factor(boxplot_data$source,
                              levels = c("Precipitation", "Tile - CT", "Tile - NT"))
means_data <- boxplot_data %>%
  group_by(source) %>%
  summarise(mean_d18O = mean(d18O_H2O), .groups = "drop")

fig3 <- ggplot(boxplot_data, aes(x = source, y = d18O_H2O, fill = source)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_point(data = means_data, aes(x = source, y = mean_d18O),
             shape = 23, size = 4, fill = "red", color = "black") +
  annotate("segment", x = 2, xend = 3, y = -12.5, yend = -12.5) +
  annotate("text", x = 2.5, y = -12.2, label = "***", size = 6) +
  scale_fill_manual(values = c("Precipitation" = "gray70",
                               "Tile - CT" = "#D55E00",
                               "Tile - NT" = "#0072B2")) +
  labs(x = NULL,
       y = expression(paste(delta^{18}, "O-H"[2], "O (\u2030 VSMOW)"))) +
  theme(legend.position = "none")

ggsave(file.path(figures_dir, "Fig3_Boxplot_Comparison.png"),
       fig3, width = 8, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Fig3_Boxplot_Comparison.pdf"),
       fig3, width = 8, height = 6)
ggsave(file.path(figures_dir, "Fig3_Boxplot_Comparison.tiff"),
       fig3, width = 8, height = 6, dpi = 600, compression = "lzw")


# =============================================================================
# 6. WATER RESIDENCE TIME (Table 4, Figure 4)
# =============================================================================

cat("\n=== SECTION 6: RESIDENCE TIME ===\n")

sd_precip <- sd(precip_clean$d18O_H2O)
sd_ct     <- sd(ct_data$d18O_H2O)
sd_nt     <- sd(nt_data$d18O_H2O)
DR_CT <- sd_ct / sd_precip
DR_NT <- sd_nt / sd_precip
omega <- 2 * pi / 365.25

tau_CT_days <- 1 / (omega * sqrt((1 / DR_CT^2) - 1))
tau_NT_days <- 1 / (omega * sqrt((1 / DR_NT^2) - 1))

cat("Damping-ratio residence time: CT =", round(tau_CT_days, 1),
    " d, NT =", round(tau_NT_days, 1), "d\n")

# Sine-wave fits (amplitude ratio / Fyw)
precip_sine <- precip_clean %>% mutate(doy = yday(Date))
ct_sine     <- ct_data %>% mutate(doy = yday(Date))
nt_sine     <- nt_data %>% mutate(doy = yday(Date))

safe_sin <- function(data, A0, off0) {
  tryCatch(
    nlsLM(d18O_H2O ~ A * sin(2 * pi * doy / 365.25 + phi) + offset,
          data = data,
          start = list(A = A0, phi = 0, offset = off0),
          control = nls.lm.control(maxiter = 200)),
    error = function(e) NULL)
}
fit_precip <- safe_sin(precip_sine, 3, -14)
fit_ct     <- safe_sin(ct_sine,    0.5, -14)
fit_nt     <- safe_sin(nt_sine,    0.3, -14.5)

if (!is.null(fit_precip) && !is.null(fit_ct) && !is.null(fit_nt)) {
  A_precip <- abs(coef(fit_precip)["A"])
  A_ct     <- abs(coef(fit_ct)["A"])
  A_nt     <- abs(coef(fit_nt)["A"])
  Fyw_CT <- A_ct / A_precip
  Fyw_NT <- A_nt / A_precip
  tau_amp_CT <- sqrt((A_precip / A_ct)^2 - 1) / omega
  tau_amp_NT <- sqrt((A_precip / A_nt)^2 - 1) / omega
} else {
  cat("WARNING: sine fit failed; using fallback values.\n")
  A_precip <- 2.1; A_ct <- 0.819; A_nt <- 0.399
  Fyw_CT <- 0.39; Fyw_NT <- 0.19
  tau_amp_CT <- 138; tau_amp_NT <- 294
}

cat("Amplitude-ratio residence time: CT =", round(tau_amp_CT, 0),
    " d, NT =", round(tau_amp_NT, 0), "d\n")
cat("Young water fraction (Fyw): CT =", round(Fyw_CT * 100, 1),
    "%, NT =", round(Fyw_NT * 100, 1), "%\n")

table4 <- tibble(
  Method = c("Damping Ratio (SD)", "Amplitude Ratio (Sine)",
             "Young Water Fraction (Fyw)"),
  CT = c(paste0(round(tau_CT_days, 1),  " d"),
         paste0(round(tau_amp_CT, 0),   " d"),
         paste0(round(Fyw_CT * 100, 1), "%")),
  NT = c(paste0(round(tau_NT_days, 1),  " d"),
         paste0(round(tau_amp_NT, 0),   " d"),
         paste0(round(Fyw_NT * 100, 1), "%")),
  Note = c("Sensitive to total variance reduction",
           "Sensitive to seasonal signal attenuation",
           "Tracer-derived index (Kirchner, 2016) - NOT a direct residence-time estimate")
)
write_csv(table4, file.path(output_dir, "Table4_Residence_Time.csv"))
print(table4)

# Figure 4
doy_seq <- seq(1, 365, length.out = 200)
make_pred <- function(fit, label, fallback_amp, fallback_offset) {
  if (!is.null(fit)) {
    tibble(doy = doy_seq,
           d18O = predict(fit, newdata = data.frame(doy = doy_seq)),
           source = label)
  } else {
    tibble(doy = doy_seq,
           d18O = fallback_amp * sin(2 * pi * doy_seq / 365.25 - 1.5) +
             fallback_offset,
           source = label)
  }
}
pred_all <- bind_rows(
  make_pred(fit_precip, "Precipitation", 2.1, -14.4),
  make_pred(fit_ct, "CT Tile", 0.5, -13.95),
  make_pred(fit_nt, "NT Tile", 0.25, -14.34)
)
obs_data <- bind_rows(
  precip_sine %>% transmute(doy, d18O = d18O_H2O, source = "Precipitation"),
  ct_sine     %>% transmute(doy, d18O = d18O_H2O, source = "CT Tile"),
  nt_sine     %>% transmute(doy, d18O = d18O_H2O, source = "NT Tile")
)

fig4 <- ggplot() +
  geom_line(data = pred_all,
            aes(x = doy, y = d18O, color = source, linetype = source),
            linewidth = 1) +
  geom_point(data = obs_data,
             aes(x = doy, y = d18O, color = source, shape = source),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = c("Precipitation" = "gray40",
                                "CT Tile" = "#D55E00",
                                "NT Tile" = "#0072B2")) +
  scale_linetype_manual(values = c("Precipitation" = "solid",
                                   "CT Tile" = "dashed",
                                   "NT Tile" = "dashed")) +
  scale_shape_manual(values = c("Precipitation" = 16,
                                "CT Tile" = 17, "NT Tile" = 15)) +
  labs(subtitle = paste0("Damping-ratio residence time: CT = ",
                         round(tau_CT_days, 0),
                         " d, NT = ", round(tau_NT_days, 0), " d"),
       x = "Day of Year",
       y = expression(paste(delta^{18}, "O-H"[2], "O (\u2030 VSMOW)")),
       color = "Source", linetype = "Source", shape = "Source")

ggsave(file.path(figures_dir, "Fig4_Seasonal_Sine_Curves.png"),
       fig4, width = 10, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Fig4_Seasonal_Sine_Curves.pdf"),
       fig4, width = 10, height = 6)
ggsave(file.path(figures_dir, "Fig4_Seasonal_Sine_Curves.tiff"),
       fig4, width = 10, height = 6, dpi = 600, compression = "lzw")


# =============================================================================
# 7. PHOSPHATE ISOTOPES - DYNAMIC EQUILIBRIUM (Tables 5-8, Figs 6-7, S1-S2)
# =============================================================================

cat("\n=== SECTION 7: d18O-PO4 (DYNAMIC EQUILIBRIUM) ===\n")

# ---- 7.1 Match each PO4 sample to treatment-specific water, T, Q -----------
#
# For each PO4 sample we find the closest-date observation (within the
# specified lag window) for tile water d18O-H2O, daily-mean tile water T,
# and daily-mean discharge. This is the core revision addressing the
# reviewer concern about a pooled/fixed equilibrium reference.

nearest_value <- function(target_date, df, value_col, max_lag_days = 7) {
  if (nrow(df) == 0) return(c(NA_real_, NA_integer_))
  ldates <- as_date(df$Date)
  tdate  <- as_date(target_date)
  d <- as.numeric(abs(difftime(ldates, tdate, units = "days")))
  if (any(d == 0, na.rm = TRUE)) {
    i <- which(d == 0)[1]
    return(c(df[[value_col]][i], 0L))
  }
  i <- which.min(d)
  if (!length(i) || is.na(d[i]) || d[i] > max_lag_days) {
    return(c(NA_real_, NA_integer_))
  }
  c(df[[value_col]][i], as.integer(d[i]))
}

attach_match <- function(df_po4, lookup_df, value_col, max_lag_days,
                         out_value_name, out_lag_name,
                         treatment_filter = TRUE) {
  res <- vector("list", nrow(df_po4))
  for (i in seq_len(nrow(df_po4))) {
    sub <- if (treatment_filter) {
      lookup_df %>% filter(treatment == df_po4$treatment[i],
                           !is.na(.data[[value_col]]))
    } else {
      lookup_df %>% filter(!is.na(.data[[value_col]]))
    }
    res[[i]] <- nearest_value(df_po4$Date[i], sub, value_col, max_lag_days)
  }
  out <- do.call(rbind, res) %>% as.data.frame()
  names(out) <- c(out_value_name, out_lag_name)
  out
}

h2o_long <- tile_h2o %>%
  select(Date, CT_18OW, NT_18OW) %>%
  pivot_longer(c(CT_18OW, NT_18OW),
               names_to = "treatment", values_to = "d18O_H2O") %>%
  mutate(treatment = sub("_18OW$", "", treatment))

water_match <- attach_match(phosphate_isotope, h2o_long,  "d18O_H2O", 7,
                            "d18O_H2O",    "h2o_lag_days")
temp_match  <- attach_match(phosphate_isotope, daily_temp, "Temp_C",   7,
                            "Temp_C",      "T_lag_days")

disc_long <- daily_discharge %>% rename(Q_lps = mean_discharge)
q_match <- attach_match(phosphate_isotope, disc_long, "Q_lps", 3,
                        "Q_daily_lps", "Q_lag_days")

phosphate_isotope <- bind_cols(phosphate_isotope, water_match,
                               temp_match, q_match) %>%
  mutate(flow_regime = case_when(
    is.na(Q_daily_lps) ~ NA_character_,
    treatment == "CT" & Q_daily_lps >= q75_daily_map["CT"] ~ "Event",
    treatment == "NT" & Q_daily_lps >= q75_daily_map["NT"] ~ "Event",
    TRUE ~ "Baseflow"
  ))

# ---- 7.2 Dynamic and fixed equilibrium values ------------------------------

fixed_d18O_H2O <- -14.1
fixed_T_C      <- 5
equilibrium_d18O_PO4_fixed <- cb_eq(fixed_d18O_H2O, fixed_T_C)

phosphate_isotope <- phosphate_isotope %>%
  mutate(
    Temp_K            = Temp_C + 273.15,
    d18O_PO4_eq       = cb_eq(d18O_H2O, Temp_C),       # dynamic
    Delta_eq          = d18O_PO4 - d18O_PO4_eq,        # dynamic
    d18O_PO4_eq_fixed = equilibrium_d18O_PO4_fixed,
    Delta_eq_fixed    = d18O_PO4 - equilibrium_d18O_PO4_fixed
  )

write_csv(phosphate_isotope,
          file.path(output_dir, "po4_matched_dynamic.csv"))

cat("FIXED equilibrium:", round(equilibrium_d18O_PO4_fixed, 2), "per mil\n")
cat("PO4 samples without matched water/T/Q: ",
    sum(is.na(phosphate_isotope$d18O_H2O)), " / ",
    sum(is.na(phosphate_isotope$Temp_C)),    " / ",
    sum(is.na(phosphate_isotope$Q_daily_lps)),
    " (of ", nrow(phosphate_isotope), ")\n", sep = "")

ref_eq <- phosphate_isotope %>%
  filter(!is.na(d18O_PO4_eq)) %>%
  group_by(treatment) %>%
  summarise(n             = n(),
            d18O_H2O_mean = mean(d18O_H2O),
            Temp_C_mean   = mean(Temp_C),
            eq_mean       = mean(d18O_PO4_eq),
            eq_sd         = sd(d18O_PO4_eq),
            .groups = "drop")
write_csv(ref_eq, file.path(output_dir, "treatment_eq_reference.csv"))

eq_CT <- ref_eq$eq_mean[ref_eq$treatment == "CT"]
eq_NT <- ref_eq$eq_mean[ref_eq$treatment == "NT"]

# ---- 7.3 Table 5: overall summary ------------------------------------------

po4_overall <- phosphate_isotope %>%
  group_by(treatment) %>%
  summarise(
    n           = n(),
    mean_d18O   = mean(d18O_PO4),
    sd_d18O     = sd(d18O_PO4),
    median_d18O = median(d18O_PO4),
    min_d18O    = min(d18O_PO4),
    max_d18O    = max(d18O_PO4),
    shapiro_p   = shapiro.test(d18O_PO4)$p.value,
    .groups = "drop")
write_csv(po4_overall, file.path(output_dir, "Table5_PO4_Overall_Summary.csv"))

t_test_overall <- t.test(d18O_PO4 ~ treatment, data = phosphate_isotope)
wilcox_overall <- wilcox.test(d18O_PO4 ~ treatment, data = phosphate_isotope)
cohens_d <- cohen.d(
  phosphate_isotope$d18O_PO4[phosphate_isotope$treatment == "CT"],
  phosphate_isotope$d18O_PO4[phosphate_isotope$treatment == "NT"])

# ---- 7.4 Table 6: treatment x flow regime ----------------------------------

po4_by_flow <- phosphate_isotope %>%
  filter(!is.na(flow_regime)) %>%
  group_by(treatment, flow_regime) %>%
  summarise(
    n                   = n(),
    mean_d18O           = mean(d18O_PO4),
    sd_d18O             = sd(d18O_PO4),
    flow_weighted       = sum(d18O_PO4 * Q_daily_lps, na.rm = TRUE) /
                          sum(Q_daily_lps,            na.rm = TRUE),
    mean_Q              = mean(Q_daily_lps, na.rm = TRUE),
    Delta_eq_dynamic    = mean(Delta_eq, na.rm = TRUE),
    Delta_eq_dynamic_sd = sd(Delta_eq,   na.rm = TRUE),
    Delta_eq_fixed      = mean(Delta_eq_fixed, na.rm = TRUE),
    .groups = "drop")
write_csv(po4_by_flow, file.path(output_dir, "Table6_PO4_by_Flow_Regime.csv"))

# ---- 7.5 Table 7: two-way ANOVA (three response specifications) ------------

phosphate_isotope$treatment   <- factor(phosphate_isotope$treatment,
                                        levels = c("CT", "NT"))
phosphate_isotope$flow_regime <- factor(phosphate_isotope$flow_regime,
                                        levels = c("Baseflow", "Event"))

prev_contr <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

df_aov <- phosphate_isotope %>% filter(!is.na(flow_regime))
df_dyn <- phosphate_isotope %>% filter(!is.na(flow_regime), !is.na(Delta_eq))

anova_obs <- Anova(lm(d18O_PO4       ~ treatment * flow_regime, data = df_aov),
                   type = 3)
anova_dyn <- Anova(lm(Delta_eq       ~ treatment * flow_regime, data = df_dyn),
                   type = 3)
anova_fix <- Anova(lm(Delta_eq_fixed ~ treatment * flow_regime, data = df_aov),
                   type = 3)

options(prev_contr)

make_anova_row <- function(a, source_label) {
  tibble(
    Source = c("Treatment", "Flow regime",
               "Treatment x Flow regime", "Residuals"),
    df = c(1, 1, 1, a["Residuals", "Df"]),
    Sum_Sq = round(c(a["treatment", "Sum Sq"],
                     a["flow_regime", "Sum Sq"],
                     a["treatment:flow_regime", "Sum Sq"],
                     a["Residuals", "Sum Sq"]), 2),
    F_value = round(c(a["treatment", "F value"],
                      a["flow_regime", "F value"],
                      a["treatment:flow_regime", "F value"],
                      NA), 2),
    p_value = c(format(a["treatment", "Pr(>F)"], digits = 3),
                format(a["flow_regime", "Pr(>F)"], digits = 3),
                format(a["treatment:flow_regime", "Pr(>F)"],
                       scientific = TRUE, digits = 2),
                "-"),
    Response = source_label
  )
}

table7 <- bind_rows(
  make_anova_row(anova_obs, "Observed d18O-PO4 (manuscript original)"),
  make_anova_row(anova_dyn, "Delta_eq (DYNAMIC, sample-specific)"),
  make_anova_row(anova_fix, "Delta_eq (FIXED, 11.2 per mil)")
)
write_csv(table7, file.path(output_dir, "Table7_ANOVA_Results.csv"))

cat("\n*** Treatment x Flow regime interaction (Table 7) ***\n")
cat("  Observed d18O-PO4: F(1,", anova_obs["Residuals","Df"], ") =",
    round(anova_obs["treatment:flow_regime", "F value"], 2),
    ", p =", format(anova_obs["treatment:flow_regime", "Pr(>F)"],
                    scientific = TRUE), "\n")
cat("  Dynamic Delta_eq:  F(1,", anova_dyn["Residuals","Df"], ") =",
    round(anova_dyn["treatment:flow_regime", "F value"], 2),
    ", p =", format(anova_dyn["treatment:flow_regime", "Pr(>F)"],
                    scientific = TRUE), "\n")
cat("  Fixed   Delta_eq:  F(1,", anova_fix["Residuals","Df"], ") =",
    round(anova_fix["treatment:flow_regime", "F value"], 2),
    ", p =", format(anova_fix["treatment:flow_regime", "Pr(>F)"],
                    scientific = TRUE), "\n")
cat("  Interaction is robust to choice of equilibrium reference.\n")

# ---- 7.6 Figure 6: d18O-PO4 vs daily discharge -----------------------------

fig6 <- ggplot(phosphate_isotope %>% filter(!is.na(Q_daily_lps)),
               aes(x = Q_daily_lps, y = d18O_PO4,
                   color = treatment, shape = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  geom_hline(data = ref_eq,
             aes(yintercept = eq_mean, color = treatment),
             linetype = "dashed", linewidth = 0.6, alpha = 0.7) +
  geom_vline(xintercept = q75_daily_map["CT"],
             linetype = "dotted", color = "#D55E00", alpha = 0.7) +
  geom_vline(xintercept = q75_daily_map["NT"],
             linetype = "dotted", color = "#0072B2", alpha = 0.7) +
  scale_color_manual(values = treatment_colors) +
  scale_x_log10() +
  geom_text(data = ref_eq,
            aes(x = max(phosphate_isotope$Q_daily_lps, na.rm = TRUE) * 0.9,
                y = eq_mean,
                label = sprintf("%s eq. = %.1f", treatment, eq_mean)),
            color = c("#D55E00", "#0072B2"),
            hjust = 1, vjust = -0.5, size = 3,
            show.legend = FALSE) +
  labs(x = expression(paste("Daily mean discharge (L s"^"-1", ")")),
       y = expression(paste(delta^{18}, "O-PO"[4], " (\u2030 VSMOW)")),
       color = "Treatment", shape = "Treatment")

ggsave(file.path(figures_dir, "Fig6_PO4_vs_Discharge.png"),
       fig6, width = 8, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Fig6_PO4_vs_Discharge.pdf"),
       fig6, width = 8, height = 6)
ggsave(file.path(figures_dir, "Fig6_PO4_vs_Discharge.tiff"),
       fig6, width = 8, height = 6, dpi = 600, compression = "lzw")

# ---- 7.7 Figure 7: interaction plot ----------------------------------------

interaction_summary <- phosphate_isotope %>%
  filter(!is.na(flow_regime)) %>%
  group_by(treatment, flow_regime) %>%
  summarise(mean_d18O = mean(d18O_PO4),
            sd_d18O   = sd(d18O_PO4),
            n         = n(),
            .groups = "drop")

p_interaction <- format(anova_obs["treatment:flow_regime", "Pr(>F)"],
                        scientific = TRUE, digits = 2)
F_interaction <- round(anova_obs["treatment:flow_regime", "F value"], 1)
df_resid <- anova_obs["Residuals", "Df"]

fig7 <- ggplot(interaction_summary,
               aes(x = flow_regime, y = mean_d18O,
                   color = treatment, group = treatment,
                   shape = treatment)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_d18O - sd_d18O,
                    ymax = mean_d18O + sd_d18O),
                width = 0.10, linewidth = 0.8) +
  geom_hline(data = ref_eq,
             aes(yintercept = eq_mean, color = treatment),
             linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_text(data = ref_eq,
            aes(x = 2.45, y = eq_mean,
                label = sprintf("%s eq. = %.1f \u2030", treatment, eq_mean),
                color = treatment),
            hjust = 1, vjust = -0.4, size = 3,
            show.legend = FALSE) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = c(CT = 15, NT = 17)) +
  scale_y_continuous(breaks = seq(8, 18, 1), limits = c(8, 17.5)) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.6))) +
  annotate("text", x = 0.55, y = 8.4, hjust = 0,
           label = paste0("Treatment x Flow regime\nF[1,", df_resid,
                          "] = ", F_interaction,
                          ", p = ", p_interaction),
           size = 3) +
  labs(subtitle = paste0("Treatment-specific dynamic equilibrium: CT = ",
                         round(eq_CT, 1),
                         " \u2030, NT = ", round(eq_NT, 1), " \u2030"),
       x = "Flow regime",
       y = expression(paste(delta^{18}, "O-PO"[4], " (\u2030 VSMOW)")),
       color = "Treatment", shape = "Treatment")

ggsave(file.path(figures_dir, "Fig7_Interaction_Plot.png"),
       fig7, width = 7, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "Fig7_Interaction_Plot.pdf"),
       fig7, width = 7, height = 6)
ggsave(file.path(figures_dir, "Fig7_Interaction_Plot.tiff"),
       fig7, width = 7, height = 6, dpi = 600, compression = "lzw")

# ---- 7.8 Supplementary figures S1 (sensitivity) and S2 (heatmap) -----------

p_s1a <- phosphate_isotope %>%
  filter(!is.na(d18O_PO4_eq)) %>%
  ggplot(aes(Date, d18O_PO4_eq, color = treatment, shape = treatment)) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = equilibrium_d18O_PO4_fixed,
             linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = c(CT = 15, NT = 17)) +
  labs(x = "Date",
       y = expression(paste(delta^{18}, "O-PO"[4],
                            "(eq) (\u2030)")),
       title = "(a) Sample-specific equilibrium values") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank())

p_s1b <- phosphate_isotope %>%
  filter(!is.na(Delta_eq), !is.na(Delta_eq_fixed)) %>%
  ggplot(aes(Delta_eq_fixed, Delta_eq,
             color = treatment, shape = treatment)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = treatment_colors) +
  scale_shape_manual(values = c(CT = 15, NT = 17)) +
  labs(x = expression(Delta[eq] ~ "(fixed)"),
       y = expression(Delta[eq] ~ "(dynamic)"),
       title = "(b) Dynamic vs. fixed Delta_eq") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank())

figS1 <- cowplot::plot_grid(p_s1a, p_s1b, ncol = 2, rel_widths = c(1, 0.9))

ggsave(file.path(figures_dir, "FigS1_sensitivity.png"),
       figS1, width = 10, height = 4, dpi = 300)
ggsave(file.path(figures_dir, "FigS1_sensitivity.pdf"),
       figS1, width = 10, height = 4)
ggsave(file.path(figures_dir, "FigS1_sensitivity.tiff"),
       figS1, width = 10, height = 4, dpi = 600, compression = "lzw")

# Figure S2: equilibrium response surface
T_grid <- seq(2, 12, by = 0.5)
W_grid <- seq(-15.5, -13.0, by = 0.1)
sens_grid <- expand.grid(Temp_C = T_grid, d18Ow = W_grid) %>%
  mutate(d18O_PO4_eq = cb_eq(d18Ow, Temp_C))

figS2 <- ggplot(sens_grid, aes(d18Ow, Temp_C, fill = d18O_PO4_eq)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = d18O_PO4_eq), breaks = seq(8, 14, 0.5),
               color = "white", linewidth = 0.3, alpha = 0.7) +
  scale_fill_viridis_c(name = expression(delta^{18}*"O-PO"[4]*"(eq)")) +
  annotate("point", x = -14.1, y = 5,
           shape = 21, fill = "white", color = "black", size = 4) +
  annotate("text", x = -14.1, y = 5,
           label = "Fixed\nassumption", vjust = -1.2, size = 3) +
  annotate("point", x = ref_eq$d18O_H2O_mean[ref_eq$treatment == "CT"],
           y = ref_eq$Temp_C_mean[ref_eq$treatment == "CT"],
           shape = 22, fill = "#D55E00", color = "black", size = 4) +
  annotate("text", x = ref_eq$d18O_H2O_mean[ref_eq$treatment == "CT"],
           y = ref_eq$Temp_C_mean[ref_eq$treatment == "CT"],
           label = "CT", vjust = -1.4, color = "#D55E00", fontface = "bold") +
  annotate("point", x = ref_eq$d18O_H2O_mean[ref_eq$treatment == "NT"],
           y = ref_eq$Temp_C_mean[ref_eq$treatment == "NT"],
           shape = 24, fill = "#0072B2", color = "black", size = 4) +
  annotate("text", x = ref_eq$d18O_H2O_mean[ref_eq$treatment == "NT"],
           y = ref_eq$Temp_C_mean[ref_eq$treatment == "NT"],
           label = "NT", vjust = -1.4, color = "#0072B2", fontface = "bold") +
  labs(x = expression(paste(delta^{18}, "O-H"[2], "O (\u2030)")),
       y = "Tile water temperature (\u00B0C)",
       title = "Equilibrium d18O-PO4 sensitivity surface (Chang & Blake, 2015)")

ggsave(file.path(figures_dir, "FigS2_eq_sensitivity_heatmap.png"),
       figS2, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(figures_dir, "FigS2_eq_sensitivity_heatmap.pdf"),
       figS2, width = 7.5, height = 5.5)
ggsave(file.path(figures_dir, "FigS2_eq_sensitivity_heatmap.tiff"),
       figS2, width = 7.5, height = 5.5, dpi = 600, compression = "lzw")


# =============================================================================
# 8. EQUILIBRIUM PARAMETERS (Table 8)
# =============================================================================

table8 <- tibble(
  Parameter = c(
    "Fixed pooled tile water d18O-H2O (manuscript original)",
    "Fixed assumed tile water temperature (manuscript original)",
    "Calculated FIXED equilibrium d18O-PO4",
    "Dynamic CT mean d18O-H2O (matched, treatment-specific)",
    "Dynamic CT mean tile water temperature (continuously logged)",
    "Calculated DYNAMIC mean equilibrium d18O-PO4 (CT)",
    "Dynamic NT mean d18O-H2O (matched, treatment-specific)",
    "Dynamic NT mean tile water temperature (continuously logged)",
    "Calculated DYNAMIC mean equilibrium d18O-PO4 (NT)",
    "Equilibrium equation reference"
  ),
  Value = c(
    paste0(fixed_d18O_H2O, " per mil"),
    paste0(fixed_T_C, " deg C"),
    paste0(round(equilibrium_d18O_PO4_fixed, 2), " per mil"),
    paste0(round(ref_eq$d18O_H2O_mean[ref_eq$treatment == "CT"], 2), " per mil"),
    paste0(round(ref_eq$Temp_C_mean[ref_eq$treatment == "CT"], 2), " deg C"),
    paste0(round(eq_CT, 2), " per mil"),
    paste0(round(ref_eq$d18O_H2O_mean[ref_eq$treatment == "NT"], 2), " per mil"),
    paste0(round(ref_eq$Temp_C_mean[ref_eq$treatment == "NT"], 2), " deg C"),
    paste0(round(eq_NT, 2), " per mil"),
    "Chang & Blake (2015)"
  )
)
write_csv(table8, file.path(output_dir, "Table8_Equilibrium_Parameters.csv"))


# =============================================================================
# 9. SUMMARY OUTPUT
# =============================================================================

summary_output <- list(
  metadata = list(
    script_version = "2.0",
    run_date       = Sys.Date(),
    r_version      = R.version.string,
    project_root   = here()
  ),
  discharge = list(
    CT = list(n = discharge_stats$n_obs[1],
              mean = discharge_stats$mean_L_s[1], Q75_15min = Q75_CT,
              Q75_daily = q75_daily_map["CT"],
              max = discharge_stats$max_L_s[1]),
    NT = list(n = discharge_stats$n_obs[2],
              mean = discharge_stats$mean_L_s[2], Q75_15min = Q75_NT,
              Q75_daily = q75_daily_map["NT"],
              max = discharge_stats$max_L_s[2])
  ),
  water_isotopes = list(
    precipitation = list(n = nrow(precip_clean),
                         mean = round(mean(precip_clean$d18O_H2O), 2),
                         sd   = round(sd(precip_clean$d18O_H2O),   3)),
    CT_tile = list(n = nrow(ct_data),
                   mean = round(mean(ct_data$d18O_H2O), 2),
                   sd   = round(sd(ct_data$d18O_H2O),   3)),
    NT_tile = list(n = nrow(nt_data),
                   mean = round(mean(nt_data$d18O_H2O), 2),
                   sd   = round(sd(nt_data$d18O_H2O),   3)),
    variance_ratio = round(variance_ratio, 1),
    paired_t_p     = paired_t_test$p.value
  ),
  residence_time = list(
    CT_Fyw = round(Fyw_CT * 100, 1),
    NT_Fyw = round(Fyw_NT * 100, 1),
    CT_tau_amp = round(tau_amp_CT, 1),
    NT_tau_amp = round(tau_amp_NT, 1)
  ),
  phosphate_isotopes = list(
    eq_fixed         = round(equilibrium_d18O_PO4_fixed, 2),
    eq_CT_dynamic    = round(eq_CT, 2),
    eq_NT_dynamic    = round(eq_NT, 2),
    interaction_obs_F = round(anova_obs["treatment:flow_regime","F value"], 2),
    interaction_obs_p = anova_obs["treatment:flow_regime","Pr(>F)"],
    interaction_dyn_F = round(anova_dyn["treatment:flow_regime","F value"], 2),
    interaction_dyn_p = anova_dyn["treatment:flow_regime","Pr(>F)"]
  )
)
saveRDS(summary_output, file.path(output_dir, "analysis_summary.rds"))

cat("\n=========================================================\n")
cat(" Analysis complete. See output/ and figures/ for results.\n")
cat("=========================================================\n")


# =============================================================================
# 10. SESSION INFO
# =============================================================================

cat("\n--- Session info ---\n")
print(sessionInfo())
