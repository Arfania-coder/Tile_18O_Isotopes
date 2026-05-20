# =============================================================================
# Craig Diagram with Inset (Local Meteoric Water Line + Evaporation Line)
#
# Companion to 01_main_analysis.R. Produces Figure 3 (or equivalent panel,
# depending on final figure numbering in the manuscript) for the AGU JGR
# Biogeosciences submission.
#
# Authors:  Arfania, H., Kayler, Z. E., Strawn, D. G., Brooks, E. S., Laan, M.
# Site:     R. J. Cook Agronomy Farm (CAF) LTAR, Pullman, WA, USA
# Script:   02_craig_diagram.R
#
# Note: Per Z. Kayler's review, the evaporation line through the tile-water
# cluster is invisible at full extent because tile-drain d18O-H2O and d2H-H2O
# span only a narrow range. A magnified inset (LMWL trace + evap line, fitted
# only to tile data) is therefore included.
#
# Inputs (under data/):
#   18Ow_Prcep.csv, 18Ow_Tile.csv
#
# Outputs (under figures/):
#   Fig_Craig_Diagram_with_Inset.png (300 dpi), .pdf, .tiff (600 dpi, LZW)
# =============================================================================


# ---- 1. Setup ---------------------------------------------------------------

rm(list = ls())

required_packages <- c("here", "tidyverse", "lubridate", "cowplot")
new_pkgs <- required_packages[!(required_packages %in%
                                  installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_packages, library, character.only = TRUE))

data_dir    <- here("data")
figures_dir <- here("figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

treatment_colors <- c("CT" = "#D55E00", "NT" = "#0072B2")

parse_mixed_date <- function(x) {
  pdt <- parse_date_time(as.character(x),
                         orders = c("mdY", "Ymd", "Ymd HMS", "Ymd HM",
                                    "mdY HMS", "mdY HM", "dmY"),
                         quiet = TRUE)
  as_date(pdt)
}


# ---- 2. Load water isotope data --------------------------------------------

cat("\n=== Loading isotope data ===\n")

precip_d18Ow <- read_csv(file.path(data_dir, "18Ow_Prcep.csv"),
                         show_col_types = FALSE) %>%
  rename(Date = Date, d18O_H2O = `18O_H2O`, d2H_H2O = Detrium_H2O) %>%
  mutate(Date = parse_mixed_date(Date)) %>%
  filter(!is.na(d18O_H2O), !is.na(d2H_H2O), !is.na(Date))

cat("Precipitation samples: n =", nrow(precip_d18Ow), "\n")

tile_h2o_raw <- read_csv(file.path(data_dir, "18Ow_Tile.csv"),
                         show_col_types = FALSE, name_repair = "minimal")
names(tile_h2o_raw)[1] <- "Date"
tile_h2o <- tile_h2o_raw %>%
  select(Date, NT_18OW, NT_Deutrium, CT_18OW, CT_Deutrium) %>%
  filter(!is.na(Date)) %>%
  mutate(Date = parse_mixed_date(Date)) %>%
  filter(!is.na(Date))

tile_dual <- tile_h2o %>%
  transmute(Date,
            CT_d18O = CT_18OW, CT_d2H = CT_Deutrium,
            NT_d18O = NT_18OW, NT_d2H = NT_Deutrium) %>%
  pivot_longer(-Date, names_to = c("treatment", ".value"),
               names_sep = "_") %>%
  filter(!is.na(d18O), !is.na(d2H))

cat("Tile-drain samples: CT n =",
    sum(tile_dual$treatment == "CT"),
    "| NT n =", sum(tile_dual$treatment == "NT"), "\n")


# ---- 3. Identify precipitation outliers (d-excess > 3 SD from mean) ---------

precip_full <- precip_d18Ow %>%
  mutate(d_excess = d2H_H2O - 8 * d18O_H2O)

dex_mean <- mean(precip_full$d_excess, na.rm = TRUE)
dex_sd   <- sd(precip_full$d_excess,   na.rm = TRUE)

precip_full <- precip_full %>%
  mutate(is_outlier = abs(d_excess - dex_mean) > 3 * dex_sd)

precip_clean   <- precip_full %>% filter(!is_outlier)
precip_outlier <- precip_full %>% filter(is_outlier)

cat("d-excess: mean =", round(dex_mean, 2), "per mil, SD =",
    round(dex_sd, 2), "per mil\n")
cat("Outliers excluded: n =", nrow(precip_outlier), "\n")


# ---- 4. Fit LMWL and tile-water evaporation line ---------------------------

lmwl_fit <- lm(d2H_H2O ~ d18O_H2O, data = precip_clean)
lmwl_a   <- coef(lmwl_fit)[2]
lmwl_b   <- coef(lmwl_fit)[1]
lmwl_r2  <- summary(lmwl_fit)$r.squared

evap_fit <- lm(d2H ~ d18O, data = tile_dual)
evap_a   <- coef(evap_fit)[2]
evap_b   <- coef(evap_fit)[1]
evap_r2  <- summary(evap_fit)$r.squared

cat(sprintf("LMWL:      dD = %.2f * d18O + %.2f  (R^2 = %.3f, n = %d)\n",
            lmwl_a, lmwl_b, lmwl_r2, nrow(precip_clean)))
cat(sprintf("Evap line: dD = %.2f * d18O + %.2f  (R^2 = %.3f, n = %d)\n",
            evap_a, evap_b, evap_r2, nrow(tile_dual)))


# ---- 5. Build main panel ----------------------------------------------------

tile_xrange <- range(tile_dual$d18O)
tile_yrange <- range(tile_dual$d2H)
box_xpad <- 0.3
box_ypad <- 2.5

seg_x_lo <- tile_xrange[1] - 0.4
seg_x_hi <- tile_xrange[2] + 0.4

all_pts <- bind_rows(
  precip_clean %>% transmute(
    d18O = d18O_H2O, d2H = d2H_H2O,
    group = sprintf("Precipitation (n = %d)", nrow(precip_clean))),
  precip_outlier %>% transmute(
    d18O = d18O_H2O, d2H = d2H_H2O,
    group = "Precip. outlier (excluded)"),
  tile_dual %>% filter(treatment == "CT") %>% transmute(
    d18O, d2H,
    group = sprintf("Tile - CT (n = %d)",
                    sum(tile_dual$treatment == "CT"))),
  tile_dual %>% filter(treatment == "NT") %>% transmute(
    d18O, d2H,
    group = sprintf("Tile - NT (n = %d)",
                    sum(tile_dual$treatment == "NT")))
)

group_levels <- c(
  sprintf("Precipitation (n = %d)", nrow(precip_clean)),
  "Precip. outlier (excluded)",
  sprintf("Tile - CT (n = %d)", sum(tile_dual$treatment == "CT")),
  sprintf("Tile - NT (n = %d)", sum(tile_dual$treatment == "NT"))
)
all_pts$group <- factor(all_pts$group, levels = group_levels)

shape_vals  <- c(21, 4, 22, 24)
fill_vals   <- c("gray70", NA, "#D55E00", "#0072B2")
color_vals  <- c("black", "black", "black", "black")
size_vals   <- c(2.8, 3.5, 2.6, 2.6)
stroke_vals <- c(0.3, 1.2, 0.3, 0.3)
names(shape_vals) <- names(fill_vals) <- names(color_vals) <-
  names(size_vals) <- names(stroke_vals) <- group_levels

rect_xmin <- tile_xrange[1] - box_xpad
rect_xmax <- tile_xrange[2] + box_xpad
rect_ymin <- tile_yrange[1] - box_ypad
rect_ymax <- tile_yrange[2] + box_ypad

craig_main <- ggplot() +
  geom_abline(slope = lmwl_a, intercept = lmwl_b,
              color = "#3B7DB8", linewidth = 0.9) +
  geom_segment(aes(x = seg_x_lo, xend = seg_x_hi,
                   y = evap_a * seg_x_lo + evap_b,
                   yend = evap_a * seg_x_hi + evap_b),
               color = "#C0392B", linetype = "dashed", linewidth = 0.9) +
  geom_point(data = all_pts,
             aes(x = d18O, y = d2H,
                 shape = group, fill = group, color = group,
                 size = group, stroke = group)) +
  scale_shape_manual(values = shape_vals,  name = NULL,
                     limits = group_levels) +
  scale_fill_manual(values = fill_vals,   name = NULL,
                    limits = group_levels, na.value = "transparent") +
  scale_color_manual(values = color_vals,  name = NULL,
                     limits = group_levels) +
  scale_size_manual(values = size_vals,    name = NULL,
                    limits = group_levels) +
  scale_discrete_manual("stroke", values = stroke_vals,
                        name = NULL, limits = group_levels) +
  annotate("rect",
           xmin = rect_xmin, xmax = rect_xmax,
           ymin = rect_ymin, ymax = rect_ymax,
           fill = NA, color = "black",
           linewidth = 0.5, linetype = "dotted") +
  labs(x = expression(paste(delta^{18}, "O-H"[2], "O (\u2030 VSMOW)")),
       y = expression(paste(delta, "D-H"[2], "O (\u2030 VSMOW)"))) +
  theme_bw(base_size = 13) +
  theme(
    text             = element_text(color = "black"),
    axis.text        = element_text(color = "black", size = 13),
    axis.title       = element_text(color = "black", size = 14),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    legend.text      = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.key.size  = unit(0.45, "cm"),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )


# ---- 6. Build magnified inset ----------------------------------------------

inset_xpad <- 0.25
inset_ypad <- 1.5
xr <- tile_xrange + c(-inset_xpad, inset_xpad)
yr <- tile_yrange + c(-inset_ypad, inset_ypad)

craig_inset <- ggplot() +
  geom_abline(slope = lmwl_a, intercept = lmwl_b,
              color = "#3B7DB8", linewidth = 0.8, alpha = 0.7) +
  geom_abline(slope = evap_a, intercept = evap_b,
              color = "#C0392B", linetype = "dashed", linewidth = 1.2) +
  geom_point(data = tile_dual %>% filter(treatment == "CT"),
             aes(d18O, d2H),
             shape = 22, fill = "#D55E00", color = "black",
             size = 3.2, stroke = 0.3) +
  geom_point(data = tile_dual %>% filter(treatment == "NT"),
             aes(d18O, d2H),
             shape = 24, fill = "#0072B2", color = "black",
             size = 3.2, stroke = 0.3) +
  coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
  annotate("text",
           x = xr[1] + 0.05 * diff(xr),
           y = yr[1] + 0.12 * diff(yr),
           label = sprintf("Evap. line slope = %.2f", evap_a),
           hjust = 0, size = 3.4,
           color = "#C0392B", fontface = "bold") +
  annotate("text",
           x = xr[2] - 0.05 * diff(xr),
           y = yr[2] - 0.08 * diff(yr),
           label = sprintf("LMWL slope = %.2f", lmwl_a),
           hjust = 1, size = 3.4,
           color = "#3B7DB8", fontface = "bold") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    text            = element_text(color = "black"),
    axis.text       = element_text(color = "black", size = 10),
    axis.ticks      = element_line(color = "black", linewidth = 0.35),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = "black",
                                   linewidth = 0.4),
    plot.margin     = margin(2, 2, 2, 2)
  )


# ---- 7. Compose final figure with connector lines --------------------------

inset_x      <- 0.12
inset_y      <- 0.52
inset_width  <- 0.36
inset_height <- 0.40

panel_x0 <- 0.075; panel_x1 <- 0.74
panel_y0 <- 0.095; panel_y1 <- 0.965

x_lim <- c(min(c(precip_clean$d18O_H2O, tile_dual$d18O,
                 precip_outlier$d18O_H2O), na.rm = TRUE) - 0.5,
           max(c(precip_clean$d18O_H2O, tile_dual$d18O,
                 precip_outlier$d18O_H2O), na.rm = TRUE) + 0.5)
y_lim <- c(min(c(precip_clean$d2H_H2O, tile_dual$d2H,
                 precip_outlier$d2H_H2O), na.rm = TRUE) - 5,
           max(c(precip_clean$d2H_H2O, tile_dual$d2H,
                 precip_outlier$d2H_H2O), na.rm = TRUE) + 5)

data_to_canvas_x <- function(x) {
  panel_x0 + (x - x_lim[1]) / diff(x_lim) * (panel_x1 - panel_x0)
}
data_to_canvas_y <- function(y) {
  panel_y0 + (y - y_lim[1]) / diff(y_lim) * (panel_y1 - panel_y0)
}

rect_left_cx   <- data_to_canvas_x(rect_xmin)
rect_right_cx  <- data_to_canvas_x(rect_xmax)
rect_top_cy    <- data_to_canvas_y(rect_ymax)

inset_right_x  <- inset_x + inset_width
inset_bottom_y <- inset_y

fig_craig <- cowplot::ggdraw(craig_main) +
  cowplot::draw_line(x = c(inset_right_x, rect_left_cx),
                     y = c(inset_bottom_y, rect_top_cy),
                     color = "gray40", size = 0.4, linetype = "dotted") +
  cowplot::draw_line(x = c(inset_x, rect_right_cx),
                     y = c(inset_bottom_y, rect_top_cy),
                     color = "gray40", size = 0.4, linetype = "dotted") +
  cowplot::draw_plot(craig_inset,
                     x = inset_x, y = inset_y,
                     width = inset_width, height = inset_height)

ggsave(file.path(figures_dir, "Fig_Craig_Diagram_with_Inset.png"),
       fig_craig, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(figures_dir, "Fig_Craig_Diagram_with_Inset.pdf"),
       fig_craig, width = 10, height = 6.5)
ggsave(file.path(figures_dir, "Fig_Craig_Diagram_with_Inset.tiff"),
       fig_craig, width = 10, height = 6.5, dpi = 600, compression = "lzw")

cat("\nCraig diagram with inset saved to figures/.\n")


# ---- 8. Caption-ready summary ---------------------------------------------

cat("\n=== Caption-ready statistics ===\n")
cat(sprintf("LMWL:      dD = %.2f * d18O + %.2f   (R^2 = %.2f, n = %d)\n",
            lmwl_a, lmwl_b, lmwl_r2, nrow(precip_clean)))
cat(sprintf("Evap line: dD = %.2f * d18O - %.2f   (R^2 = %.2f, n = %d)\n",
            evap_a, abs(evap_b), evap_r2, nrow(tile_dual)))
cat(sprintf("Tile d18O range: %.2f to %.2f per mil\n",
            tile_xrange[1], tile_xrange[2]))
cat(sprintf("Tile d2H range:  %.2f to %.2f per mil\n",
            tile_yrange[1], tile_yrange[2]))
if (nrow(precip_outlier) > 0) {
  cat("Excluded precipitation outlier(s):\n")
  print(precip_outlier %>% select(Date, d18O_H2O, d2H_H2O, d_excess))
}
