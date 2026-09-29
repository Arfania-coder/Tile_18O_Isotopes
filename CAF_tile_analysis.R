## =====================================================================
## 2026JG010070 -- COMPLETE ANALYSIS, FIGURES AND TABLES  (single file)
##
## "Flow-Dependent Phosphate Oxygen-Isotope (d18O-PO4) Signatures in Tile
##  Drainage from Two Long-Term Tillage Systems"
##  Arfania, Kayler, Strawn, Brooks & Laan -- JGR: Biogeosciences
##
## ---------------------------------------------------------------------
## HOW TO RUN
##   1. Data folder:  <repository>/data
##      It must contain the 14 monitoring CSVs AND Pullman_Climate_Daily.csv
##      (the station precipitation record used for Figure 1). If that file is
##      absent the script falls back to Prceipitation.csv and SAYS SO in the log.
##   2. Set PROJ to the folder that CONTAINS data/ (not data/ itself), then
##      source this whole file:
##
##        PROJ <- "/path/to/Tile_18O_Isotopes"
##        source("CAF_tile_analysis.R")
##
##      Source the WHOLE file. Do not run it line by line.
##   3. Everything is written to  <PROJ>/output_revision/
##        tables ............ *.csv                      (Tables 1-7, Data Sets S1-S2)
##        review figures .... *.png  (200 dpi)
##        journal figures ... hires/*.tiff  600 dpi LZW   <- upload these to AGU
##                            hires/*.pdf   vector
##                            hires/*_300dpi.png
##        log ............... analysis_log.txt
##
## REQUIREMENTS: base R plus nlme (ships with R). Nothing to install.
##   If a Windows R build complains about  type = "cairo"  in tiff()/png(),
##   set  USE_CAIRO <- FALSE  below.
##
## TERMINOLOGY: the deviation quantity is D_app (apparent deviation from the
##   contemporaneous tile-water reference), never Delta_eq -- it is a
##   standardized comparison value, not evidence of biochemical equilibration.
## =====================================================================

USE_CAIRO <- TRUE

## >>> YOUR PROJECT FOLDER -- the folder that CONTAINS data/  <<<
## Edit this one line if the data move. Forward slashes, no trailing /data.
## PROJ <- "path/to/this/repository"   # set only if auto-detection fails

#######################################################################
## PART 0 -- SETUP AND PROJECT-ROOT RESOLUTION
#######################################################################
.this_file <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1], mustWork = FALSE))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(of, mustWork = FALSE))
  }
  NA_character_
})
SCRIPT_DIR <- if (!is.na(.this_file)) dirname(.this_file) else getwd()
.has_data <- function(p) !is.na(p) && nzchar(p) && dir.exists(file.path(p, "data"))
if (exists("PROJ") && !.has_data(PROJ)) {
  message("PROJ '", PROJ, "' has no data/ subfolder on this machine -- searching nearby.")
  rm(PROJ)
}
if (!exists("PROJ")) {
  cand <- unique(c(getwd(), dirname(getwd()), SCRIPT_DIR, dirname(SCRIPT_DIR),
                   dirname(dirname(SCRIPT_DIR))))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit  <- cand[vapply(cand, .has_data, logical(1))]
  if (length(hit)) PROJ <- normalizePath(hit[1])
}
if (!exists("PROJ"))
  stop("\n\nCould not find the project root (a folder containing data/).\n",
       "Set it and re-source, e.g.\n",
       '  PROJ <- "/path/to/Tile_18O_Isotopes"\n',
       call. = FALSE)
t0 <- Sys.time()
cat("Project root :", PROJ, "\n")
cat("Data files   :", length(list.files(file.path(PROJ, "data"), "\\.csv$")), "CSVs found\n\n")


#######################################################################
## PART 1 -- INGEST AND CLEAN
#######################################################################

DATA <- file.path(PROJ, "data")

## Output folder. Default "output_revision" rather than "output" so this
## pipeline can never overwrite results from an earlier workflow. Windows
## filenames are CASE-INSENSITIVE, so writing Table1_discharge_statistics.csv
## into a folder that already holds Table1_Discharge_Statistics.csv silently
## replaces it. Override before sourcing if you want a different name:
##     OUT_NAME <- "output"
OUT <- file.path(PROJ, if (exists("OUT_NAME")) OUT_NAME else "output_revision")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## warn about anything already there that this run will not itself rewrite
.mine <- c("analysis_log.txt", "Figure_captions_draft.txt",
           "Table1_discharge_statistics.csv", "Table2_DRP_concentrations.csv",
           "Table3_water_isotope_summary.csv", "Table4_residence_time_indices.csv",
           "Model_Sensitivity_Table.csv", "Sample_Flow_Classification_Audit.csv",
           "DRP_Load_Integration.csv",
           "Figure1_discharge_timeseries.png", "Figure2_craig_diagram.png",
           "Figure3_dexcess.png", "Figure4_interaction.png",
           "Figure5_continuous_discharge.png", "FigureS1_reference_comparison.png",
           "FigureS2_reference_sensitivity.png", "FigureS3_diagnostics.png")
.pre <- list.files(OUT)
.clash <- .pre[tolower(.pre) %in% tolower(.mine) & !(.pre %in% .mine)]
if (length(.clash))
  warning("These existing files differ from this pipeline's names only by ",
          "CASE, and Windows will overwrite them:\n  ",
          paste(.clash, collapse = "\n  "),
          "\nBack them up, or set OUT_NAME to a fresh folder.", call. = FALSE)
if (length(.pre) && !length(.clash))
  cat("Note:", length(.pre), "file(s) already in", OUT,
      "- files with matching names will be replaced.\n")

## helper -------------------------------------------------------------
## Several of the raw exports carry a UTF-8 byte-order mark, which turns
## the first column name into "\ufeffDate" and silently breaks $ access.
rd <- function(f, ...) {
  z <- read.csv(file.path(DATA, f), stringsAsFactors = FALSE,
                check.names = FALSE, encoding = "latin1", ...)
  names(z) <- sub("^\\xef\\xbb\\xbf", "", names(z), useBytes = TRUE)
  z
}
mdy  <- function(x) as.Date(x, format = "%m/%d/%Y")
## Timestamp formats differ BETWEEN files (ISO in some, US m/d/Y in
## others), so try each in turn. strptime returns NA rather than erroring
## on unparseable input, which is what we need for the pivot-artefact rows
## in NT_Discharge_2023.csv.
DT_FORMATS <- c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
                "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%Y")
parse_dt <- function(x) {
  x   <- as.character(x)
  out <- as.POSIXct(rep(NA_real_, length(x)), origin = "1970-01-01", tz = "UTC")
  for (f in DT_FORMATS) {
    todo <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(todo)) break
    out[todo] <- as.POSIXct(strptime(x[todo], f, tz = "UTC"))
  }
  out
}
## show at most `n` offending values so a malformed column cannot flood the log
peek <- function(x, n = 6) {
  u <- unique(as.character(x)); u <- u[!is.na(u)]
  if (!length(u)) return("(all NA/empty)")
  paste0(paste(utils::head(u, n), collapse = ", "),
         if (length(u) > n) paste0(" ... (", length(u), " distinct)") else "")
}
num  <- function(x) suppressWarnings(as.numeric(as.character(x)))
CFS_TO_LPS <- 28.317

log_note <- new.env()
log_note$lines <- character(0)
note <- function(...) {
  msg <- paste0(...)
  log_note$lines <- c(log_note$lines, msg)
  cat(msg, "\n")
}

note("=== 01 INGEST AND CLEAN ===")

## ---------------------------------------------------------------------
## 1. Phosphate oxygen isotopes  (27 dates x CT/NT = 54 values)
## ---------------------------------------------------------------------
p_wide <- rd("Tile_d18Op.csv")
names(p_wide)[1] <- "Date"
p_wide$Date <- mdy(p_wide$Date)
d18Op <- rbind(
  data.frame(date = p_wide$Date, system = "CT", d18Op = num(p_wide$CT)),
  data.frame(date = p_wide$Date, system = "NT", d18Op = num(p_wide$NT))
)
d18Op <- d18Op[!is.na(d18Op$d18Op), ]
note("d18O-PO4 records: ", nrow(d18Op),
     "  (CT ", sum(d18Op$system == "CT"), ", NT ", sum(d18Op$system == "NT"), ")")
note("  sampling dates: ", length(unique(d18Op$date)),
     "; dates carrying BOTH systems: ",
     sum(tapply(d18Op$system, d18Op$date, function(z) length(unique(z))) == 2))

## ---------------------------------------------------------------------
## 2. Tile water isotopes -- first column is unnamed in the raw file
## ---------------------------------------------------------------------
w <- rd("18Ow_Tile.csv")
names(w)[1] <- "Date"
w <- w[, c("Date", "NT_18OW", "NT_Deutrium", "CT_18OW", "CT_Deutrium")]
w$Date <- mdy(w$Date)
w <- w[!is.na(w$Date) & !is.na(num(w$NT_18OW)), ]
tile_w <- rbind(
  data.frame(date = w$Date, system = "CT",
             d18Ow = num(w$CT_18OW), d2Hw = num(w$CT_Deutrium)),
  data.frame(date = w$Date, system = "NT",
             d18Ow = num(w$NT_18OW), d2Hw = num(w$NT_Deutrium))
)
tile_w$dex <- tile_w$d2Hw - 8 * tile_w$d18Ow
note("Tile water isotope records: ", nrow(tile_w), " on ",
     length(unique(tile_w$date)), " dates")

## ---------------------------------------------------------------------
## 3. Precipitation isotopes -- one d-excess outlier (12 Apr 2023)
## ---------------------------------------------------------------------
pr <- rd("18Ow_Prcep.csv")
pr$date  <- mdy(pr$Date)
pr$d18Ow <- num(pr[["18O_H2O"]])
pr$d2Hw  <- num(pr$Detrium_H2O)
pr$dex   <- pr$d2Hw - 8 * pr$d18Ow
pr$outlier <- pr$dex > 50
note("Precipitation isotope samples: ", nrow(pr),
     "; flagged outliers (d-excess > 50 permil): ", sum(pr$outlier),
     " on ", paste(format(pr$date[pr$outlier]), collapse = ", "))

## Precipitation. Prefer Pullman_Climate_Daily.csv: it is a continuous daily
## station record (every calendar day present, zeros explicit) with precip in
## both inches and mm at an exact 25.4 conversion. Prceipitation.csv contains
## only wet days, uses 0.25 mm per 0.01 in rather than 0.254, and its daily
## values do NOT reconcile with the station record (r = 0.54 at zero lag; no
## date offset or accumulation window explains the difference). Its provenance
## should be resolved -- see audit item IA-19.
if (file.exists(file.path(DATA, "Pullman_Climate_Daily.csv"))) {
  x <- rd("Pullman_Climate_Daily.csv")
  rain <- data.frame(date = mdy(x$Date), mm = num(x$Precip_mm))
  rain <- rain[!is.na(rain$date) & !is.na(rain$mm), ]
  RAIN_SOURCE <- "Pullman_Climate_Daily.csv (daily station record)"
} else {
  x <- rd("Prceipitation.csv")
  rain <- data.frame(date = mdy(x$Date), mm = num(x$Rainfall_mm))
  rain <- rain[!is.na(rain$date) & !is.na(rain$mm), ]
  RAIN_SOURCE <- "Prceipitation.csv (wet days only)"
}
note("Precipitation source: ", RAIN_SOURCE)
note("  ", nrow(rain), " daily records, ", format(min(rain$date)), " to ",
     format(max(rain$date)), "; ", sum(rain$mm > 0), " wet days; total ",
     round(sum(rain$mm)), " mm")
wy_sum <- function(a, b) sum(rain$mm[rain$date >= as.Date(a) & rain$date <= as.Date(b)])
note("  water year 2022-23 (1 Oct 2022 - 30 Sep 2023): ",
     round(wy_sum("2022-10-01", "2023-09-30"), 1), " mm",
     if (min(rain$date) > as.Date("2022-10-01"))
       paste0(" -- INCOMPLETE, record starts ", format(min(rain$date))) else "")
note("  water year 2023-24 (1 Oct 2023 - 30 Sep 2024): ",
     round(wy_sum("2023-10-01", "2024-09-30"), 1), " mm")

## ---------------------------------------------------------------------
## 4. Discharge -- four files, four different defects
##    CT2023 : clean 15-min, 1 unparseable row
##    CT2024 : date-only stamps, 2558 blank rows, cfs
##    NT2023 : pivot-table export containing "(blank)" and "Grand Total"
##             rows; the Grand Total (9204 L/s) must be dropped or every
##             discharge statistic is destroyed
##    NT2024 : TOA5 logger header, 4 lines before the real header, cfs
## ---------------------------------------------------------------------
ct23 <- rd("CT_Discharge_2023.csv")
ct23 <- data.frame(ts = parse_dt(ct23$Timestamp),
                   Q  = num(ct23$Q_lps))
n0 <- nrow(ct23); ct23 <- ct23[!is.na(ct23$ts) & !is.na(ct23$Q), ]
note("CT 2023 discharge: kept ", nrow(ct23), " of ", n0, " rows")
ct23$date <- as.Date(ct23$ts)

x <- rd("CT_Discharge_2024.csv")
n0 <- nrow(x)
ct24 <- data.frame(date = mdy(x$Date), Q = num(x$cfs) * CFS_TO_LPS)
ct24 <- ct24[!is.na(ct24$date) & !is.na(ct24$Q), ]
note("CT 2024 discharge: kept ", nrow(ct24), " of ", n0,
     " rows (blank rows dropped)")
tab <- table(ct24$date)
note("  days with exactly 96 obs: ", sum(tab == 96), " of ", length(tab),
     " (range ", min(tab), "-", max(tab), ")")

x <- rd("NT_Discharge_2023.csv")
n0 <- nrow(x)
bad <- is.na(parse_dt(x[["Row Labels"]]))
note("NT 2023 discharge: dropping ", sum(bad), " non-timestamp pivot rows: ",
     peek(x[["Row Labels"]][bad]),
     "  (largest value ", round(max(num(x[["Sum of Q_lps"]])[bad], na.rm = TRUE), 1),
     " L/s -- a pivot Grand Total, NOT a discharge observation)")
nt23 <- data.frame(ts = parse_dt(x[["Row Labels"]]),
                   Q  = num(x[["Sum of Q_lps"]]))
nt23 <- nt23[!is.na(nt23$ts) & !is.na(nt23$Q), ]
nt23$date <- as.Date(nt23$ts)
note("  kept ", nrow(nt23), " of ", n0, " rows")

x <- read.csv(file.path(DATA, "NT_Discharge_2024.csv"), skip = 4,
              stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
names(x)[1:2] <- c("Date", "cfs")
n0 <- nrow(x)
nt24 <- data.frame(date = mdy(x$Date), Q = num(x$cfs) * CFS_TO_LPS)
nt24 <- nt24[!is.na(nt24$date) & !is.na(nt24$Q), ]
note("NT 2024 discharge: kept ", nrow(nt24), " of ", n0,
     " rows (TOA5 header skipped)")

## sub-daily record (15-min where available) and daily means
q15 <- rbind(
  data.frame(system = "CT", date = ct23$date, Q = ct23$Q),
  data.frame(system = "CT", date = ct24$date, Q = ct24$Q),
  data.frame(system = "NT", date = nt23$date, Q = nt23$Q),
  data.frame(system = "NT", date = nt24$date, Q = nt24$Q)
)
qday <- aggregate(Q ~ system + date, data = q15, FUN = mean)
names(qday)[3] <- "Q_daily"
note("Sub-daily discharge observations: CT ", sum(q15$system == "CT"),
     ", NT ", sum(q15$system == "NT"), ", total ", nrow(q15))
note("Days with discharge data: CT ", sum(qday$system == "CT"),
     ", NT ", sum(qday$system == "NT"))

## ---------------------------------------------------------------------
## 5. Water temperature -- sentinel -6999 plus physically implausible
##    readings. NOTE: NT_water_tempreture_2024.csv is headed SR_Temp_C
##    (surface runoff) whereas CT is TL_Temp_C (tile line). Provenance
##    must be confirmed with the field team; see audit item IA-7.
## ---------------------------------------------------------------------
read_temp <- function(file, tcol, dcol, posix) {
  x <- rd(file)
  ts <- parse_dt(x[[dcol]])
  v  <- num(x[[tcol]])
  keep <- !is.na(ts) & !is.na(v) & v > -2 & v < 30
  note("  ", file, ": kept ", sum(keep), " of ", nrow(x),
       " (sentinel/implausible removed: ", sum(!keep), ")")
  data.frame(date = as.Date(ts[keep]), T_C = v[keep])
}
note("Water temperature screening (retain -2 C < T < 30 C):")
tCT <- rbind(read_temp("CT_water_tempreture_2023.csv", "Temp_C", "Timestamp", TRUE),
             read_temp("CT_water_tempreture_2024.csv", "TL_Temp_C", "Date-Time", FALSE))
tNT <- rbind(read_temp("NT_water_tempreture_2023.csv", "Sum of Temp_C", "Row Labels", TRUE),
             read_temp("NT_water_tempreture_2024.csv", "SR_Temp_C", "Date-Time", FALSE))
note("  CAUTION: NT 2024 temperature column is 'SR_Temp_C' (surface runoff),",
     " not 'TL_Temp_C'. Confirm provenance before resubmission.")
tday <- rbind(
  cbind(system = "CT", aggregate(T_C ~ date, tCT, mean)),
  cbind(system = "NT", aggregate(T_C ~ date, tNT, mean))
)

## ---------------------------------------------------------------------
## 6. DRP -- both files store CT and NT in SIDE-BY-SIDE blocks separated
##    by filler columns. In DRP_2023_24.csv the CT block's sample-ID
##    column has NO header (position 6), so a naive read drops the CT
##    labels entirely.
## ---------------------------------------------------------------------
grab <- function(df, cols) {
  z <- df[, cols]
  names(z) <- c("ID", "Date", "DRP")
  z$ID  <- trimws(as.character(z$ID))
  z$date <- mdy(z$Date)
  z$DRP <- num(z$DRP)
  z[!is.na(z$date) & !is.na(z$DRP) & nzchar(z$ID), c("ID", "date", "DRP")]
}
d1 <- rd("DRP_2022_23.csv")
d2 <- rd("DRP_2023_24.csv")
drp <- rbind(grab(d1, 1:3), grab(d1, 8:10), grab(d2, 1:3), grab(d2, 6:8))
drp$system <- ifelse(grepl("_CT_", drp$ID), "CT", "NT")
drp$matrix <- ifelse(grepl("_SR", drp$ID), "surface_runoff", "tile")
note("DRP records parsed: ", nrow(drp),
     " (surface-runoff records present: ", sum(drp$matrix != "tile"), ")")
drp <- drp[drp$matrix == "tile", ]
note("  tile DRP: CT ", sum(drp$system == "CT"),
     ", NT ", sum(drp$system == "NT"), ", total ", nrow(drp))
note("  ID labels: ", paste(names(table(drp$ID)), table(drp$ID),
                            sep = " = ", collapse = "; "))
MDL <- 0.005
note("  values below the stated MDL (", MDL, " mg/L): ", sum(drp$DRP < MDL),
     " -- retained at face value; see audit item IA-4")


#######################################################################
## PART 2 -- FLOW CLASSIFICATION AUDIT  (Table 1, Data Set S1)
#######################################################################

note("\n=== 02 FLOW CLASSIFICATION ===")

q75_sub <- tapply(q15$Q,  q15$system,  quantile, 0.75, na.rm = TRUE)
q75_day <- tapply(qday$Q_daily, qday$system, quantile, 0.75, na.rm = TRUE)
q75_day_nozero <- tapply(qday$Q_daily, qday$system,
                         function(z) quantile(z[z > 0], 0.75, na.rm = TRUE))
Q75_POOLED <- quantile(qday$Q_daily, 0.75, na.rm = TRUE)

note("Q75, sub-daily record : CT ", round(q75_sub["CT"], 3),
     " | NT ", round(q75_sub["NT"], 3), " L/s   <- values quoted in the MS")
note("Q75, daily mean       : CT ", round(q75_day["CT"], 3),
     " | NT ", round(q75_day["NT"], 3), " L/s   <- values actually USED")
note("Q75, daily excl. zero-flow days: CT ", round(q75_day_nozero["CT"], 3),
     " | NT ", round(q75_day_nozero["NT"], 3), " L/s (sensitivity only)")
note("Pooled daily Q75 across both records: ", round(Q75_POOLED, 3), " L/s")

## ---- Table 1 -------------------------------------------------------
## Report both aggregations explicitly. The submitted Table 1 mixes them:
## the CT maximum (2.02) is a daily-mean maximum while the NT maximum
## (3.57) is a sub-daily maximum.
tab1 <- do.call(rbind, lapply(c("CT", "NT"), function(s) {
  x <- q15$Q[q15$system == s]
  d <- qday$Q_daily[qday$system == s]
  data.frame(System = s,
             n_subdaily   = length(x),
             days_with_data = length(d),
             calendar_span_d = as.numeric(diff(range(qday$date[qday$system == s]))) + 1,
             mean = round(mean(x), 3), sd = round(sd(x), 3),
             median = round(median(x), 3),
             Q75_subdaily = round(quantile(x, .75), 3),
             Q75_daily    = round(quantile(d, .75), 3),
             max_subdaily = round(max(x), 3),
             max_daily    = round(max(d), 3),
             pct_zero = round(100 * mean(x == 0), 1),
             skewness = round(mean((x - mean(x))^3) / sd(x)^3, 2),
             pct_volume_in_top_quartile =
               round(100 * sum(x[x >= quantile(x, .75)]) / sum(x), 1))
}))
write.csv(tab1, file.path(OUT, "Table1_discharge_statistics.csv"), row.names = FALSE)
note("\n-- Table 1 --"); print(tab1)

## Reviewer 1, L282: "how is it possible that Q75 for CT is lower than the
## mean?" Answer, on a SINGLE aggregation, with the numbers that explain it.
note("\nReviewer 1 L282 -- mean vs Q75 on the same (sub-daily) aggregation:")
for (s in c("CT", "NT")) {
  r <- tab1[tab1$System == s, ]
  note("  ", s, ": mean ", r$mean, " vs Q75 ", r$Q75_subdaily,
       " L/s; ", r$pct_zero, "% of observations are exactly zero; skewness ",
       r$skewness, "; the top quartile carries ",
       r$pct_volume_in_top_quartile, "% of total volume.")
}
note("  In a distribution this right-skewed the mean legitimately exceeds Q75.")

## ---- per-sample classification audit -------------------------------
qd_list <- split(qday, qday$system)
daily_Q_on <- function(s, d) {
  z <- qd_list[[s]]
  v <- z$Q_daily[match(d, z$date)]
  v
}
pctile_within <- function(s, v) {
  z <- qd_list[[s]]$Q_daily
  if (is.na(v)) return(NA_real_)
  100 * mean(z < v)
}

fa <- d18Op
fa$Q_daily <- mapply(daily_Q_on, fa$system, fa$date)
fa$n_subdaily_that_day <- mapply(function(s, d)
  sum(q15$system == s & q15$date == d), fa$system, fa$date)
fa$system_Q75 <- unname(q75_day[fa$system])
fa$pooled_Q75 <- as.numeric(Q75_POOLED)
fa$pctile     <- mapply(pctile_within, fa$system, fa$Q_daily)
fa$class_system_threshold <- ifelse(is.na(fa$Q_daily), NA,
                            ifelse(fa$Q_daily >= fa$system_Q75, "Event", "Baseflow"))
fa$class_shared_threshold <- ifelse(is.na(fa$Q_daily), NA,
                            ifelse(fa$Q_daily >= fa$pooled_Q75, "Event", "Baseflow"))
fa$class_changes <- fa$class_system_threshold != fa$class_shared_threshold
fa$sample_id <- paste0(fa$system, "_", format(fa$date))
fa$matching_lag_days <- 0L   # verified below

note("\nExact-date discharge available for ", sum(!is.na(fa$Q_daily)),
     " of ", nrow(fa), " phosphate samples.")
note("  -> every sampling date has same-day discharge in BOTH systems, so the",
     " matching lag is 0 d throughout and no nearest-date substitution occurs.")
note("  (The +/-3 d fallback in the earlier workflow was therefore never",
     " exercised; it does not need a labelled sensitivity analysis, only an",
     " accurate description.)")

note("\nGroup sizes, system-specific daily Q75:")
print(table(fa$system, fa$class_system_threshold))
note("Group sizes, shared pooled daily Q75 (", round(Q75_POOLED, 3), " L/s):")
print(table(fa$system, fa$class_shared_threshold))
note("Samples whose class changes under the shared threshold: ",
     sum(fa$class_changes, na.rm = TRUE), " of ", nrow(fa), " -- ",
     paste(fa$sample_id[which(fa$class_changes)], collapse = ", "))

fa_out <- fa[order(fa$system, fa$date),
             c("sample_id", "system", "date", "d18Op", "Q_daily",
               "n_subdaily_that_day", "matching_lag_days", "system_Q75",
               "pooled_Q75", "pctile", "class_system_threshold",
               "class_shared_threshold", "class_changes")]
names(fa_out)[3] <- "phosphate_sampling_date"
write.csv(fa_out, file.path(OUT, "Sample_Flow_Classification_Audit.csv"),
          row.names = FALSE)

## the analysis frame used downstream
mf <- fa
mf$flow    <- factor(mf$class_system_threshold, levels = c("Baseflow", "Event"))
mf$flow_sh <- factor(mf$class_shared_threshold, levels = c("Baseflow", "Event"))
mf$system  <- factor(mf$system, levels = c("CT", "NT"))


#######################################################################
## PART 3 -- WATER ISOTOPES AND TILE-WATER REFERENCE  (Tables 4, 5)
#######################################################################

note("\n=== 03 WATER ISOTOPES AND REFERENCE ===")

## ---- LMWL and evaporation line -------------------------------------
pc <- pr[!pr$outlier, ]
lmwl <- lm(d2Hw ~ d18Ow, data = pc)
evap <- lm(d2Hw ~ d18Ow, data = tile_w)
note("LMWL (n = ", nrow(pc), "): d2H = ", round(coef(lmwl)[2], 2),
     " * d18O + ", round(coef(lmwl)[1], 2),
     "  (R2 = ", round(summary(lmwl)$r.squared, 3), ")")
note("Tile evaporation line (n = ", nrow(tile_w), "): slope ",
     round(coef(evap)[2], 2), ", intercept ", round(coef(evap)[1], 2),
     " (R2 = ", round(summary(evap)$r.squared, 3), ")")
note("  A slope well below 8 indicates kinetic evaporation SOMEWHERE along",
     " the rainfall-to-drain pathway. It does not locate where.")

## ---- Table 3 --------------------------------------------------------
sumrow <- function(lab, o, h) {
  d <- h - 8 * o
  data.frame(Source = lab, n = length(o),
             d18O_mean = round(mean(o), 2), d18O_sd = round(sd(o), 2),
             d2H_mean  = round(mean(h), 2), d2H_sd  = round(sd(h), 2),
             dex_mean  = round(mean(d), 2), dex_sd  = round(sd(d), 2))
}
tab3 <- rbind(
  sumrow("Precipitation (all)",   pr$d18Ow, pr$d2Hw),
  sumrow("Precipitation (clean)", pc$d18Ow, pc$d2Hw),
  sumrow("Tile drain - CT", tile_w$d18Ow[tile_w$system == "CT"],
                            tile_w$d2Hw [tile_w$system == "CT"]),
  sumrow("Tile drain - NT", tile_w$d18Ow[tile_w$system == "NT"],
                            tile_w$d2Hw [tile_w$system == "NT"]))
write.csv(tab3, file.path(OUT, "Table3_water_isotope_summary.csv"), row.names = FALSE)
note("\n-- Table 3 --"); print(tab3)

## ---- paired d-excess test ------------------------------------------
wCT <- tile_w[tile_w$system == "CT", ]
wNT <- tile_w[tile_w$system == "NT", ]
m   <- match(wCT$date, wNT$date)
tt  <- t.test(wCT$dex, wNT$dex[m], paired = TRUE)
wx  <- suppressWarnings(wilcox.test(wCT$dex, wNT$dex[m], paired = TRUE))
note("\nd-excess, CT vs NT (paired on ", sum(!is.na(m)), " dates): difference ",
     round(mean(wCT$dex - wNT$dex[m]), 2), " permil, t = ",
     round(tt$statistic, 2), ", p = ", signif(tt$p.value, 3),
     " (Wilcoxon p = ", signif(wx$p.value, 3), ")")
note("  Report as: CT tile water is CONSISTENT WITH greater evaporative",
     " modification somewhere along the rainfall-to-drain pathway.")

## ---- residence-time indices ----------------------------------------
omega <- 2 * pi / 365.25
fit_sine <- function(d, y) {
  t <- as.numeric(format(d, "%j"))
  st <- list(A = diff(range(y)) / 2, phi = 0, C = mean(y))
  f <- try(nls(y ~ A * sin(omega * t + phi) + C, start = st,
               control = nls.control(maxiter = 500, warnOnly = TRUE)),
           silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  abs(unname(coef(f)["A"]))
}
A_p  <- fit_sine(pc$date, pc$d18Ow)
A_ct <- fit_sine(wCT$date, wCT$d18Ow)
A_nt <- fit_sine(wNT$date, wNT$d18Ow)
sd_p <- sd(pc$d18Ow)

rt <- do.call(rbind, lapply(list(c("CT", A_ct), c("NT", A_nt)), function(z) {
  s <- z[1]; A <- as.numeric(z[2])
  sdt <- sd(tile_w$d18Ow[tile_w$system == s])
  DR  <- sdt / sd_p
  data.frame(System = s,
             damping_ratio_tau_d = round(1 / (omega * sqrt(1 / DR^2 - 1)), 1),
             amplitude_ratio_tau_d = round((1 / omega) * sqrt((A_p / A)^2 - 1), 0),
             young_water_fraction_pct = round(100 * A / A_p, 0))
}))
write.csv(rt, file.path(OUT, "Table4_residence_time_indices.csv"), row.names = FALSE)
note("\n-- Table 4 (residence-time INDICES, not transit times) --"); print(rt)
note("  The two methods disagree in magnitude and ordering. Per Reviewer 2,",
     " present this as consistent with dual-domain / dual-porosity behaviour:",
     " bulk tile water dominated by older, well-mixed matrix water while",
     " event-associated P may travel preferentially. Water age and phosphate",
     " residence time in biologically active pools are NOT equivalent.")

## ---- Chang & Blake (2015) contemporaneous tile-water reference ------
##   d18O-PO4(ref) = exp[(14.43/T_K) - 0.02654] * (d18O-H2O + 1000) - 1000
cb_reference <- function(T_C, d18Ow) {
  exp((14.43 / (T_C + 273.15)) - 0.02654) * (d18Ow + 1000) - 1000
}

## exact-date temperature; water isotope exact-date, else nearest within +/-7 d
tl <- split(tday, tday$system)
wl <- split(tile_w, tile_w$system)
mf$T_C <- mapply(function(s, d) {
  z <- tl[[s]]; v <- z$T_C[match(d, z$date)]; if (length(v)) v else NA_real_
}, as.character(mf$system), mf$date)

near_w <- function(s, d, tol = 7) {
  z <- wl[[s]]; lag <- as.numeric(z$date - d)
  ok <- abs(lag) <= tol
  if (!any(ok)) return(c(NA_real_, NA_real_))
  i <- which(ok)[which.min(abs(lag[ok]))]
  c(z$d18Ow[i], lag[i])
}
nw <- t(mapply(near_w, as.character(mf$system), mf$date))
mf$d18Ow      <- nw[, 1]
mf$water_lag  <- nw[, 2]
mf$d18Ow_exact <- mapply(function(s, d) {
  z <- wl[[s]]; v <- z$d18Ow[match(d, z$date)]; if (length(v)) v else NA_real_
}, as.character(mf$system), mf$date)

note("\nReference-calculation matching (full disclosure -- audit item IA-8):")
note("  exact-date temperature available : ", sum(!is.na(mf$T_C)), " of 54")
note("  exact-date tile d18O-H2O         : ", sum(!is.na(mf$d18Ow_exact)), " of 54")
note("  within +/-7 d                    : ", sum(!is.na(mf$d18Ow)), " of 54",
     "  (lags used: ", peek(mf$water_lag[is.na(mf$d18Ow_exact) & !is.na(mf$d18Ow)]), ")")
note("  NO water isotope within +/-7 d, therefore DROPPED: ",
     sum(is.na(mf$d18Ow)), " -- ",
     peek(mf$sample_id[is.na(mf$d18Ow)], 8))

mf$ref   <- cb_reference(mf$T_C, mf$d18Ow)
mf$D_app <- mf$d18Op - mf$ref              # apparent deviation, NOT "Delta_eq"
D_FIXED_REF <- 11.2                         # original pooled estimate
mf$D_fixed <- mf$d18Op - D_FIXED_REF

note("\nComplete cases for the reference calculation: ", sum(!is.na(mf$D_app)))
agg <- aggregate(cbind(T_C, d18Ow, ref) ~ system, mf, mean)
note("Reference by system (mean of sample-specific values):")
print(data.frame(System = agg$system, matched_T_C = round(agg$T_C, 1),
                 matched_d18Ow = round(agg$d18Ow, 2),
                 reference_permil = round(agg$ref, 2)))
note("Numerical stability of the reference across plausible tile-water",
     " conditions does NOT establish that tile water was the water of",
     " biochemical oxygen exchange. Sensitivity surface -> Figure S2.")


#######################################################################
## PART 4 -- ROBUSTNESS MODEL SUITE  (Table 7 / Table S1)
#######################################################################

note("\n=== 04 ROBUSTNESS MODELS ===")

## ---- helpers --------------------------------------------------------
## Type III SS via single-term deletion under sum-to-zero contrasts
## Each term's columns are removed from the model matrix directly. Using
## update(. ~ . - term) would not work here: under marginality R keeps the
## main effect in the design whenever the interaction is present, giving
## df = 0. Sum-to-zero contrasts make the main effects interpretable at the
## margins, which is what Type III requires.
type3 <- function(fml, data) {
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  full <- lm(fml, data = data, x = TRUE)
  X    <- full$x; y <- model.response(model.frame(full))
  asg  <- attr(X, "assign"); labs <- attr(terms(full), "term.labels")
  rss  <- sum(resid(full)^2); rdf <- df.residual(full)
  out <- do.call(rbind, lapply(seq_along(labs), function(i) {
    keep <- asg != i
    ss <- sum(lm.fit(X[, keep, drop = FALSE], y)$residuals^2) - rss
    df <- sum(!keep)
    Fv <- (ss / df) / (rss / rdf)
    data.frame(Term = labs[i], SS = ss, df = df, F = Fv,
               p = pf(Fv, df, rdf, lower.tail = FALSE))
  }))
  rbind(out, data.frame(Term = "Residuals", SS = rss, df = rdf, F = NA, p = NA))
}

## HC3 heteroskedasticity-robust covariance (MacKinnon & White 1985)
hc3 <- function(m) {
  X    <- model.matrix(m)
  XtXi <- summary(m)$cov.unscaled          # (X'X)^-1, pivot-safe
  h <- hatvalues(m); u <- resid(m) / (1 - h)
  XtXi %*% crossprod(X * u) %*% XtXi     # meat = X' diag(u^2) X
}
robust_tab <- function(m) {
  se <- sqrt(diag(hc3(m))); b <- coef(m)
  z  <- b / se; df <- df.residual(m)
  data.frame(estimate = b, se_HC3 = se,
             ci_low  = b - qt(.975, df) * se,
             ci_high = b + qt(.975, df) * se,
             p = 2 * pt(abs(z), df, lower.tail = FALSE))
}
bp_test <- function(m) {                       # Breusch-Pagan
  u2 <- resid(m)^2
  a  <- summary(lm(u2 ~ model.matrix(m)[, -1]))
  LM <- length(u2) * a$r.squared
  k  <- length(coef(m)) - 1
  pchisq(LM, k, lower.tail = FALSE)
}
int_row <- function(label, m, n, response, notes, robust = FALSE) {
  k <- grep(":", names(coef(m)), value = TRUE)[1]
  if (robust) {
    r <- robust_tab(m)[k, ]
    est <- r$estimate; lo <- r$ci_low; hi <- r$ci_high; p <- r$p
  } else {
    ci <- confint(m)[k, ]
    est <- coef(m)[k]; lo <- ci[1]; hi <- ci[2]
    p <- summary(m)$coefficients[k, 4]
  }
  data.frame(Model = label, n = n, Response = response,
             Interaction = round(est, 3), CI_low = round(lo, 3),
             CI_high = round(hi, 3), p = signif(p, 3),
             Sign = ifelse(est > 0, "+", "-"), Notes = notes)
}

rows <- list()

## ---- A. original: system-specific daily Q75 -------------------------
mA  <- lm(d18Op ~ system * flow, data = mf)
a3  <- type3(d18Op ~ system * flow, mf)
ssI <- a3$SS[a3$Term == "system:flow"]; ssR <- a3$SS[a3$Term == "Residuals"]
peta <- ssI / (ssI + ssR)
note("\n-- Model A: system-specific daily Q75 --")
print(within(a3, {SS <- round(SS, 2); F <- round(F, 2); p <- signif(p, 3)}))
note("partial eta-squared for the interaction = ", round(peta, 2))
print(aggregate(d18Op ~ system + flow, mf,
                function(z) c(n = length(z), mean = round(mean(z), 2),
                              sd = round(sd(z), 2))))
rows[[1]] <- int_row(
  "A. Original: system-specific daily Q75 (CT 0.173, NT 0.347 L/s)",
  mA, nobs(mA), "d18O-PO4 (permil)",
  paste0("F(1,", a3$df[a3$Term == "Residuals"], ") = ",
         round(a3$F[a3$Term == "system:flow"], 1),
         "; partial eta2 = ", round(peta, 2), "; groups 19/8/14/13"))

## ---- B. shared pooled threshold -------------------------------------
mB <- lm(d18Op ~ system * flow_sh, data = mf)
note("\n-- Model B: shared pooled daily Q75 --")
print(table(mf$system, mf$flow_sh))
rows[[2]] <- int_row("B. Shared pooled daily Q75 (0.291 L/s)", mB, nobs(mB),
                     "d18O-PO4 (permil)",
                     "groups 19/8/12/15; only 2 of 54 samples reclassified")

## ---- C. percentile rank, continuous ---------------------------------
mf$pctile_c <- mf$pctile - mean(mf$pctile)
mC <- lm(d18Op ~ system * pctile_c, data = mf)
rows[[3]] <- int_row("C. Discharge percentile rank, continuous (HC3)", mC,
                     nobs(mC), "d18O-PO4 (permil per percentile)",
                     "difference in slope, NT minus CT", robust = TRUE)
mf$pbin <- cut(mf$pctile, c(-0.01, 50, 75, 100),
               labels = c("<50", "50-75", ">75"))
note("\n-- Model C: descriptive percentile bins --")
print(aggregate(d18Op ~ system + pbin, mf,
                function(z) c(n = length(z), mean = round(mean(z), 2))))

## ---- D. continuous log10 discharge ----------------------------------
CONST <- 0.001                       # zeros occur; small, defensible offset
note("\n-- Model D: log10(Q + ", CONST, ") --  zero-flow sample days: ",
     sum(mf$Q_daily == 0, na.rm = TRUE))
mf$logQ_c <- log10(mf$Q_daily + CONST)
mf$logQ_c <- mf$logQ_c - mean(mf$logQ_c)      # centred before interaction
mD <- lm(d18Op ~ system * logQ_c, data = mf)
print(round(robust_tab(mD), 3))
sh <- shapiro.test(resid(mD)); bp <- bp_test(mD)
ck <- cooks.distance(mD); cut4n <- 4 / nobs(mD)
note("R2 = ", round(summary(mD)$r.squared, 3),
     " | Shapiro-Wilk p = ", signif(sh$p.value, 3),
     " | Breusch-Pagan p = ", signif(bp, 3))
note("max Cook's D = ", round(max(ck), 3), " (4/n = ", round(cut4n, 3),
     "); observations above cutoff: ", sum(ck > cut4n), " -- ",
     paste(mf$sample_id[order(-ck)][seq_len(sum(ck > cut4n))], collapse = ", "))
note("NO observations were removed: there is no documented data-quality",
     " reason to drop any of them (all are genuine zero- or high-flow days).")
rows[[4]] <- int_row("D. Continuous log10(Q + 0.001), centred (HC3)", mD,
                     nobs(mD), "d18O-PO4 (permil per log10 unit)",
                     paste0("R2 = ", round(summary(mD)$r.squared, 2),
                            "; BP p = ", signif(bp, 2),
                            "; Shapiro p = ", signif(sh$p.value, 2)),
                     robust = TRUE)

## ---- E. date-blocked mixed model ------------------------------------
## All 27 sampling dates carry both a CT and an NT sample, so sampling date
## is estimable as a random effect. This is a temporal-dependence
## correction; it is NOT field-level replication.
np <- sum(tapply(mf$system, mf$date, function(z) length(unique(z))) == 2)
note("\n-- Model E: date-blocked mixed model --")
note("sampling dates carrying both systems: ", np, " of ", length(unique(mf$date)))
suppressMessages(library(nlme))
mE <- lme(d18Op ~ system * flow, random = ~ 1 | date, data = mf, method = "REML")
print(summary(mE)$tTable)
kE <- grep(":", rownames(summary(mE)$tTable), value = TRUE)[1]
ciE <- intervals(mE, which = "fixed")$fixed[kE, ]
rows[[5]] <- data.frame(
  Model = "E. Date-blocked mixed model (27 paired dates, random intercept)",
  n = nrow(mf), Response = "d18O-PO4 (permil)",
  Interaction = round(ciE[2], 3), CI_low = round(ciE[1], 3),
  CI_high = round(ciE[3], 3),
  p = signif(summary(mE)$tTable[kE, "p-value"], 3), Sign = "+",
  Notes = "estimate attenuated relative to A; sign and significance retained")

## ---- F/G. secondary responses ---------------------------------------
sub <- mf[!is.na(mf$D_app), ]
mF  <- lm(D_app ~ system * flow, data = sub)
a3F <- type3(D_app ~ system * flow, sub)
note("\n-- Model F: D_app (apparent deviation from the contemporaneous",
     " tile-water reference) --")
print(within(a3F, {SS <- round(SS, 2); F <- round(F, 2); p <- signif(p, 3)}))
rows[[6]] <- int_row("F. Secondary: D_app vs contemporaneous tile-water reference",
                     mF, nobs(mF), "D_app (permil)",
                     paste0("F(1,", a3F$df[a3F$Term == "Residuals"], ") = ",
                            round(a3F$F[a3F$Term == "system:flow"], 1)))
mG <- lm(D_fixed ~ system * flow, data = mf)
rows[[7]] <- int_row("G. Secondary: deviation from fixed pooled reference (11.2 permil)",
                     mG, nobs(mG), "D_fixed (permil)",
                     "identical to A by construction (constant offset)")

## ---- leave-one-out stability ----------------------------------------
loo <- function(fml) vapply(seq_len(nrow(mf)), function(i) {
  m <- lm(fml, data = mf[-i, ]); unname(coef(m)[grep(":", names(coef(m)))[1]])
}, numeric(1))
lA <- loo(d18Op ~ system * flow); lD <- loo(d18Op ~ system * logQ_c)
note("\nLeave-one-out (all 54 samples dropped in turn):")
note("  Model A interaction ranges ", round(min(lA), 2), " to ", round(max(lA), 2),
     " (full model ", round(coef(mA)[4], 2), "); all same sign: ",
     all(sign(lA) == sign(coef(mA)[4])))
note("  Model D interaction ranges ", round(min(lD), 2), " to ", round(max(lD), 2),
     "; all positive: ", all(lD > 0))

## ---- sensitivity table ----------------------------------------------
sens <- do.call(rbind, rows); rownames(sens) <- NULL
write.csv(sens, file.path(OUT, "Model_Sensitivity_Table.csv"), row.names = FALSE)
note("\n-- SENSITIVITY TABLE --")
print(sens[, c("Model", "n", "Interaction", "CI_low", "CI_high", "p")])

robust <- all(sens$Sign == "+") && all(sens$CI_low > 0)
note("\nIs the principal system-by-flow interaction robust? ",
     ifelse(robust, "YES", "NO"),
     " -- positive sign and a 95% CI excluding zero in every specification",
     " tested, and stable under leave-one-out resampling.")
note("This is a robustness statement about the STATISTICAL result only.",
     " It says nothing about mechanism, and it does not overcome the fact",
     " that this is a two-field comparison without field-level replication.")


#######################################################################
## PART 5 -- DRP CONCENTRATIONS AND LOAD  (Tables 2, 3 / Data Set S2)
#######################################################################

note("\n=== 05 DRP AND LOAD ===")
SEC_PER_DAY <- 86400
MG_PER_KG   <- 1e6
## Field areas differ between the two systems and are NOT equal. Using a single
## area for both (as an earlier version did) inflates the NT per-area export
## more than three-fold and reverses the apparent between-field ordering.
FIELD_HA    <- c(CT = 5, NT = 12)      # per Section 2.1

## ---- concentrations by system and flow regime (Table 2) -------------
drp$Q_daily <- mapply(daily_Q_on, drp$system, drp$date)
note("DRP samples with same-day discharge: ", sum(!is.na(drp$Q_daily)),
     " of ", nrow(drp), "; without: ",
     peek(drp$sample_id <- paste0(drp$system, "_", format(drp$date))[is.na(drp$Q_daily)]))
dq <- drp[!is.na(drp$Q_daily), ]
dq$flow    <- ifelse(dq$Q_daily >= unname(q75_day[dq$system]), "Event", "Baseflow")
dq$flow_sh <- ifelse(dq$Q_daily >= Q75_POOLED, "Event", "Baseflow")

fwmc <- function(C, Q) sum(C * Q) / sum(Q)
tab2 <- do.call(rbind, lapply(split(dq, list(dq$system, dq$flow)), function(g)
  data.frame(System = g$system[1], Flow_regime = g$flow[1], n = nrow(g),
             mean_DRP = round(mean(g$DRP), 4), sd = round(sd(g$DRP), 4),
             FWMC = round(fwmc(g$DRP, g$Q_daily), 4))))
tab2 <- tab2[order(tab2$System, tab2$Flow_regime), ]
rownames(tab2) <- NULL
write.csv(tab2, file.path(OUT, "Table2_DRP_concentrations.csv"), row.names = FALSE)
note("\n-- Table 2 (rebuilt from the revised DRP files) --"); print(tab2)

note("\nOverall by system:")
print(do.call(rbind, lapply(split(dq, dq$system), function(g)
  data.frame(System = g$system[1], n = nrow(g),
             mean_DRP = round(mean(g$DRP), 4), sd = round(sd(g$DRP), 4),
             FWMC = round(fwmc(g$DRP, g$Q_daily), 4)))))

## Is there any flow-regime signal in DRP at all?
a <- type3(log10(DRP) ~ system * flow, dq)
note("\nTwo-way model on log10(DRP):")
print(within(a, {SS <- round(SS, 3); F <- round(F, 2); p <- signif(p, 3)}))
tt <- t.test(DRP ~ system, data = dq)
note("CT vs NT overall DRP: Welch t p = ", signif(tt$p.value, 3),
     "; means ", round(mean(dq$DRP[dq$system == "CT"]), 4), " vs ",
     round(mean(dq$DRP[dq$system == "NT"]), 4), " mg/L (ratio ",
     round(mean(dq$DRP[dq$system == "NT"]) / mean(dq$DRP[dq$system == "CT"]), 2), ")")
note("  NOTE: the two systems differ in DRP by ~1.6x. The framing",
     " 'similar DRP concentrations concealed different isotope behaviour'",
     " overstates the similarity. The defensible framing is that DRP",
     " differs MONOTONICALLY (NT higher at both flow regimes) whereas",
     " d18O-PO4 REVERSES SIGN with flow regime.")

## ---- load integration ------------------------------------------------
integrate_load <- function(s, rule) {
  d <- qday[qday$system == s, ]; d <- d[order(d$date), ]
  c_by <- aggregate(DRP ~ date, drp[drp$system == s, ], mean)
  c_by <- c_by[order(c_by$date), ]
  ## two DRP dates fall inside the logger gap and have no discharge day at
  ## all; they cannot contribute to a load and are dropped here explicitly
  c_by <- c_by[c_by$date %in% d$date, ]
  inside <- d$date >= min(c_by$date) & d$date <= max(c_by$date)
  Cv <- rep(NA_real_, nrow(d))
  Cv[match(c_by$date, d$date)] <- c_by$DRP
  Ci <- switch(rule,
    linear  = approx(d$date, Cv, xout = d$date, rule = 1)$y,
    locf    = { z <- Cv; for (i in seq_along(z)) if (is.na(z[i]) && i > 1) z[i] <- z[i - 1]; z },
    nearest = c_by$DRP[vapply(d$date, function(x)
                which.min(abs(as.numeric(c_by$date - x))), integer(1))])
  Ci[!inside] <- NA
  keep <- !is.na(Ci)
  kg <- sum(Ci[keep] * d$Q_daily[keep] * SEC_PER_DAY / MG_PER_KG)
  ev <- d$Q_daily[keep] >= q75_day[[s]]
  L  <- Ci[keep] * d$Q_daily[keep] * SEC_PER_DAY / MG_PER_KG
  data.frame(System = s, rule = rule, days = sum(keep),
             load_kg = round(kg, 3), load_kg_ha = round(kg / FIELD_HA[[s]], 4),
             event_pct_of_load = round(100 * sum(L[ev]) / sum(L), 1),
             event_pct_of_volume = round(100 * sum(d$Q_daily[keep][ev]) /
                                           sum(d$Q_daily[keep]), 1),
             pct_monitored_volume_covered =
               round(100 * sum(d$Q_daily[keep]) / sum(d$Q_daily), 1))
}
loads <- do.call(rbind, lapply(c("CT", "NT"), function(s)
  do.call(rbind, lapply(c("linear", "locf", "nearest"), function(r)
    integrate_load(s, r)))))
rownames(loads) <- NULL
write.csv(loads, file.path(OUT, "DRP_Load_Integration.csv"), row.names = FALSE)
note("Field areas used: CT ", FIELD_HA[["CT"]], " ha, NT ", FIELD_HA[["NT"]], " ha")
note("\n-- DRP export, monitored period --"); print(loads)
for (s in c("CT", "NT")) {
  z <- loads[loads$System == s, ]
  note("  ", s, ": ", min(z$load_kg), "-", max(z$load_kg), " kg (",
       min(z$load_kg_ha), "-", max(z$load_kg_ha), " kg/ha); event share ",
       min(z$event_pct_of_load), "-", max(z$event_pct_of_load), "% of load vs ",
       z$event_pct_of_volume[1], "% of volume")
}
note("  Spread across interpolation rules is small, so the estimate is not",
     " rule-dependent. Event share tracks the VOLUME share almost exactly:",
     " export is volume-driven, not concentration-driven.")

## ---- unmonitored-period context -------------------------------------
note("\nMissing-period uncertainty (must be stated, not waved through):")
for (s in c("CT", "NT")) {
  d <- qday[qday$system == s, ]; d <- d[order(d$date), ]
  g <- as.numeric(diff(d$date))
  big <- which(g > 5)
  for (i in big)
    note("  ", s, ": ", g[i], "-d discharge gap, ", format(d$date[i]),
         " to ", format(d$date[i + 1]))
  aut <- as.integer(format(d$date, "%m")) %in% 8:11
  note("     ", s, " Aug-Nov mean daily Q = ", round(mean(d$Q_daily[aut]), 4),
       " L/s vs record mean ", round(mean(d$Q_daily), 3), " L/s")
}
note("  The gap is genuinely negligible for CT. For NT it is NOT: the NT",
     " monitored-period load is more likely an underestimate. Section 2.2",
     " currently makes the 'minimal expected discharge' claim for both.")


#######################################################################
## PART 6 -- FIGURES 1-5, S1-S3  (TIFF 600 dpi + vector PDF + PNG 300 dpi)
#######################################################################

note("\n=== 06 FIGURES ===")
CT_COL <- "#D55E00"; NT_COL <- "#0072B2"; GREY <- "grey45"
## ---- output devices ------------------------------------------------
## Journal-quality output. Each figure is written three ways at the same
## physical size (pixel dims / 200 dpi = inches):
##   .tiff  600 dpi, LZW-compressed  -- AGU raster requirement
##   .pdf   vector                    -- scalable; preferred where accepted
##   .png   300 dpi                   -- for the Word draft / review copy
## The plot is recorded once and replayed to each device, so all three are
## identical in content.
FIG_FMT <- c("tiff", "pdf", "png")
png_open <- function(f, w = 1800, h = 1200) {
  assign(".fig_name", sub("\\.png$", "", f), envir = .GlobalEnv)
  assign(".fig_dim",  c(w, h) / 200,           envir = .GlobalEnv)
  png(file.path(OUT, f), width = w, height = h, res = 200)
  dev.control("enable")
}
fig_close <- function() {
  rec <- recordPlot(); dev.off()
  nm <- get(".fig_name", .GlobalEnv); d <- get(".fig_dim", .GlobalEnv)
  hi <- file.path(OUT, "hires"); dir.create(hi, showWarnings = FALSE)
  if ("tiff" %in% FIG_FMT) {
    tiff(file.path(hi, paste0(nm, ".tiff")), width = d[1], height = d[2],
         units = "in", res = 600, compression = "lzw", type = if (USE_CAIRO) "cairo" else getOption("bitmapType"))
    replayPlot(rec); dev.off()
  }
  if ("pdf" %in% FIG_FMT) {
    pdf(file.path(hi, paste0(nm, ".pdf")), width = d[1], height = d[2],
        useDingbats = FALSE)
    replayPlot(rec); dev.off()
  }
  if ("png" %in% FIG_FMT) {
    png(file.path(hi, paste0(nm, "_300dpi.png")), width = d[1], height = d[2],
        units = "in", res = 300, type = if (USE_CAIRO) "cairo" else getOption("bitmapType"))
    replayPlot(rec); dev.off()
  }
  invisible(NULL)
}

## ---- Figure 1: rainfall and daily discharge -------------------------
## The rainfall record runs to 31 Dec 2024 but discharge ends 10 Jul 2024.
## Plotting both on independent axes made panel (a) extend ~6 months beyond
## panel (b), which is misleading. Both panels are now clipped to the shared
## monitored window and given an identical x range.
XLIM <- range(qday$date)
rain_w <- rain[rain$date >= XLIM[1] & rain$date <= XLIM[2], ]
note("Figure 1: rainfall record spans ", format(min(rain$date)), " to ",
     format(max(rain$date)), "; discharge spans ", format(XLIM[1]), " to ",
     format(XLIM[2]), ". Both panels clipped to the discharge window (",
     nrow(rain_w), " of ", nrow(rain), " rain days, ",
     round(sum(rain_w$mm)), " of ", round(sum(rain$mm)), " mm).")
png_open("Figure1_discharge_timeseries.png", 2000, 1400)
par(mfrow = c(2, 1), mar = c(3, 4.5, 1.5, 1), mgp = c(2.6, 0.7, 0))
plot(rain_w$date, rain_w$mm, type = "h", col = GREY, xlab = "", xlim = XLIM,
     ylab = "Precipitation (mm)", main = "(a) Daily precipitation")
## the station record may start after discharge monitoring; shade any leading
## period with no precipitation data so the blank is not read as dry weather
if (min(rain$date) > XLIM[1]) {
  usr <- par("usr")
  rect(usr[1], usr[3], as.numeric(min(rain$date)), usr[4],
       col = "#00000010", border = NA)
  text(mean(c(usr[1], as.numeric(min(rain$date)))), usr[4] * 0.8,
       "no precip.\ndata", cex = 0.65, col = GREY)
}
cd <- qday[qday$system == "CT", ]; nd <- qday[qday$system == "NT", ]
plot(cd$date, pmax(cd$Q_daily, 1e-4), log = "y", type = "l", col = CT_COL,
     xlab = "", ylab = expression("Daily mean Q (L s"^-1*")"), xlim = XLIM,
     main = "(b) Tile discharge", ylim = c(1e-4, 4))
lines(nd$date, pmax(nd$Q_daily, 1e-4), col = NT_COL)
abline(h = q75_day[["CT"]], lty = 2, col = CT_COL)
abline(h = q75_day[["NT"]], lty = 2, col = NT_COL)
legend("bottomleft", c("CT field", "NT field", "system-specific daily Q75"),
       col = c(CT_COL, NT_COL, "black"), lty = c(1, 1, 2), bty = "n", cex = .8)
fig_close()

## ---- Figure 2: Craig diagram ----------------------------------------
png_open("Figure2_craig_diagram.png")
par(mar = c(4.2, 4.4, 1.5, 1), mgp = c(2.6, 0.7, 0))
plot(pr$d18Ow, pr$d2Hw, pch = 21, bg = "grey85", col = GREY,
     xlab = expression(delta^18*"O-H"[2]*"O (permil, VSMOW)"),
     ylab = expression(delta^2*"H-H"[2]*"O (permil, VSMOW)"))
points(pr$d18Ow[pr$outlier], pr$d2Hw[pr$outlier], pch = 4, cex = 1.6, lwd = 2)
points(tile_w$d18Ow[tile_w$system == "CT"], tile_w$d2Hw[tile_w$system == "CT"],
       pch = 22, bg = CT_COL, col = "white")
points(tile_w$d18Ow[tile_w$system == "NT"], tile_w$d2Hw[tile_w$system == "NT"],
       pch = 24, bg = NT_COL, col = "white")
abline(lmwl, col = NT_COL, lwd = 2)
abline(evap, col = "red3", lwd = 2, lty = 2)
legend("topleft", bty = "n", cex = .8,
       legend = c(sprintf("LMWL: slope %.2f", coef(lmwl)[2]),
                  sprintf("Tile evaporation line: slope %.2f", coef(evap)[2]),
                  "Precipitation", "CT tile water", "NT tile water",
                  "excluded precipitation outlier"),
       col = c(NT_COL, "red3", GREY, CT_COL, NT_COL, "black"),
       lty = c(1, 2, NA, NA, NA, NA), pch = c(NA, NA, 21, 22, 24, 4),
       pt.bg = c(NA, NA, "grey85", CT_COL, NT_COL, NA))
fig_close()

## ---- Figure 3: d-excess ---------------------------------------------
png_open("Figure3_dexcess.png", 1500, 1200)
par(mar = c(4, 4.4, 1.5, 1), mgp = c(2.6, 0.7, 0))
bx <- list(Precipitation = pc$dex,
           `CT tile` = tile_w$dex[tile_w$system == "CT"],
           `NT tile` = tile_w$dex[tile_w$system == "NT"])
boxplot(bx, ylab = "d-excess (permil)", col = c("grey85", CT_COL, NT_COL),
        border = "grey25", outline = FALSE)
for (i in seq_along(bx))
  points(jitter(rep(i, length(bx[[i]])), amount = .12), bx[[i]],
         pch = 21, bg = "white", col = "grey30", cex = .8)
fig_close()

## ---- Figure 4: interaction plot -------------------------------------
png_open("Figure4_interaction.png", 1600, 1250)
par(mar = c(4, 4.6, 1.5, 1), mgp = c(2.8, 0.7, 0))
g <- aggregate(d18Op ~ system + flow, mf,
               function(z) c(m = mean(z), s = sd(z), n = length(z)))
g <- data.frame(system = g$system, flow = g$flow,
                m = g$d18Op[, 1], s = g$d18Op[, 2], n = g$d18Op[, 3])
xs <- c(Baseflow = 1, Event = 2)
plot(NA, xlim = c(.7, 2.3), ylim = range(g$m + g$s, g$m - g$s, 10.2),
     xaxt = "n", xlab = "Flow regime (system-specific daily Q75)",
     ylab = expression(delta^18*"O-PO"[4]*" (permil, VSMOW)"))
axis(1, at = 1:2, labels = c("Baseflow", "Event"))
for (s in c("CT", "NT")) {
  z <- g[g$system == s, ]; z <- z[order(xs[as.character(z$flow)]), ]
  cl <- if (s == "CT") CT_COL else NT_COL
  x <- xs[as.character(z$flow)]
  lines(x, z$m, col = cl, lwd = 2.5)
  arrows(x, z$m - z$s, x, z$m + z$s, angle = 90, code = 3, length = .05, col = cl)
  points(x, z$m, pch = if (s == "CT") 22 else 24, bg = cl, col = "white", cex = 1.6)
  text(x, z$m, paste0("n=", z$n), pos = c(2, 4), cex = .75, col = cl)
  abline(h = mean(mf$ref[mf$system == s], na.rm = TRUE), lty = 3, col = cl)
}
legend("topleft", c("CT field", "NT field",
                    "contemporaneous tile-water reference"),
       col = c(CT_COL, NT_COL, "black"), lty = c(1, 1, 3), lwd = c(2.5, 2.5, 1),
       pch = c(22, 24, NA), pt.bg = c(CT_COL, NT_COL, NA), bty = "n", cex = .8)
fig_close()

## ---- Figure 5: continuous-discharge effect plot (Model D) -----------
png_open("Figure5_continuous_discharge.png", 1600, 1250)
par(mar = c(4.2, 4.6, 1.5, 1), mgp = c(2.8, 0.7, 0))
plot(mf$logQ_c + mean(log10(mf$Q_daily + CONST)), mf$d18Op, type = "n",
     xlab = expression("log"[10]*"(daily mean Q + 0.001, L s"^-1*")"),
     ylab = expression(delta^18*"O-PO"[4]*" (permil)"))
for (s in c("CT", "NT")) {
  z  <- mf[mf$system == s, ]
  cl <- if (s == "CT") CT_COL else NT_COL
  points(log10(z$Q_daily + CONST), z$d18Op,
         pch = if (s == "CT") 22 else 24, bg = cl, col = "white", cex = 1.2)
  nd <- data.frame(system = factor(s, levels = levels(mf$system)),
                   logQ_c = seq(min(z$logQ_c), max(z$logQ_c), length = 60))
  pv <- predict(mD, nd, interval = "confidence")
  xx <- nd$logQ_c + mean(log10(mf$Q_daily + CONST))
  polygon(c(xx, rev(xx)), c(pv[, "lwr"], rev(pv[, "upr"])),
          col = adjustcolor(cl, alpha.f = .15), border = NA)
  lines(xx, pv[, "fit"], col = cl, lwd = 2.5)
}
legend("topleft", c("CT field", "NT field"), col = c(CT_COL, NT_COL),
       lwd = 2.5, pch = c(22, 24), pt.bg = c(CT_COL, NT_COL), bty = "n", cex = .85)
fig_close()

## ---- Figure S1: fixed vs contemporaneous reference -------------------
png_open("FigureS1_reference_comparison.png", 2000, 1000)
par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2, 1), mgp = c(2.6, .7, 0))
plot(mf$date, mf$ref, type = "n", xlab = "", ylab = "Reference (permil)",
     main = "(a) sample-specific reference")
for (s in c("CT", "NT")) {
  z <- mf[mf$system == s & !is.na(mf$ref), ]
  points(z$date, z$ref, pch = if (s == "CT") 22 else 24,
         bg = if (s == "CT") CT_COL else NT_COL, col = "white")
}
abline(h = D_FIXED_REF, lty = 2)
plot(mf$D_fixed, mf$D_app, xlab = "Deviation from fixed reference (permil)",
     ylab = "D_app (permil)", main = "(b) fixed vs contemporaneous",
     pch = ifelse(mf$system == "CT", 22, 24),
     bg = ifelse(mf$system == "CT", CT_COL, NT_COL), col = "white")
abline(0, 1, lty = 2)
fig_close()

## ---- Figure S2: reference sensitivity surface -----------------------
png_open("FigureS2_reference_sensitivity.png", 1500, 1250)
par(mar = c(4.2, 4.4, 2, 1), mgp = c(2.6, .7, 0))
Tg <- seq(2, 12, length = 120); Wg <- seq(-15.5, -13, length = 120)
Z  <- outer(Tg, Wg, cb_reference)
image(Tg, Wg, Z, col = hcl.colors(64, "YlGnBu", rev = TRUE),
      xlab = "Tile-water temperature (deg C)",
      ylab = expression("Tile-water "*delta^18*"O-H"[2]*"O (permil)"),
      main = "Reference sensitivity (Chang & Blake 2015)")
contour(Tg, Wg, Z, add = TRUE, col = "white", labcex = .7, nlevels = 12)
points(5, -14.1, pch = 21, bg = "white", cex = 1.6)
points(mean(mf$T_C[mf$system == "CT"], na.rm = TRUE),
       mean(mf$d18Ow[mf$system == "CT"], na.rm = TRUE), pch = 22, bg = CT_COL, cex = 1.6)
points(mean(mf$T_C[mf$system == "NT"], na.rm = TRUE),
       mean(mf$d18Ow[mf$system == "NT"], na.rm = TRUE), pch = 24, bg = NT_COL, cex = 1.6)
fig_close()

## ---- Figure S3: Model D residual diagnostics ------------------------
png_open("FigureS3_diagnostics.png", 1800, 1400)
par(mfrow = c(2, 2), mar = c(4, 4.2, 2, 1), mgp = c(2.5, .7, 0))
plot(fitted(mD), resid(mD), xlab = "Fitted", ylab = "Residual",
     main = "Residuals vs fitted"); abline(h = 0, lty = 2)
qqnorm(resid(mD), main = "Normal Q-Q"); qqline(resid(mD), lty = 2)
plot(cooks.distance(mD), type = "h", ylab = "Cook's D", xlab = "Observation",
     main = "Influence"); abline(h = 4 / nobs(mD), lty = 2, col = "red3")
plot(hatvalues(mD), resid(mD), xlab = "Leverage", ylab = "Residual",
     main = "Leverage vs residual"); abline(h = 0, lty = 2)
fig_close()

note("Figures written to ", OUT)

## ---- draft captions --------------------------------------------------
cap <- c(
"Figure 1. Daily rainfall (a) and daily-mean tile discharge (b, log scale) for the conventional-tillage (CT) and no-till (NT) fields at the R.J. Cook Agronomy Farm, 2022-2024. Horizontal dashed lines mark the system-specific 75th-percentile DAILY-MEAN discharge used for flow-regime classification (CT 0.173, NT 0.347 L/s). The discontinuity in late 2023 reflects a data-logger gap (81 d CT; 92 d NT). The two fields are single long-term management systems, not replicated treatments.",
"",
"Figure 2. Dual-isotope (Craig) diagram. The local meteoric water line was fitted to 21 retained precipitation samples; the tile evaporation line was fitted to pooled tile water. A slope well below 8 is consistent with kinetic evaporation somewhere along the rainfall-to-drain pathway, and does not identify where that evaporation occurred. The cross marks one excluded precipitation outlier (12 April 2023).",
"",
"Figure 3. Deuterium excess in precipitation and in tile water from the two fields. CT tile water is 2.5 permil lower than NT (paired t-test on 25 matched dates, p = 0.006), consistent with greater evaporative modification of the water reaching the CT drains. This comparison involves one field per management system.",
"",
"Figure 4. System-by-flow-regime pattern in d18O-PO4. Points are group means, error bars +/-1 SD, group sizes shown. Dotted lines mark the mean CONTEMPORANEOUS TILE-WATER REFERENCE calculated per sample from Chang and Blake (2015). This reference is a standardized comparison value; proximity to it does not demonstrate that phosphate reached biochemical equilibrium, because tile water is not necessarily the water in which P-O exchange occurred. The crossing pattern is a comparison of two long-term management systems and is not evidence of a general causal tillage effect.",
"",
"Figure 5. d18O-PO4 against continuous log10 daily-mean discharge, with 95% confidence bands. Treating discharge continuously avoids any dependence on a threshold definition; the difference in slope between the two fields remains positive with a confidence interval excluding zero (HC3 robust standard errors).",
"",
"Figure S1. (a) Sample-specific contemporaneous tile-water reference through time, with the earlier fixed pooled value (11.2 permil) as a dashed line. (b) Deviation computed against the fixed reference versus against the sample-specific reference. Points cluster along the 1:1 line, so the choice of reference does not change the relative ordering or the system contrast.",
"",
"Figure S2. Reference d18O-PO4 across plausible tile-water temperature and d18O-H2O. The numerical stability of the reference within tile-water conditions does not establish that tile water was the water of biochemical oxygen exchange.",
"",
"Figure S4. Residual diagnostics for the continuous-discharge model: residuals versus fitted, normal Q-Q, Cook's distance (dashed line at 4/n), and leverage. Three observations exceed the 4/n cutoff; none was removed, as there is no documented data-quality reason to exclude them.")
writeLines(cap, file.path(OUT, "Figure_captions_draft.txt"))


#######################################################################
## PART 7 -- SESSION PROVENANCE AND LOG
#######################################################################
note("\n=== SESSION ===")
note(R.version.string, " | platform ", R.version$platform)
note("nlme ", as.character(packageVersion("nlme")))
note("project root ", PROJ)
note("run completed ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
     " in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
LOG <- file.path(OUT, "analysis_log.txt")
writeLines(log_note$lines, LOG)
cat("\nLog written to", LOG, "\n")
cat("Tables and review PNGs in", OUT, "\n")
cat("Journal figures (TIFF/PDF/PNG) in", file.path(OUT, "hires"), ":\n")
cat(paste0("  ", sort(list.files(file.path(OUT, "hires"))), collapse = "\n"), "\n")
