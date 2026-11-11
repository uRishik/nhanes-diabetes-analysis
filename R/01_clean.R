# 01_clean.R
# Merge NHANES 2021-2023 files, define outcomes and covariates, and save an
# analysis-ready data set. All participants are kept (not just the analytic
# sample) so that subpopulation analyses use the full survey design.

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
})

read_raw <- function(name) read_xpt(file.path("data/raw", paste0(name, ".xpt")))

demo <- read_raw("DEMO_L")
ghb  <- read_raw("GHB_L")
diq  <- read_raw("DIQ_L")
bmx  <- read_raw("BMX_L")

dat <- demo |>
  left_join(ghb, by = "SEQN") |>
  left_join(diq, by = "SEQN") |>
  left_join(bmx, by = "SEQN")

dat <- dat |>
  mutate(
    # ---- outcomes -----------------------------------------------------------
    hba1c = LBXGH,
    # DIQ010: 1 = yes, 2 = no, 3 = borderline; 7/9 = refused/don't know
    dx_diabetes = case_when(
      DIQ010 == 1 ~ 1L,
      DIQ010 %in% c(2, 3) ~ 0L,
      TRUE ~ NA_integer_
    ),
    undx_diabetes = case_when(
      dx_diabetes == 0 & hba1c >= 6.5 ~ 1L,
      dx_diabetes == 0 & hba1c <  6.5 ~ 0L,
      dx_diabetes == 1 ~ 0L,
      TRUE ~ NA_integer_
    ),
    total_diabetes = as.integer(dx_diabetes == 1 | undx_diabetes == 1),

    # ---- covariates ---------------------------------------------------------
    age = RIDAGEYR,
    age_group = cut(age, c(19, 39, 59, Inf),
                    labels = c("20-39", "40-59", "60+")),
    sex = factor(RIAGENDR, 1:2, c("Male", "Female")),
    race_eth = factor(
      case_when(
        RIDRETH3 %in% c(1, 2) ~ "Hispanic",
        RIDRETH3 == 3 ~ "Non-Hispanic White",
        RIDRETH3 == 4 ~ "Non-Hispanic Black",
        RIDRETH3 == 6 ~ "Non-Hispanic Asian",
        RIDRETH3 == 7 ~ "Other/Multiracial"
      ),
      levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic",
                 "Non-Hispanic Asian", "Other/Multiracial")
    ),
    bmi = BMXBMI,
    bmi_cat = cut(bmi, c(-Inf, 25, 30, 35, Inf), right = FALSE,
                  labels = c("<25", "25-29.9", "30-34.9", ">=35")),
    pir_cat = cut(INDFMPIR, c(-Inf, 1.3, 3.5, Inf), right = FALSE,
                  labels = c("<1.3", "1.3-3.49", ">=3.5")),
    education = factor(
      case_when(
        DMDEDUC2 %in% c(1, 2) ~ "Less than high school",
        DMDEDUC2 == 3 ~ "High school/GED",
        DMDEDUC2 == 4 ~ "Some college",
        DMDEDUC2 == 5 ~ "College graduate+"
      ),
      levels = c("College graduate+", "Some college", "High school/GED",
                 "Less than high school")
    ),

    # ---- survey design ------------------------------------------------------
    # HbA1c is measured in the phlebotomy subsample, so the phlebotomy weight
    # (WTPH2YR) is the correct weight. People outside that subsample get 0.
    wt = if_else(is.na(WTPH2YR), 0, WTPH2YR),

    # ---- analytic sample flag -----------------------------------------------
    pregnant = !is.na(RIDEXPRG) & RIDEXPRG == 1,
    in_sample = age >= 20 & !pregnant & wt > 0 & !is.na(hba1c) &
                !is.na(dx_diabetes)
  )

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(dat, "data/processed/nhanes_2021_2023.rds")

# Sample-flow numbers for the README
flow <- tibble::tibble(
  step = c("NHANES 2021-2023 participants",
           "Age >= 20",
           "Not pregnant",
           "In phlebotomy subsample with HbA1c",
           "Answered diabetes question (analytic sample)"),
  n = c(
    nrow(dat),
    sum(dat$age >= 20),
    sum(dat$age >= 20 & !dat$pregnant),
    sum(dat$age >= 20 & !dat$pregnant & dat$wt > 0 & !is.na(dat$hba1c)),
    sum(dat$in_sample)
  )
)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(flow, "output/tables/sample_flow.csv", row.names = FALSE)
print(flow)
