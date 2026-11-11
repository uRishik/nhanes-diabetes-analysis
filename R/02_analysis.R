# 02_analysis.R
# Survey-weighted analysis of diabetes prevalence and its correlates among
# US adults, NHANES August 2021-August 2023.

suppressPackageStartupMessages({
  library(survey)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

options(survey.lonely.psu = "adjust")

dat <- readRDS("data/processed/nhanes_2021_2023.rds")

# Build the design on ALL participants, then subset. Subsetting the design
# (rather than the data) keeps the full stratum/PSU structure, which gives
# correct standard errors for subpopulation (domain) estimates.
des_full <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt,
  nest = TRUE, data = dat
)
des <- subset(des_full, in_sample)

fmt_pct <- function(x) sprintf("%.1f", 100 * x)

# ---------------------------------------------------------------------------
# 1. Overall prevalence: weighted vs. unweighted
# ---------------------------------------------------------------------------
prev_ci <- function(formula, design) {
  est <- svyciprop(formula, design, method = "logit")
  ci <- attr(est, "ci")
  c(est = as.numeric(est), lo = ci[[1]], hi = ci[[2]])
}

overall <- bind_rows(
  tibble(outcome = "Total diabetes",       !!!prev_ci(~total_diabetes, des)),
  tibble(outcome = "Diagnosed diabetes",   !!!prev_ci(~dx_diabetes, des)),
  tibble(outcome = "Undiagnosed diabetes", !!!prev_ci(~undx_diabetes, des))
)

samp <- filter(dat, in_sample)
overall$unweighted <- c(mean(samp$total_diabetes),
                        mean(samp$dx_diabetes),
                        mean(samp$undx_diabetes))

overall_out <- overall |>
  transmute(
    Outcome = outcome,
    `Weighted % (95% CI)` = sprintf("%s (%s-%s)", fmt_pct(est), fmt_pct(lo), fmt_pct(hi)),
    `Unweighted %` = fmt_pct(unweighted)
  )
write.csv(overall_out, "output/tables/overall_prevalence.csv", row.names = FALSE)
print(overall_out)

# Share of diabetes that is undiagnosed
des_dm <- subset(des, total_diabetes == 1)
undx_share <- prev_ci(~undx_diabetes, des_dm)
cat(sprintf("\nUndiagnosed share of all diabetes: %s%% (%s-%s)\n",
            fmt_pct(undx_share["est"]), fmt_pct(undx_share["lo"]), fmt_pct(undx_share["hi"])))
write.csv(
  tibble(measure = "Undiagnosed share of total diabetes",
         estimate = undx_share["est"], lower = undx_share["lo"], upper = undx_share["hi"]),
  "output/tables/undiagnosed_share.csv", row.names = FALSE
)

# ---------------------------------------------------------------------------
# 2. Table 1: weighted characteristics by diabetes status
# ---------------------------------------------------------------------------
vars <- c(age_group = "Age group", sex = "Sex", race_eth = "Race/ethnicity",
          bmi_cat = "BMI (kg/m2)", pir_cat = "Income-to-poverty ratio",
          education = "Education")

table1 <- lapply(names(vars), function(v) {
  f <- as.formula(paste0("~", v))
  by_dm <- svyby(f, ~total_diabetes, des, svymean, na.rm = TRUE)
  overall_v <- svymean(f, des, na.rm = TRUE)
  levs <- levels(dat[[v]])
  n_by <- table(samp[[v]])
  tibble(
    Variable = vars[[v]],
    Level = levs,
    `Unweighted n` = as.integer(n_by[levs]),
    `Overall %` = fmt_pct(coef(overall_v)[paste0(v, levs)]),
    `No diabetes %` = fmt_pct(unlist(by_dm[by_dm$total_diabetes == 0, paste0(v, levs)])),
    `Diabetes %` = fmt_pct(unlist(by_dm[by_dm$total_diabetes == 1, paste0(v, levs)]))
  )
}) |> bind_rows()

age_mean <- svyby(~age, ~total_diabetes, des, svymean)
table1 <- bind_rows(
  tibble(Variable = "Age, mean (years)", Level = "",
         `Unweighted n` = nrow(samp),
         `Overall %` = sprintf("%.1f", coef(svymean(~age, des))),
         `No diabetes %` = sprintf("%.1f", age_mean$age[age_mean$total_diabetes == 0]),
         `Diabetes %` = sprintf("%.1f", age_mean$age[age_mean$total_diabetes == 1])),
  table1
)
write.csv(table1, "output/tables/table1_characteristics.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Prevalence by subgroup
# ---------------------------------------------------------------------------
subgroup_prev <- lapply(names(vars), function(v) {
  levs <- levels(dat[[v]])
  lapply(levs, function(l) {
    d <- subset(des, get(v) == l)
    r <- prev_ci(~total_diabetes, d)
    tibble(Variable = vars[[v]], Level = l, est = r["est"], lo = r["lo"], hi = r["hi"])
  }) |> bind_rows()
}) |> bind_rows()

write.csv(
  subgroup_prev |> mutate(across(c(est, lo, hi), \(x) round(100 * x, 1))),
  "output/tables/prevalence_by_subgroup.csv", row.names = FALSE
)

# ---------------------------------------------------------------------------
# 4. Survey-weighted logistic regression
# ---------------------------------------------------------------------------
# Complete-case design for the adjusted model
des_cc <- subset(des, !is.na(bmi_cat) & !is.na(pir_cat) & !is.na(education))

m_unadj <- svyglm(total_diabetes ~ bmi_cat, design = des_cc, family = quasibinomial())
m_adj <- svyglm(
  total_diabetes ~ bmi_cat + age_group + sex + race_eth + pir_cat + education,
  design = des_cc, family = quasibinomial()
)

# Inference uses the DESIGN degrees of freedom (#PSUs - #strata = 15), as
# recommended by NCHS. svyglm's default residual df subtracts the number of
# model parameters from the design df, which leaves ~0 df for a model with
# 14 coefficients and makes t-based CIs undefined.
design_df <- degf(des_cc)

tidy_or <- function(m, label, df = design_df) {
  b <- coef(m)
  se <- sqrt(diag(vcov(m)))
  t_crit <- qt(0.975, df)
  tibble(
    term = names(b),
    estimate = exp(b),
    conf.low = exp(b - t_crit * se),
    conf.high = exp(b + t_crit * se),
    p.value = 2 * pt(abs(b / se), df, lower.tail = FALSE)
  ) |>
    filter(term != "(Intercept)") |>
    mutate(model = label)
}
ors <- bind_rows(tidy_or(m_unadj, "Unadjusted"), tidy_or(m_adj, "Adjusted"))

pretty_term <- function(t) {
  t |>
    sub("^bmi_cat", "BMI: ", x = _) |>
    sub("^age_group", "Age: ", x = _) |>
    sub("^sex", "Sex: ", x = _) |>
    sub("^race_eth", "Race/eth: ", x = _) |>
    sub("^pir_cat", "Income/poverty: ", x = _) |>
    sub("^education", "Education: ", x = _)
}

or_table <- ors |>
  mutate(Term = pretty_term(term),
         `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
         p = format.pval(p.value, digits = 2, eps = 0.001)) |>
  select(model, Term, `OR (95% CI)`, p)
write.csv(or_table, "output/tables/logistic_regression_ORs.csv", row.names = FALSE)
print(or_table, n = Inf)

cat(sprintf("\nComplete-case n for adjusted model: %d of %d\n",
            nrow(des_cc$variables), nrow(des$variables)))

# Design-based Wald test for overall BMI effect in the adjusted model
bmi_test <- regTermTest(m_adj, ~bmi_cat, method = "Wald", df = design_df)
print(bmi_test)
writeLines(capture.output(print(bmi_test)), "output/tables/bmi_wald_test.txt")

# Model discrimination (weighted AUC via design-weighted concordance is not
# built into survey; report unweighted AUC as a rough diagnostic)
p_hat <- fitted(m_adj)
y <- des_cc$variables$total_diabetes
auc <- {
  r <- rank(p_hat); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cat(sprintf("Unweighted AUC (adjusted model): %.3f\n", auc))
writeLines(sprintf("%.3f", auc), "output/tables/adjusted_model_auc.txt")

# ---------------------------------------------------------------------------
# 5. Figures
# ---------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title.position = "plot"))

# Figure 1: prevalence by BMI category and age group
by_bmi_age <- svyby(~total_diabetes, ~bmi_cat + age_group, des,
                    svyciprop, vartype = "ci", method = "logit", na.rm = TRUE)
p1 <- ggplot(by_bmi_age,
             aes(bmi_cat, total_diabetes, colour = age_group, group = age_group)) +
  geom_pointrange(aes(ymin = ci_l, ymax = ci_u),
                  position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_manual(values = c("#7fb3d5", "#2e86c1", "#1b4f72")) +
  labs(title = "Diabetes prevalence rises with BMI in every age group",
       subtitle = "US adults, NHANES 2021-2023; survey-weighted with 95% CIs",
       x = "BMI category (kg/m²)", y = "Total diabetes prevalence",
       colour = "Age group")
ggsave("output/figures/fig1_prevalence_bmi_age.png", p1, width = 8, height = 5, dpi = 300, bg = "white")

# Figure 2: forest plot of adjusted ORs
forest <- ors |>
  filter(model == "Adjusted") |>
  mutate(Term = pretty_term(term),
         Term = factor(Term, levels = rev(Term)))
p2 <- ggplot(forest, aes(estimate, Term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2,
                orientation = "y", colour = "#1b4f72") +
  geom_point(size = 2.5, colour = "#1b4f72") +
  scale_x_log10() +
  labs(title = "Adjusted odds ratios for total diabetes",
       subtitle = "Survey-weighted logistic regression; reference groups: BMI <25, age 20-39,\nmale, non-Hispanic White, income/poverty <1.3, college graduate+",
       x = "Odds ratio (log scale, 95% CI)", y = NULL)
ggsave("output/figures/fig2_forest_adjusted_or.png", p2, width = 8, height = 6.5, dpi = 300, bg = "white")

# Figure 3: diagnosed vs undiagnosed by age group
dx_age <- svyby(~dx_diabetes + undx_diabetes, ~age_group, des, svymean) |>
  as_tibble() |>
  select(age_group, Diagnosed = dx_diabetes, Undiagnosed = undx_diabetes) |>
  pivot_longer(-age_group, names_to = "status", values_to = "prev")
p3 <- ggplot(dx_age, aes(age_group, prev, fill = status)) +
  geom_col(width = 0.6) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(Diagnosed = "#1b4f72", Undiagnosed = "#e59866")) +
  labs(title = "Diagnosed and undiagnosed diabetes by age",
       subtitle = "Undiagnosed = no prior diagnosis but HbA1c ≥ 6.5%",
       x = "Age group", y = "Weighted prevalence", fill = NULL)
ggsave("output/figures/fig3_dx_undx_by_age.png", p3, width = 7, height = 5, dpi = 300, bg = "white")

message("Analysis complete. Tables in output/tables, figures in output/figures.")
