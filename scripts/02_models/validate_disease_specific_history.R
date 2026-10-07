#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glmmTMB)
})

data <- read_csv("output/revised/analysis_data.csv", show_col_types = FALSE)
required <- c("country", "disease", "year", "cases", "population", "coverage_10",
              "prior", "gdp", "uhc", "urban", "health", "time", "hist_10", "gap_10")
stopifnot(all(required %in% names(data)), !anyDuplicated(data[c("country", "disease", "year")]))

base <- cases ~ coverage_10 + prior + gdp + uhc + urban + health + time +
  offset(log(population)) + (1 | country)
fit_checked <- function(formula, rows) {
  model <- glmmTMB(formula, data = rows, family = nbinom2,
                   control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000)))
  if (model$fit$convergence != 0 || !isTRUE(model$sdr$pdHess)) {
    model <- glmmTMB(formula, data = rows, family = nbinom2,
                     control = glmmTMBControl(optimizer = optim,
                                              optArgs = list(method = "BFGS"),
                                              optCtrl = list(maxit = 4000)))
  }
  stopifnot(model$fit$convergence == 0, isTRUE(model$sdr$pdHess), is.finite(AIC(model)))
  model
}

predictions <- list()
for (disease_name in sort(unique(data$disease))) {
  rows <- filter(data, disease == disease_name)
  for (cutoff in c(2016L, 2019L)) {
    train <- filter(rows, year <= cutoff)
    test <- filter(rows, year > cutoff, year <= if (cutoff == 2016L) 2019L else 2023L)
    stopifnot(nrow(train) > 0, nrow(test) > 0,
              all(is.finite(test$cases)), all(test$population > 0))
    for (model_name in c("M1_history", "M3_cohort_time")) {
      exposure <- if (model_name == "M1_history") "hist_10" else "gap_10"
      model <- fit_checked(update(base, paste("~ . +", exposure)), train)
      predicted <- predict(model, newdata = test, type = "response", allow.new.levels = TRUE)
      dispersion <- sigma(model)
      stopifnot(length(predicted) == nrow(test), all(is.finite(predicted)), all(predicted >= 0))
      predictions[[length(predictions) + 1L]] <- tibble(
        disease = disease_name, cutoff, model = model_name,
        country = test$country, year = test$year, actual = test$cases,
        prediction = predicted, dispersion, n_train = nrow(train))
    }
  }
}

predictions <- bind_rows(predictions)
metrics <- predictions %>% group_by(disease, cutoff, model) %>% summarise(
  n_train = first(n_train), n_test = n(),
  MAE = mean(abs(actual - prediction)),
  RMSE = sqrt(mean((actual - prediction)^2)),
  log_MSE = mean((log1p(actual) - log1p(prediction))^2),
  mean_log_predictive_density = mean(dnbinom(actual, mu = prediction,
                                            size = dispersion, log = TRUE)),
  .groups = "drop")
comparison <- metrics %>% select(disease, cutoff, model, RMSE, MAE, log_MSE,
                                 mean_log_predictive_density) %>%
  tidyr::pivot_wider(names_from = model,
                     values_from = c(RMSE, MAE, log_MSE, mean_log_predictive_density)) %>%
  mutate(RMSE_change_pct = 100 * (RMSE_M3_cohort_time / RMSE_M1_history - 1),
         MAE_change_pct = 100 * (MAE_M3_cohort_time / MAE_M1_history - 1))

dir.create("output/revised/disease_specific_holdout", recursive = TRUE, showWarnings = FALSE)
dir.create("output/revised/supplementary_tables", recursive = TRUE, showWarnings = FALSE)
write_csv(predictions, "output/revised/disease_specific_holdout/predictions.csv")
write_csv(metrics, "output/revised/disease_specific_holdout/metrics.csv")
write_csv(comparison, "output/revised/disease_specific_holdout/comparison.csv")
write_csv(metrics, "output/revised/supplementary_tables/TableS46_disease_specific_time_holdout.csv")
print(comparison, n = Inf)
