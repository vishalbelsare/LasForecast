#### Prepare FRED-MD data for inflation forecasting application
#
# Target (Y): infl_cpi — monthly CPI inflation, 100 * diff(log(CPIAUCSL))
# Predictors (X): all other FRED-MD variables (transformed per McCracken & Ng 2016)
#   - UNRATE is kept in levels (raw) rather than first-differenced, because it is
#     a near-unit-root process relevant to predictive regression methodology
#   - CPIAUCSL is excluded from predictors (used to construct Y)
#
# Outputs:
#   fredmd    — tibble: date, infl_cpi, UNRATE (levels), transformed predictors
#   data_raw  — tibble: date, infl_cpi, UNRATE (levels), untransformed series
#   vars_all  — character vector of predictor variable names (excluding date, infl_cpi)

devtools::load_all()

# Set data range
start_date <- "1959-11-01" # Include two more periods before 1960-1-1 to allow lags and prediction
end_date <- "2025-04-01"

# Read data
data_raw <- fbi::fredmd("data-raw/2025-08-md.csv", date_start = as.Date(start_date), date_end = as.Date(end_date), transform = FALSE)
data_trans <- fbi::fredmd("data-raw/2025-08-md.csv", date_start = as.Date(start_date), date_end = as.Date(end_date))
varlist <- fbi::fredmd_description

## Select Variables based on the varlist provided by fredmd_description
vars_all <- intersect(colnames(data_trans), varlist$fred)

data_trans <- data_trans %>%
  as_tibble() %>%
  select(all_of(c("date", vars_all)))

data_raw <- data_raw %>%
  as_tibble() %>%
  select(all_of(c("date", vars_all)))

data_trans <- data_trans %>%
  filter(date >= start_date & date <= end_date) %>%
  select_if(~ !any(is.na(.)))

data_raw <- data_raw %>%
  filter(date >= start_date & date <= end_date) %>%
  select_if(~ !any(is.na(.))) %>%
  mutate(infl_cpi = c(NA, diff(log(CPIAUCSL)) * 100))


data_trans_full <- data_trans %>%
  select(-UNRATE) %>%
  bind_cols(data_raw %>% select(UNRATE, infl_cpi)) %>%
  select(date, infl_cpi, UNRATE, everything()) %>%
  select(-CPIAUCSL) %>%
  slice(-1)

data_raw_full <- data_raw %>%
  select(date, infl_cpi, UNRATE, everything()) %>%
  select(-CPIAUCSL) %>%
  slice(-1)


fredmd <- data_trans_full
data_raw <- data_raw_full
vars_all <- setdiff(colnames(fredmd), c("date", "infl_cpi"))

usethis::use_data(fredmd, overwrite = TRUE)
usethis::use_data(data_raw, overwrite = TRUE)
usethis::use_data(vars_all, overwrite = TRUE)
