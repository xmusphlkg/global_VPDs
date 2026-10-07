# ==============================================================================
# 01_data_download.R
# 从 WHO Immunization Data Portal 和 UN WPP 2024 下载并整理项目数据
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
})

dir.create(here("data", "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)

# 下载地址（访问日期和完整来源记录见 data/source.md）
coverage_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/coverage-data.xlsx"
)
cases_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/reported-cases-data.xlsx"
)
population_url <- paste0(
  "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/",
  "CSV_FILES/WPP2024_TotalPopulationBySex.csv.gz"
)
incidence_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/incidence-rate-data.xlsx"
)
introduction_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/vaccine-introduction-data.xlsx"
)
schedule_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/vaccine-schedule-data.xlsx"
)
other_indicators_url <- paste0(
  "https://srhdpeuwpubsa-geecgzbpd5h0fueu.z01.azurefd.net/",
  "whdh/WIISE/export/other-indicators-data.xlsx"
)
mr_reporting_url <- "https://immunizationdata.who.int/docs/librariesprovider21/measles-and-rubella/403-table-web-reporting-data.xlsx?sfvrsn=9409ffb9_10"
mr_epi_url <- "https://immunizationdata.who.int/docs/librariesprovider21/measles-and-rubella/404-table-web-epi-curve-data.xlsx?sfvrsn=5922ebf7_20"
coverage_survey_url <- "https://immunizationdata.who.int/docs/librariesprovider21/additional-datasets/coverage-survey-data.xlsx?sfvrsn=cafe05bd_1"

download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) {
    message("下载: ", url)
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  } else {
    message("已存在，跳过下载: ", dest)
  }
}

coverage_file <- here("data", "raw", "who_immunization_portal_coverage.xlsx")
cases_file <- here("data", "raw", "who_immunization_portal_reported_cases.xlsx")
population_file <- here("data", "raw", "WPP2024_TotalPopulationBySex.csv.gz")
incidence_file <- here("data", "raw", "who_immunization_portal_incidence_rate.xlsx")
introduction_file <- here("data", "raw", "who_vaccine-introduction-data.xlsx")
schedule_file <- here("data", "raw", "who_vaccine-schedule-data.xlsx")
other_indicators_file <- here("data", "raw", "who_other-indicators-data.xlsx")
mr_reporting_file <- here("data", "raw", "who_403-table-web-reporting-data.xlsx")
mr_epi_file <- here("data", "raw", "who_404-table-web-epi-curve-data.xlsx")
coverage_survey_file <- here("data", "raw", "who_coverage-survey-data.xlsx")

download_if_missing(coverage_url, coverage_file)
download_if_missing(cases_url, cases_file)
download_if_missing(population_url, population_file)
download_if_missing(incidence_url, incidence_file)
download_if_missing(introduction_url, introduction_file)
download_if_missing(schedule_url, schedule_file)
download_if_missing(other_indicators_url, other_indicators_file)
download_if_missing(mr_reporting_url, mr_reporting_file)
download_if_missing(mr_epi_url, mr_epi_file)
download_if_missing(coverage_survey_url, coverage_survey_file)

# 项目使用的疫苗-疾病映射
default_vaccine_disease_mapping <- tribble(
  ~disease,    ~disease_code, ~vaccine, ~vaccine_code,
  "Measles",   "MEASLES",    "MCV1",  "MCV1",
  "Measles",   "MEASLES",    "MCV2",  "MCV2",
  "Rubella",   "RUBELLA",   "RCV1",  "RCV1",
  "Pertussis", "PERTUSSIS", "DTP3",  "DTP3",
  "Diphtheria","DIPHTHERIA","DTP3",  "DTP3"
)
mapping_file <- here("data", "processed", "vaccine_disease_mapping.csv")
if (!file.exists(mapping_file)) {
  write_csv(default_vaccine_disease_mapping, mapping_file)
} else {
  message("保留已有疫苗-疾病映射: ", mapping_file)
}

# WHO 免疫数据门户覆盖率导出。门户使用 DTPCV3 表示 DTP3。
coverage <- read_excel(coverage_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(
    group = GROUP,
    country = CODE,
    year = as.integer(YEAR),
    vaccine = recode(ANTIGEN, DTPCV3 = "DTP3", .default = ANTIGEN),
    coverage_category = COVERAGE_CATEGORY,
    coverage = as.numeric(COVERAGE)
  ) %>%
  filter(
    group == "COUNTRIES",
    year >= 2000, year <= 2025,
    vaccine %in% c("MCV1", "MCV2", "RCV1", "DTP3"),
    coverage_category == "WUENIC",
    !is.na(coverage)
  ) %>%
  select(country, year, vaccine, coverage) %>%
  arrange(country, year, vaccine)
write_csv(coverage, here("data", "raw", "coverage_data.csv"))

# 更完整的 WUENIC 覆盖率（所有抗原）和国家报告覆盖率。
coverage_all <- read_excel(coverage_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(
    group = GROUP, country = CODE, country_name = NAME,
    year = as.integer(YEAR),
    antigen = recode(ANTIGEN, DTPCV1 = "DTP1", DTPCV3 = "DTP3", .default = ANTIGEN),
    antigen_description = ANTIGEN_DESCRIPTION,
    coverage_category = COVERAGE_CATEGORY,
    coverage_category_description = COVERAGE_CATEGORY_DESCRIPTION,
    target_number = as.numeric(TARGET_NUMBER), doses = as.numeric(DOSES),
    coverage = as.numeric(COVERAGE)
  ) %>%
  filter(group == "COUNTRIES", year >= 2000, year <= 2025,
         !is.na(coverage))
write_csv(
  coverage_all %>% filter(coverage_category == "WUENIC") %>%
    transmute(country, year, vaccine = antigen, coverage),
  here("data", "raw", "coverage_data.csv")
)
write_csv(
  coverage_all %>% filter(coverage_category == "WUENIC",
                          antigen %in% c("MCV1", "MCV2", "RCV1", "DTP3")) %>%
    transmute(country, year, vaccine = antigen, coverage),
  here("data", "raw", "coverage_core_data.csv")
)
write_csv(
  coverage_all %>% filter(coverage_category == "WUENIC"),
  here("data", "raw", "coverage_wuenic_all.csv")
)
write_csv(
  coverage_all %>% filter(coverage_category %in% c("ADMIN", "OFFICIAL")),
  here("data", "raw", "coverage_reported_all.csv")
)
write_csv(
  coverage_all %>%
    filter(coverage_category %in% c("ADMIN", "OFFICIAL"),
           antigen %in% c("MCV1", "MCV2", "RCV1", "DTP3")),
  here("data", "raw", "coverage_core_reported.csv")
)

# WHO 免疫数据门户报告病例导出。空白 CASES 保留为缺失值并排除。
cases <- read_excel(cases_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(
    group = GROUP,
    country = CODE,
    year = as.integer(YEAR),
    disease = recode(
      DISEASE,
      MEASLES = "Measles",
      RUBELLA = "Rubella",
      PERTUSSIS = "Pertussis",
      DIPHTHERIA = "Diphtheria",
      .default = NA_character_
    ),
    cases = as.numeric(CASES)
  ) %>%
  filter(
    group == "COUNTRIES",
    year >= 2000, year <= 2025,
    !is.na(disease), !is.na(cases)
  ) %>%
  select(country, year, disease, cases) %>%
  arrange(country, year, disease)
write_csv(cases, here("data", "raw", "case_data.csv"))

# UN WPP 2024 总人口（7 月 1 日，Medium variant）。原始单位为千人，转换为人。
population <- read_csv(gzfile(population_file), show_col_types = FALSE,
                       name_repair = "minimal") %>%
  filter(
    LocTypeName == "Country/Area",
    Variant == "Medium",
    Time >= 2000, Time <= 2025,
    !is.na(ISO3_code), ISO3_code != ""
  ) %>%
  transmute(
    country = ISO3_code,
    year = as.integer(Time),
    population = as.numeric(PopTotal) * 1000
  ) %>%
  filter(!is.na(population)) %>%
  arrange(country, year)
write_csv(population, here("data", "raw", "population_data.csv"))

# WHO 免疫门户直接提供的发病率，可用于检查 cases / population 自行计算的结果。
incidence <- read_excel(incidence_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(
    group = GROUP, country = CODE, country_name = NAME,
    year = as.integer(YEAR), disease_code = DISEASE,
    disease = DISEASE_DESCRIPTION, denominator = DENOMINATOR,
    incidence_rate = as.numeric(INCIDENCE_RATE)
  ) %>%
  filter(group == "COUNTRIES", year >= 2000, year <= 2025,
         disease_code %in% c("MEASLES", "RUBELLA", "PERTUSSIS", "DIPHTHERIA"))
write_csv(incidence, here("data", "raw", "incidence_rate_data.csv"))

# 疫苗引入与国家免疫程序：用于区分“尚未引入”与“覆盖率缺失”。
introduction <- read_excel(introduction_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(country = ISO_3_CODE, country_name = COUNTRYNAME,
            who_region = WHO_REGION, year = as.integer(YEAR),
            vaccine = DESCRIPTION, introduced = INTRO) %>%
  filter(year >= 2000, year <= 2025)
write_csv(introduction, here("data", "raw", "vaccine_introduction_data.csv"))
write_csv(
  introduction %>% filter(str_detect(
    vaccine, regex("Measles|Rubella|pertussis|Pneumococcal|Polio|Hepatitis|Hib|Rotavirus|Japanese|Meningococcal|Yellow Fever",
                   ignore_case = TRUE))),
  here("data", "raw", "vaccine_introduction_core.csv")
)

schedule <- read_excel(schedule_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(country = ISO_3_CODE, country_name = COUNTRYNAME,
            who_region = WHO_REGION, year = as.integer(YEAR),
            vaccine_code = VACCINECODE, vaccine = VACCINE_DESCRIPTION,
            schedule_round = as.numeric(SCHEDULEROUNDS),
            target_population = TARGETPOP_DESCRIPTION,
            geographic_area = GEOAREA, age_administered = AGEADMINISTERED,
            source_comment = SOURCECOMMENT) %>%
  filter(year >= 2000, year <= 2025)
write_csv(schedule, here("data", "raw", "vaccine_schedule_data.csv"))
write_csv(
  schedule %>% filter(str_detect(vaccine, regex("measles|rubella|pertussis|diphtheria|dtp|mcv|rcv", ignore_case = TRUE))),
  here("data", "raw", "vaccine_schedule_core.csv")
)

# 系统绩效和监测质量指标：选择与接种缺口、延迟接种和供应中断直接相关的指标。
other_indicators <- read_excel(other_indicators_file, sheet = "Data", .name_repair = "minimal") %>%
  transmute(country = ISO_3_CODE, country_name = COUNTRYNAME,
            who_region = WHO_REGION, year = as.integer(YEAR),
            indicator_code = INDCODE, description = DESCRIPTION,
            indicator_category = INDCAT_DESCRIPTION,
            indicator_sort = INDSORT, value = VALUE) %>%
  filter(year >= 2000, year <= 2025)
selected_indicator_codes <- c(
  "CALC_DTP1_DTP3", "CALC_DTP1_MCV1", "CALC_MCV1_MCV2",
  "CALC_DIST_80DTP3_P", "CALC_DIST_90DTP3_P", "CALC_DIST_49DTP3_P",
  "CALC_DISTR_NotReporting_DTP3_P", "CALC_DIST_80MCV1_P",
  "CALC_DIST_90MCV1_P", "CALC_DIST_49MCV1_P",
  "CALC_DISTR_NotReporting_MCV1_P", "CALC_DIST_80MCV2_P",
  "CALC_DIST_90MCV2_P", "CALC_DIST_49MCV2_P",
  "CALC_DISTR_NotReporting_MCV2_P", "SUPPLY_STOCKOUT",
  "SUPPLY_ANYDIST_DTPCV", "SUPPLY_ANYDIST_MEASLES_CV",
  "SUPPLY_NATIONAL_DTPCV", "SUPPLY_NATIONAL_MEASLES_CV",
  "SUPPLY_MM_DTPCV", "SUPPLY_MM_MEASLES_CV",
  "IS_RECORD_DELAYEDDOSE", "SAFETY_AEFI_CASES"
)
write_csv(
  other_indicators %>% filter(indicator_code %in% selected_indicator_codes),
  here("data", "raw", "monitoring_system_indicators.csv")
)

# WUENIC 使用的调查数据，可用于评估估计值与调查证据的一致性。
# Excel 第二个工作表第一行是标题，第二行才是字段名。
coverage_survey <- read_excel(coverage_survey_file, sheet = 2, skip = 1,
                              .name_repair = "minimal") %>%
  rename(
    country = ISO3,
    survey_name = surveyNameEnglish,
    survey_name_production = surveyNameProduction,
    collect_begin = collectBegin,
    collect_end = collectEnd,
    cohort_year = cohortYear,
    vaccine = vaccine,
    evidence = evidence,
    cards_seen = cardsSeen,
    validity = validity,
    age_vaccination = ageVaccination,
    age_interview = ageInterview,
    coverage = coverage,
    sample_size = Sample_Size
  ) %>%
  mutate(cohort_year = as.integer(cohort_year), coverage = as.numeric(coverage),
         sample_size = as.numeric(sample_size)) %>%
  filter(!is.na(country), country != "", !is.na(cohort_year))
write_csv(coverage_survey, here("data", "raw", "coverage_survey_data.csv"))

# 麻疹/风疹月度病例分类和排除病例，用于监测质量敏感性分析。
mr_epi <- read_excel(mr_epi_file, sheet = 2, .name_repair = "minimal")
names(mr_epi) <- c(
  "region", "country_name", "country", "year", "month",
  "measles_suspect", "measles_clinical", "measles_epi_linked",
  "measles_lab_confirmed", "measles_total", "rubella_clinical",
  "rubella_epi_linked", "rubella_lab_confirmed", "rubella_total",
  "discarded"
)
write_csv(mr_epi, here("data", "raw", "measles_rubella_epi_curve_data.csv"))

# 数据清单（原始文件和规范化快照）
data_inventory <- tribble(
  ~dataset, ~filename, ~source, ~required, ~downloaded,
  "WUENIC 覆盖率（原始）", "who_immunization_portal_coverage.xlsx", "WHO Immunization Data Portal", "是", file.exists(coverage_file),
  "报告病例（原始）", "who_immunization_portal_reported_cases.xlsx", "WHO Immunization Data Portal", "是", file.exists(cases_file),
  "人口（原始）", "WPP2024_TotalPopulationBySex.csv.gz", "UN WPP 2024", "是", file.exists(population_file),
  "疫苗覆盖率（规范化）", "coverage_data.csv", "由 WHO 原始导出筛选", "是", file.exists(here("data", "raw", "coverage_data.csv")),
  "病例报告（规范化）", "case_data.csv", "由 WHO 原始导出筛选", "是", file.exists(here("data", "raw", "case_data.csv")),
  "人口（规范化）", "population_data.csv", "由 UN WPP 2024 筛选", "是", file.exists(here("data", "raw", "population_data.csv")),
  "所有抗原 WUENIC 覆盖率", "coverage_wuenic_all.csv", "由 WHO 覆盖率导出筛选", "是", file.exists(here("data", "raw", "coverage_wuenic_all.csv")),
  "核心抗原 WUENIC 覆盖率", "coverage_core_data.csv", "由 WHO 覆盖率导出筛选", "是", file.exists(here("data", "raw", "coverage_core_data.csv")),
  "国家报告覆盖率", "coverage_reported_all.csv", "由 WHO 覆盖率导出筛选", "是", file.exists(here("data", "raw", "coverage_reported_all.csv")),
  "WHO 发病率", "incidence_rate_data.csv", "由 WHO 发病率导出筛选", "是", file.exists(here("data", "raw", "incidence_rate_data.csv")),
  "疫苗引入", "vaccine_introduction_data.csv", "WHO Immunization Data Portal", "是", file.exists(here("data", "raw", "vaccine_introduction_data.csv")),
  "国家免疫程序", "vaccine_schedule_data.csv", "WHO Immunization Data Portal", "是", file.exists(here("data", "raw", "vaccine_schedule_data.csv")),
  "监测和系统指标", "monitoring_system_indicators.csv", "由 WHO 其他指标导出筛选", "是", file.exists(here("data", "raw", "monitoring_system_indicators.csv")),
  "覆盖率调查", "coverage_survey_data.csv", "WHO/UNICEF WUENIC 调查数据库", "否", file.exists(here("data", "raw", "coverage_survey_data.csv")),
  "麻疹风疹月度监测", "measles_rubella_epi_curve_data.csv", "WHO 麻疹和风疹监测", "否", file.exists(here("data", "raw", "measles_rubella_epi_curve_data.csv")),
  "WPP 出生人数", "births_data.csv", "UN WPP 2024（由外部数据整理流程下载整理）", "是", file.exists(here("data", "raw", "births_data.csv")),
  "WPP 单岁年龄人口", "population_by_single_age.csv", "UN WPP 2024（由外部数据整理流程下载整理）", "是", file.exists(here("data", "raw", "population_by_single_age.csv")),
  "WPP 年龄组人口", "population_age_groups.csv", "由 WPP 单岁人口聚合", "是", file.exists(here("data", "raw", "population_age_groups.csv")),
  "国家—年份人口和出生队列", "country_year_demography.csv", "WPP 出生人数与年龄组人口整合", "是", file.exists(here("data", "raw", "country_year_demography.csv")),
  "国家层面协变量", "country_covariates.csv", "World Bank WDI/WHO UHC", "否", file.exists(here("data", "raw", "country_covariates.csv")),
  "统一国家范围", "analysis_country_universe.csv", "WPP 与 WHO 数据交集", "是", file.exists(here("data", "processed", "analysis_country_universe.csv"))
)
write_csv(data_inventory, here("data", "raw", "data_inventory.csv"))

message("\n数据下载和整理完成：")
message("  覆盖率: ", nrow(coverage), " 行")
message("  病例: ", nrow(cases), " 行")
message("  人口: ", nrow(population), " 行")
message("  详细来源: data/source.md")
