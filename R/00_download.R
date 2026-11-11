# 00_download.R
# Download the NHANES August 2021-August 2023 public-use files used in this
# project directly from the CDC/NCHS website.

files <- c(
  DEMO = "DEMO_L",  # demographics, survey design variables, exam weights
  GHB  = "GHB_L",   # glycohemoglobin (HbA1c)
  DIQ  = "DIQ_L",   # diabetes questionnaire (self-reported diagnosis)
  BMX  = "BMX_L"    # body measures (BMI, waist)
)

base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles"
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

for (f in files) {
  dest <- file.path("data/raw", paste0(f, ".xpt"))
  if (!file.exists(dest)) {
    message("Downloading ", f, " ...")
    download.file(file.path(base_url, paste0(f, ".xpt")), dest, mode = "wb", quiet = TRUE)
  }
}
message("Raw files in data/raw: ", paste(list.files("data/raw"), collapse = ", "))
