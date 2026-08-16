# ==========================================================================
# 01_download_purpleair_outdoor.R
#
# Download and merge monthly PurpleAir outdoor files (CONUS, India, South America, Africa)
# from the EPA PurpleAir API / RSIG extracts into annual per-region RData objects,
# then combine the regional outdoor files into the pooled `outdoor` object.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #2 lines 12-257, #4 lines 508-752, #6 lines 1000-1244, #8 lines 1492-1736, #10 lines 1985-2001).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #2 (Rmd lines 12-257) ------------------------

setwd("/Users/priyanka/Downloads/PurpleAir_EPA/")
# PurpleAir Monthly Data Merger
# Memory-efficient approach for large files

library(data.table)  # Fast and memory-efficient
library(dplyr)       # For data manipulation
library(lubridate)   # For date handling

# Set up parameters
year <- "2024"  # Change to your year
month <- "01"   # Change to your month (use 01, 02, etc.)

# Define file patterns
file_patterns <- list(
  pm25_atm_a = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_atm_a.txt"),
  pm25_atm_b = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_atm_b.txt"),
  pm25_cf_1_a = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_cf_1_a.txt"),
  pm25_cf_1_b = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_cf_1_b.txt"),  
  humidity = paste0(year, "_", month, "_purpleair_conus_outdoor_humidity.txt"),
  pressure = paste0(year, "_", month, "_purpleair_conus_outdoor_pressure.txt"),
  temperature = paste0(year, "_", month, "_purpleair_conus_outdoor_temperature.txt")
)

# Function to read and process each file
read_purpleair_file <- function(filepath, data_type) {
  cat("Reading:", filepath, "\n")
  
  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read with data.table for speed and memory efficiency
  dt <- fread(filepath, 
              header = TRUE,
              stringsAsFactors = FALSE,
              showProgress = TRUE)
  
  # Define required merge columns
  required_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # Check if all required columns exist
  missing_cols <- required_cols[!required_cols %in% names(dt)]
  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns in", filepath, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Ensure Timestamp(UTC) is properly formatted
  if ("Timestamp(UTC)" %in% names(dt)) {
    dt[["Timestamp(UTC)"]] <- as.POSIXct(dt[["Timestamp(UTC)"]])
  }
  
  # Define columns that should be kept as-is (merge keys + common metadata)
  # These will be the same across files, so we keep them from the first file only
  keep_as_is <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # For the first file (pm25_atm_a), keep common metadata columns
  if (data_type == "pm25_atm_a") {
    # Keep elevation, count, and note from the first file
    if ("ELEVATION(m)" %in% names(dt)) keep_as_is <- c(keep_as_is, "ELEVATION(m)")
    if ("COUNT(-)" %in% names(dt)) keep_as_is <- c(keep_as_is, "COUNT(-)")
    if ("NOTE" %in% names(dt)) keep_as_is <- c(keep_as_is, "NOTE")
  } else {
    # For other files, remove these columns to avoid duplicates
    cols_to_remove <- c("ELEVATION(m)", "COUNT(-)", "NOTE")
    existing_cols_to_remove <- cols_to_remove[cols_to_remove %in% names(dt)]
    if (length(existing_cols_to_remove) > 0) {
      dt[, (existing_cols_to_remove) := NULL]
    }
  }
  
  # Rename measurement columns to be specific to data type
  measurement_cols <- setdiff(names(dt), keep_as_is)
  
  # Rename measurement columns with data type suffix
  if (length(measurement_cols) > 0) {
    new_names <- paste0(measurement_cols, "_", data_type)
    setnames(dt, measurement_cols, new_names)
  }
  
  cat("Loaded", nrow(dt), "rows from", filepath, "\n")
  cat("Key columns present:", paste(required_cols, collapse = ", "), "\n")
  cat("Measurement columns renamed:", length(measurement_cols), "columns\n")
  return(dt)
}

# Read all files
cat("Starting to read PurpleAir files...\n")
data_list <- list()

for (data_type in names(file_patterns)) {
  filepath <- file_patterns[[data_type]]
  data_list[[data_type]] <- read_purpleair_file(filepath, data_type)
}

# Remove NULL entries (files that couldn't be read)
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
  stop("No files could be read. Please check file paths and names.")
}

# Memory-efficient merging approach
cat("Merging files by Timestamp(UTC), LONGITUDE(deg), LATITUDE(deg), STATION(-)...\n")

# Define merge columns
merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")

# Start with the first dataset
merged_data <- data_list[[1]]

# Progressively merge with remaining datasets
if (length(data_list) > 1) {
  for (i in 2:length(data_list)) {
    cat("Merging dataset", i, "of", length(data_list), "\n")
    
    # Use data.table merge for efficiency - no suffixes needed since columns are pre-renamed
    merged_data <- merge(merged_data, 
                        data_list[[i]], 
                        by = merge_cols, 
                        all = TRUE)  # Full outer join
    
    # Clean up memory
    gc()
  }
}

# Alternative approach using dplyr (comment out data.table approach above if using this)
# merged_data <- data_list %>%
#   reduce(full_join, by = merge_cols)

# Final data summary
cat("\nMerged data summary:\n")
cat("Dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n")
cat("Timestamp range:", as.character(min(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), 
    "to", as.character(max(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), "\n")
cat("Unique stations:", length(unique(merged_data[["STATION(-)"]])), "\n")
cat("Latitude range:", round(min(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Longitude range:", round(min(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Column names:\n")
print(names(merged_data))

# Check for missing values
cat("\nMissing values by column:\n")
missing_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(missing_counts[missing_counts > 0])

# Save merged data (optional)
output_filename <- paste0(year, "_", month, "_purpleair_merged.csv")
cat("\nSaving merged data to:", output_filename, "\n")
fwrite(merged_data, output_filename, row.names = FALSE)

cat("Process completed successfully!\n")

# Optional: Create a function for processing multiple months
process_multiple_months <- function(year, months) {
  all_monthly_data <- list()
  
  # Define merge columns for the function scope
  merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  for (month in months) {
    cat("\n=== Processing", year, "-", month, "===\n")
    
    # Update file patterns for current month
    current_patterns <- list(
      pm25_atm_a = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_atm_a.txt"),
      pm25_atm_b = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_atm_b.txt"),
pm25_cf_1_a = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_cf_1_a.txt"),
      pm25_cf_1_b = paste0(year, "_", month, "_purpleair_conus_outdoor_pm25_cf_1_b.txt"),
      humidity = paste0(year, "_", month, "_purpleair_conus_outdoor_humidity.txt"),
      pressure = paste0(year, "_", month, "_purpleair_conus_outdoor_pressure.txt"),
      temperature = paste0(year, "_", month, "_purpleair_conus_outdoor_temperature.txt")
    )
    
    # Read and merge monthly files
    monthly_list <- list()
    for (data_type in names(current_patterns)) {
      filepath <- current_patterns[[data_type]]
      monthly_list[[data_type]] <- read_purpleair_file(filepath, data_type)
    }
    
    # Remove NULL entries
    monthly_list <- monthly_list[!sapply(monthly_list, is.null)]
    
    if (length(monthly_list) > 0) {
      # Merge monthly data
      monthly_merged <- monthly_list[[1]]
      if (length(monthly_list) > 1) {
        for (i in 2:length(monthly_list)) {
          monthly_merged <- merge(monthly_merged, 
                                 monthly_list[[i]], 
                                 by = merge_cols, 
                                 all = TRUE)
        }
      }
      
      all_monthly_data[[paste0(year, "_", month)]] <- monthly_merged
    }
    
    # Clean up memory after each month
    rm(monthly_list, monthly_merged)
    gc()
  }
  
  # Combine all months
  if (length(all_monthly_data) > 0) {
    cat("\nCombining all months...\n")
    final_data <- rbindlist(all_monthly_data, fill = TRUE)
    
    # Save annual merged file
    annual_filename <- paste0(year, "_purpleair_annual_merged.csv")
    fwrite(final_data, annual_filename, row.names = FALSE)
    cat("Annual data saved to:", annual_filename, "\n")
    
    return(final_data)
  }
}

#Example usage for multiple months:
 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
 annual_data_2018 <- process_multiple_months("2018", months_to_process)
 
annual_data_2019 <- process_multiple_months("2019", months_to_process)
  
annual_data_2020 <- process_multiple_months("2020", months_to_process)

annual_data_2021 <- process_multiple_months("2021", months_to_process)

annual_data_2022 <- process_multiple_months("2022", months_to_process)

annual_data_2023 <- process_multiple_months("2023", months_to_process)

annual_data_2024 <- process_multiple_months("2024", months_to_process)

 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08")
annual_data_2025 <- process_multiple_months("2025", months_to_process)

purpleair_outdoor_conus<-rbind(annual_data_2018, annual_data_2019, annual_data_2020, annual_data_2021, annual_data_2022, annual_data_2023, annual_data_2024, annual_data_2025)
save(purpleair_outdoor_conus, file="/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_conus_outdoor.RData")

# ---- source chunk #4 (Rmd lines 508-752) ------------------------

# PurpleAir Monthly Data Merger
# Memory-efficient approach for large files
setwd("/Users/priyanka/Downloads/PurpleAir_EPA/")
library(data.table)  # Fast and memory-efficient
library(dplyr)       # For data manipulation
library(lubridate)   # For date handling

# Set up parameters
year <- "2024"  # Change to your year
month <- "01"   # Change to your month (use 01, 02, etc.)

# Define file patterns
file_patterns <- list(
  pm25_atm_a = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_atm_a.txt"),
  pm25_atm_b = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_atm_b.txt"),
    pm25_cf_1_a = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_cf_1_a.txt"),
  pm25_cf_1_b = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_cf_1_b.txt"),
  humidity = paste0(year, "_", month, "_purpleair_india_outdoor_humidity.txt"),
  pressure = paste0(year, "_", month, "_purpleair_india_outdoor_pressure.txt"),
  temperature = paste0(year, "_", month, "_purpleair_india_outdoor_temperature.txt")
)

# Function to read and process each file
read_purpleair_file <- function(filepath, data_type) {
  cat("Reading:", filepath, "\n")
  
  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read with data.table for speed and memory efficiency
  dt <- fread(filepath, 
              header = TRUE,
              stringsAsFactors = FALSE,
              showProgress = TRUE)
  
  # Define required merge columns
  required_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # Check if all required columns exist
  missing_cols <- required_cols[!required_cols %in% names(dt)]
  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns in", filepath, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Ensure Timestamp(UTC) is properly formatted
  if ("Timestamp(UTC)" %in% names(dt)) {
    dt[["Timestamp(UTC)"]] <- as.POSIXct(dt[["Timestamp(UTC)"]])
  }
  
  # Define columns that should be kept as-is (merge keys + common metadata)
  # These will be the same across files, so we keep them from the first file only
  keep_as_is <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # For the first file (pm25_atm_a), keep common metadata columns
  if (data_type == "pm25_atm_a") {
    # Keep elevation, count, and note from the first file
    if ("ELEVATION(m)" %in% names(dt)) keep_as_is <- c(keep_as_is, "ELEVATION(m)")
    if ("COUNT(-)" %in% names(dt)) keep_as_is <- c(keep_as_is, "COUNT(-)")
    if ("NOTE" %in% names(dt)) keep_as_is <- c(keep_as_is, "NOTE")
  } else {
    # For other files, remove these columns to avoid duplicates
    cols_to_remove <- c("ELEVATION(m)", "COUNT(-)", "NOTE")
    existing_cols_to_remove <- cols_to_remove[cols_to_remove %in% names(dt)]
    if (length(existing_cols_to_remove) > 0) {
      dt[, (existing_cols_to_remove) := NULL]
    }
  }
  
  # Rename measurement columns to be specific to data type
  measurement_cols <- setdiff(names(dt), keep_as_is)
  
  # Rename measurement columns with data type suffix
  if (length(measurement_cols) > 0) {
    new_names <- paste0(measurement_cols, "_", data_type)
    setnames(dt, measurement_cols, new_names)
  }
  
  cat("Loaded", nrow(dt), "rows from", filepath, "\n")
  cat("Key columns present:", paste(required_cols, collapse = ", "), "\n")
  cat("Measurement columns renamed:", length(measurement_cols), "columns\n")
  return(dt)
}

# Read all files
cat("Starting to read PurpleAir files...\n")
data_list <- list()

for (data_type in names(file_patterns)) {
  filepath <- file_patterns[[data_type]]
  data_list[[data_type]] <- read_purpleair_file(filepath, data_type)
}

# Remove NULL entries (files that couldn't be read)
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
  stop("No files could be read. Please check file paths and names.")
}

# Memory-efficient merging approach
cat("Merging files by Timestamp(UTC), LONGITUDE(deg), LATITUDE(deg), STATION(-)...\n")

# Define merge columns
merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")

# Start with the first dataset
merged_data <- data_list[[1]]

# Progressively merge with remaining datasets
if (length(data_list) > 1) {
  for (i in 2:length(data_list)) {
    cat("Merging dataset", i, "of", length(data_list), "\n")
    
    # Use data.table merge for efficiency - no suffixes needed since columns are pre-renamed
    merged_data <- merge(merged_data, 
                        data_list[[i]], 
                        by = merge_cols, 
                        all = TRUE)  # Full outer join
    
    # Clean up memory
    gc()
  }
}

# Alternative approach using dplyr (comment out data.table approach above if using this)
# merged_data <- data_list %>%
#   reduce(full_join, by = merge_cols)

# Final data summary
cat("\nMerged data summary:\n")
cat("Dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n")
cat("Timestamp range:", as.character(min(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), 
    "to", as.character(max(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), "\n")
cat("Unique stations:", length(unique(merged_data[["STATION(-)"]])), "\n")
cat("Latitude range:", round(min(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Longitude range:", round(min(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Column names:\n")
print(names(merged_data))

# Check for missing values
cat("\nMissing values by column:\n")
missing_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(missing_counts[missing_counts > 0])

# Save merged data (optional)
output_filename <- paste0(year, "_", month, "_purpleair_merged.csv")
cat("\nSaving merged data to:", output_filename, "\n")
fwrite(merged_data, output_filename, row.names = FALSE)

cat("Process completed successfully!\n")

# Optional: Create a function for processing multiple months
process_multiple_months <- function(year, months) {
  all_monthly_data <- list()
  
  # Define merge columns for the function scope
  merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  for (month in months) {
    cat("\n=== Processing", year, "-", month, "===\n")
    
    # Update file patterns for current month
    current_patterns <- list(
      pm25_atm_a = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_atm_a.txt"),
      pm25_atm_b = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_atm_b.txt"),
      pm25_cf_1_a = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_cf_1_a.txt"),
      pm25_cf_1_b = paste0(year, "_", month, "_purpleair_india_outdoor_pm25_cf_1_b.txt"),
      humidity = paste0(year, "_", month, "_purpleair_india_outdoor_humidity.txt"),
      pressure = paste0(year, "_", month, "_purpleair_india_outdoor_pressure.txt"),
      temperature = paste0(year, "_", month, "_purpleair_india_outdoor_temperature.txt")
    )
    
    # Read and merge monthly files
    monthly_list <- list()
    for (data_type in names(current_patterns)) {
      filepath <- current_patterns[[data_type]]
      monthly_list[[data_type]] <- read_purpleair_file(filepath, data_type)
    }
    
    # Remove NULL entries
    monthly_list <- monthly_list[!sapply(monthly_list, is.null)]
    
    if (length(monthly_list) > 0) {
      # Merge monthly data
      monthly_merged <- monthly_list[[1]]
      if (length(monthly_list) > 1) {
        for (i in 2:length(monthly_list)) {
          monthly_merged <- merge(monthly_merged, 
                                 monthly_list[[i]], 
                                 by = merge_cols, 
                                 all = TRUE)
        }
      }
      
      all_monthly_data[[paste0(year, "_", month)]] <- monthly_merged
    }
    
    # Clean up memory after each month
    rm(monthly_list, monthly_merged)
    gc()
  }
  
  # Combine all months
  if (length(all_monthly_data) > 0) {
    cat("\nCombining all months...\n")
    final_data <- rbindlist(all_monthly_data, fill = TRUE)
    
    # Save annual merged file
    annual_filename <- paste0(year, "_purpleair_annual_merged.csv")
    fwrite(final_data, annual_filename, row.names = FALSE)
    cat("Annual data saved to:", annual_filename, "\n")
    
    return(final_data)
  }
}

#Example usage for multiple months:
 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
 annual_data_2018 <- process_multiple_months("2018", months_to_process)
 
annual_data_2019 <- process_multiple_months("2019", months_to_process)
  
annual_data_2020 <- process_multiple_months("2020", months_to_process)

annual_data_2021 <- process_multiple_months("2021", months_to_process)

annual_data_2022 <- process_multiple_months("2022", months_to_process)

annual_data_2023 <- process_multiple_months("2023", months_to_process)

annual_data_2024 <- process_multiple_months("2024", months_to_process)

 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08")
annual_data_2025 <- process_multiple_months("2025", months_to_process)

purpleair_outdoor_india<-rbind(annual_data_2018, annual_data_2019, annual_data_2020, annual_data_2021, annual_data_2022, annual_data_2023, annual_data_2024, annual_data_2025)
save(purpleair_outdoor_india, file="/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_india_outdoor.RData")

# ---- source chunk #6 (Rmd lines 1000-1244) ------------------------

# PurpleAir Monthly Data Merger
# Memory-efficient approach for large files
setwd("/Users/priyanka/Downloads/PurpleAir_EPA/")
library(data.table)  # Fast and memory-efficient
library(dplyr)       # For data manipulation
library(lubridate)   # For date handling

# Set up parameters
year <- "2024"  # Change to your year
month <- "01"   # Change to your month (use 01, 02, etc.)

# Define file patterns
file_patterns <- list(
  pm25_atm_a = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_atm_a.txt"),
  pm25_atm_b = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_atm_b.txt"),
    pm25_cf_1_a = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_cf_1_a.txt"),
  pm25_cf_1_b = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_cf_1_b.txt"),
  humidity = paste0(year, "_", month, "_purpleair_south_america_outdoor_humidity.txt"),
  pressure = paste0(year, "_", month, "_purpleair_south_america_outdoor_pressure.txt"),
  temperature = paste0(year, "_", month, "_purpleair_south_america_outdoor_temperature.txt")
)

# Function to read and process each file
read_purpleair_file <- function(filepath, data_type) {
  cat("Reading:", filepath, "\n")
  
  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read with data.table for speed and memory efficiency
  dt <- fread(filepath, 
              header = TRUE,
              stringsAsFactors = FALSE,
              showProgress = TRUE)
  
  # Define required merge columns
  required_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # Check if all required columns exist
  missing_cols <- required_cols[!required_cols %in% names(dt)]
  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns in", filepath, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Ensure Timestamp(UTC) is properly formatted
  if ("Timestamp(UTC)" %in% names(dt)) {
    dt[["Timestamp(UTC)"]] <- as.POSIXct(dt[["Timestamp(UTC)"]])
  }
  
  # Define columns that should be kept as-is (merge keys + common metadata)
  # These will be the same across files, so we keep them from the first file only
  keep_as_is <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # For the first file (pm25_atm_a), keep common metadata columns
  if (data_type == "pm25_atm_a") {
    # Keep elevation, count, and note from the first file
    if ("ELEVATION(m)" %in% names(dt)) keep_as_is <- c(keep_as_is, "ELEVATION(m)")
    if ("COUNT(-)" %in% names(dt)) keep_as_is <- c(keep_as_is, "COUNT(-)")
    if ("NOTE" %in% names(dt)) keep_as_is <- c(keep_as_is, "NOTE")
  } else {
    # For other files, remove these columns to avoid duplicates
    cols_to_remove <- c("ELEVATION(m)", "COUNT(-)", "NOTE")
    existing_cols_to_remove <- cols_to_remove[cols_to_remove %in% names(dt)]
    if (length(existing_cols_to_remove) > 0) {
      dt[, (existing_cols_to_remove) := NULL]
    }
  }
  
  # Rename measurement columns to be specific to data type
  measurement_cols <- setdiff(names(dt), keep_as_is)
  
  # Rename measurement columns with data type suffix
  if (length(measurement_cols) > 0) {
    new_names <- paste0(measurement_cols, "_", data_type)
    setnames(dt, measurement_cols, new_names)
  }
  
  cat("Loaded", nrow(dt), "rows from", filepath, "\n")
  cat("Key columns present:", paste(required_cols, collapse = ", "), "\n")
  cat("Measurement columns renamed:", length(measurement_cols), "columns\n")
  return(dt)
}

# Read all files
cat("Starting to read PurpleAir files...\n")
data_list <- list()

for (data_type in names(file_patterns)) {
  filepath <- file_patterns[[data_type]]
  data_list[[data_type]] <- read_purpleair_file(filepath, data_type)
}

# Remove NULL entries (files that couldn't be read)
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
  stop("No files could be read. Please check file paths and names.")
}

# Memory-efficient merging approach
cat("Merging files by Timestamp(UTC), LONGITUDE(deg), LATITUDE(deg), STATION(-)...\n")

# Define merge columns
merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")

# Start with the first dataset
merged_data <- data_list[[1]]

# Progressively merge with remaining datasets
if (length(data_list) > 1) {
  for (i in 2:length(data_list)) {
    cat("Merging dataset", i, "of", length(data_list), "\n")
    
    # Use data.table merge for efficiency - no suffixes needed since columns are pre-renamed
    merged_data <- merge(merged_data, 
                        data_list[[i]], 
                        by = merge_cols, 
                        all = TRUE)  # Full outer join
    
    # Clean up memory
    gc()
  }
}

# Alternative approach using dplyr (comment out data.table approach above if using this)
# merged_data <- data_list %>%
#   reduce(full_join, by = merge_cols)

# Final data summary
cat("\nMerged data summary:\n")
cat("Dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n")
cat("Timestamp range:", as.character(min(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), 
    "to", as.character(max(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), "\n")
cat("Unique stations:", length(unique(merged_data[["STATION(-)"]])), "\n")
cat("Latitude range:", round(min(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Longitude range:", round(min(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Column names:\n")
print(names(merged_data))

# Check for missing values
cat("\nMissing values by column:\n")
missing_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(missing_counts[missing_counts > 0])

# Save merged data (optional)
output_filename <- paste0(year, "_", month, "_purpleair_merged.csv")
cat("\nSaving merged data to:", output_filename, "\n")
fwrite(merged_data, output_filename, row.names = FALSE)

cat("Process completed successfully!\n")

# Optional: Create a function for processing multiple months
process_multiple_months <- function(year, months) {
  all_monthly_data <- list()
  
  # Define merge columns for the function scope
  merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  for (month in months) {
    cat("\n=== Processing", year, "-", month, "===\n")
    
    # Update file patterns for current month
    current_patterns <- list(
      pm25_atm_a = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_atm_a.txt"),
      pm25_atm_b = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_atm_b.txt"),
      pm25_cf_1_a = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_cf_1_a.txt"),
      pm25_cf_1_b = paste0(year, "_", month, "_purpleair_south_america_outdoor_pm25_cf_1_b.txt"),
      humidity = paste0(year, "_", month, "_purpleair_south_america_outdoor_humidity.txt"),
      pressure = paste0(year, "_", month, "_purpleair_south_america_outdoor_pressure.txt"),
      temperature = paste0(year, "_", month, "_purpleair_south_america_outdoor_temperature.txt")
    )
    
    # Read and merge monthly files
    monthly_list <- list()
    for (data_type in names(current_patterns)) {
      filepath <- current_patterns[[data_type]]
      monthly_list[[data_type]] <- read_purpleair_file(filepath, data_type)
    }
    
    # Remove NULL entries
    monthly_list <- monthly_list[!sapply(monthly_list, is.null)]
    
    if (length(monthly_list) > 0) {
      # Merge monthly data
      monthly_merged <- monthly_list[[1]]
      if (length(monthly_list) > 1) {
        for (i in 2:length(monthly_list)) {
          monthly_merged <- merge(monthly_merged, 
                                 monthly_list[[i]], 
                                 by = merge_cols, 
                                 all = TRUE)
        }
      }
      
      all_monthly_data[[paste0(year, "_", month)]] <- monthly_merged
    }
    
    # Clean up memory after each month
    rm(monthly_list, monthly_merged)
    gc()
  }
  
  # Combine all months
  if (length(all_monthly_data) > 0) {
    cat("\nCombining all months...\n")
    final_data <- rbindlist(all_monthly_data, fill = TRUE)
    
    # Save annual merged file
    annual_filename <- paste0(year, "_purpleair_annual_merged.csv")
    fwrite(final_data, annual_filename, row.names = FALSE)
    cat("Annual data saved to:", annual_filename, "\n")
    
    return(final_data)
  }
}

#Example usage for multiple months:
 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
 annual_data_2018 <- process_multiple_months("2018", months_to_process)
 
annual_data_2019 <- process_multiple_months("2019", months_to_process)
  
annual_data_2020 <- process_multiple_months("2020", months_to_process)

annual_data_2021 <- process_multiple_months("2021", months_to_process)

annual_data_2022 <- process_multiple_months("2022", months_to_process)

annual_data_2023 <- process_multiple_months("2023", months_to_process)

annual_data_2024 <- process_multiple_months("2024", months_to_process)

 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08")
annual_data_2025 <- process_multiple_months("2025", months_to_process)

purpleair_outdoor_south_america<-rbind(annual_data_2018, annual_data_2019, annual_data_2020, annual_data_2021, annual_data_2022, annual_data_2023, annual_data_2024, annual_data_2025)
save(purpleair_outdoor_south_america, file="/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_south_america_outdoor.RData")

# ---- source chunk #8 (Rmd lines 1492-1736) ------------------------

# PurpleAir Monthly Data Merger
# Memory-efficient approach for large files
setwd("/Users/priyanka/Downloads/PurpleAir_EPA/")
library(data.table)  # Fast and memory-efficient
library(dplyr)       # For data manipulation
library(lubridate)   # For date handling

# Set up parameters
year <- "2024"  # Change to your year
month <- "01"   # Change to your month (use 01, 02, etc.)

# Define file patterns
file_patterns <- list(
  pm25_atm_a = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_atm_a.txt"),
  pm25_atm_b = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_atm_b.txt"),
    pm25_cf_1_a = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_cf_1_a.txt"),
  pm25_cf_1_b = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_cf_1_b.txt"),
  humidity = paste0(year, "_", month, "_purpleair_africa_outdoor_humidity.txt"),
  pressure = paste0(year, "_", month, "_purpleair_africa_outdoor_pressure.txt"),
  temperature = paste0(year, "_", month, "_purpleair_africa_outdoor_temperature.txt")
)

# Function to read and process each file
read_purpleair_file <- function(filepath, data_type) {
  cat("Reading:", filepath, "\n")
  
  # Check if file exists
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(NULL)
  }
  
  # Read with data.table for speed and memory efficiency
  dt <- fread(filepath, 
              header = TRUE,
              stringsAsFactors = FALSE,
              showProgress = TRUE)
  
  # Define required merge columns
  required_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # Check if all required columns exist
  missing_cols <- required_cols[!required_cols %in% names(dt)]
  if (length(missing_cols) > 0) {
    warning(paste("Missing required columns in", filepath, ":", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Ensure Timestamp(UTC) is properly formatted
  if ("Timestamp(UTC)" %in% names(dt)) {
    dt[["Timestamp(UTC)"]] <- as.POSIXct(dt[["Timestamp(UTC)"]])
  }
  
  # Define columns that should be kept as-is (merge keys + common metadata)
  # These will be the same across files, so we keep them from the first file only
  keep_as_is <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  # For the first file (pm25_atm_a), keep common metadata columns
  if (data_type == "pm25_atm_a") {
    # Keep elevation, count, and note from the first file
    if ("ELEVATION(m)" %in% names(dt)) keep_as_is <- c(keep_as_is, "ELEVATION(m)")
    if ("COUNT(-)" %in% names(dt)) keep_as_is <- c(keep_as_is, "COUNT(-)")
    if ("NOTE" %in% names(dt)) keep_as_is <- c(keep_as_is, "NOTE")
  } else {
    # For other files, remove these columns to avoid duplicates
    cols_to_remove <- c("ELEVATION(m)", "COUNT(-)", "NOTE")
    existing_cols_to_remove <- cols_to_remove[cols_to_remove %in% names(dt)]
    if (length(existing_cols_to_remove) > 0) {
      dt[, (existing_cols_to_remove) := NULL]
    }
  }
  
  # Rename measurement columns to be specific to data type
  measurement_cols <- setdiff(names(dt), keep_as_is)
  
  # Rename measurement columns with data type suffix
  if (length(measurement_cols) > 0) {
    new_names <- paste0(measurement_cols, "_", data_type)
    setnames(dt, measurement_cols, new_names)
  }
  
  cat("Loaded", nrow(dt), "rows from", filepath, "\n")
  cat("Key columns present:", paste(required_cols, collapse = ", "), "\n")
  cat("Measurement columns renamed:", length(measurement_cols), "columns\n")
  return(dt)
}

# Read all files
cat("Starting to read PurpleAir files...\n")
data_list <- list()

for (data_type in names(file_patterns)) {
  filepath <- file_patterns[[data_type]]
  data_list[[data_type]] <- read_purpleair_file(filepath, data_type)
}

# Remove NULL entries (files that couldn't be read)
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
  stop("No files could be read. Please check file paths and names.")
}

# Memory-efficient merging approach
cat("Merging files by Timestamp(UTC), LONGITUDE(deg), LATITUDE(deg), STATION(-)...\n")

# Define merge columns
merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")

# Start with the first dataset
merged_data <- data_list[[1]]

# Progressively merge with remaining datasets
if (length(data_list) > 1) {
  for (i in 2:length(data_list)) {
    cat("Merging dataset", i, "of", length(data_list), "\n")
    
    # Use data.table merge for efficiency - no suffixes needed since columns are pre-renamed
    merged_data <- merge(merged_data, 
                        data_list[[i]], 
                        by = merge_cols, 
                        all = TRUE)  # Full outer join
    
    # Clean up memory
    gc()
  }
}

# Alternative approach using dplyr (comment out data.table approach above if using this)
# merged_data <- data_list %>%
#   reduce(full_join, by = merge_cols)

# Final data summary
cat("\nMerged data summary:\n")
cat("Dimensions:", nrow(merged_data), "rows x", ncol(merged_data), "columns\n")
cat("Timestamp range:", as.character(min(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), 
    "to", as.character(max(merged_data[["Timestamp(UTC)"]], na.rm = TRUE)), "\n")
cat("Unique stations:", length(unique(merged_data[["STATION(-)"]])), "\n")
cat("Latitude range:", round(min(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LATITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Longitude range:", round(min(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), 
    "to", round(max(merged_data[["LONGITUDE(deg)"]], na.rm = TRUE), 3), "\n")
cat("Column names:\n")
print(names(merged_data))

# Check for missing values
cat("\nMissing values by column:\n")
missing_counts <- sapply(merged_data, function(x) sum(is.na(x)))
print(missing_counts[missing_counts > 0])

# Save merged data (optional)
output_filename <- paste0(year, "_", month, "_purpleair_merged.csv")
cat("\nSaving merged data to:", output_filename, "\n")
fwrite(merged_data, output_filename, row.names = FALSE)

cat("Process completed successfully!\n")

# Optional: Create a function for processing multiple months
process_multiple_months <- function(year, months) {
  all_monthly_data <- list()
  
  # Define merge columns for the function scope
  merge_cols <- c("Timestamp(UTC)", "LONGITUDE(deg)", "LATITUDE(deg)", "STATION(-)")
  
  for (month in months) {
    cat("\n=== Processing", year, "-", month, "===\n")
    
    # Update file patterns for current month
    current_patterns <- list(
      pm25_atm_a = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_atm_a.txt"),
      pm25_atm_b = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_atm_b.txt"),
      pm25_cf_1_a = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_cf_1_a.txt"),
      pm25_cf_1_b = paste0(year, "_", month, "_purpleair_africa_outdoor_pm25_cf_1_b.txt"),
      humidity = paste0(year, "_", month, "_purpleair_africa_outdoor_humidity.txt"),
      pressure = paste0(year, "_", month, "_purpleair_africa_outdoor_pressure.txt"),
      temperature = paste0(year, "_", month, "_purpleair_africa_outdoor_temperature.txt")
    )
    
    # Read and merge monthly files
    monthly_list <- list()
    for (data_type in names(current_patterns)) {
      filepath <- current_patterns[[data_type]]
      monthly_list[[data_type]] <- read_purpleair_file(filepath, data_type)
    }
    
    # Remove NULL entries
    monthly_list <- monthly_list[!sapply(monthly_list, is.null)]
    
    if (length(monthly_list) > 0) {
      # Merge monthly data
      monthly_merged <- monthly_list[[1]]
      if (length(monthly_list) > 1) {
        for (i in 2:length(monthly_list)) {
          monthly_merged <- merge(monthly_merged, 
                                 monthly_list[[i]], 
                                 by = merge_cols, 
                                 all = TRUE)
        }
      }
      
      all_monthly_data[[paste0(year, "_", month)]] <- monthly_merged
    }
    
    # Clean up memory after each month
    rm(monthly_list, monthly_merged)
    gc()
  }
  
  # Combine all months
  if (length(all_monthly_data) > 0) {
    cat("\nCombining all months...\n")
    final_data <- rbindlist(all_monthly_data, fill = TRUE)
    
    # Save annual merged file
    annual_filename <- paste0(year, "_purpleair_annual_merged.csv")
    fwrite(final_data, annual_filename, row.names = FALSE)
    cat("Annual data saved to:", annual_filename, "\n")
    
    return(final_data)
  }
}

#Example usage for multiple months:
 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
 annual_data_2018 <- process_multiple_months("2018", months_to_process)
 
annual_data_2019 <- process_multiple_months("2019", months_to_process)
  
annual_data_2020 <- process_multiple_months("2020", months_to_process)

annual_data_2021 <- process_multiple_months("2021", months_to_process)

annual_data_2022 <- process_multiple_months("2022", months_to_process)

annual_data_2023 <- process_multiple_months("2023", months_to_process)

annual_data_2024 <- process_multiple_months("2024", months_to_process)

 months_to_process <- c("01", "02", "03", "04", "05", "06", "07", "08")
annual_data_2025 <- process_multiple_months("2025", months_to_process)

purpleair_outdoor_africa<-rbind(annual_data_2018, annual_data_2019, annual_data_2020, annual_data_2021, annual_data_2022, annual_data_2023, annual_data_2024, annual_data_2025)
save(purpleair_outdoor_africa, file="/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_africa_outdoor.RData")

# ---- source chunk #10 (Rmd lines 1985-2001) ------------------------

load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_india_outdoor.RData")
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_conus_outdoor.RData")
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_africa_outdoor.RData")
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_south_america_outdoor.RData")

purpleair_outdoor_conus$region<-"CONUS"
purpleair_outdoor_india$region<-"India"
purpleair_outdoor_africa$region<-"Africa"
purpleair_outdoor_south_america$region<-"SouthAmerica"

outdoor<-rbind(purpleair_outdoor_conus, purpleair_outdoor_india, purpleair_outdoor_south_america, purpleair_outdoor_africa)

outdoor<-outdoor %>%  dplyr::mutate(datetime = floor_date(`Timestamp(UTC)`, "hour"))

save(outdoor, file="/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor.RData")
