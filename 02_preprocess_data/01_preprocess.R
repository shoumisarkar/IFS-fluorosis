# Load necessary libraries
library(tidyverse)
library(readxl)
library(writexl)

# Load data
fri_scores <- read.csv("path/to/Fluorosis/RealData/Dental Exam Data/permanent_fluorosis_hypoplasia_opacity_zone.txt", sep="")
six_weeks_to_8_5_years_data <- read_excel("path/to/Fluorosis/RealData/Six Monthly Questionnaire Data/SixWeeks_8.5Yrs_DS.xlsx")
post_8_5_years_data <- read.csv("path/to/Fluorosis/RealData/Six Monthly Questionnaire Data/SixMonthQuest_DS.csv")

# Preprocess FRI Scores data
filtered_fri_scores <- fri_scores %>%
  filter(Tooth %in% 7:10) %>%
  select(1:7)

# Convert to long format
fri_scores_long <- filtered_fri_scores %>%
  pivot_longer(cols = starts_with("FRI_Score"),
               names_to = "Year",
               values_to = "FRI_Score",
               names_pattern = "FRI_Score(\\d+)") %>%
  mutate(Year = as.numeric(Year))

# Explore questionnaire data
num_cases_age_10 <- length(which(post_8_5_years_data$Qidy == 10)); num_cases_age_10
unique_ids_count <- length(unique(post_8_5_years_data$id)); unique_ids_count
# unequal values: patients are not exactly of age 10 when coming in for an `Age 10` examination;
# hence, we have to use the closest age to determine which questionnaire the data corresponds to. 
# E.g. age 10 would mean the answers are for the `Age 9` examination.

# Find closest scheduled exam age
scheduled_exam_years <- c(9, 13, 17, 23)
find_closest_value <- function(x, values) values[which.min(abs(values - x))]

# Process post-8.5-years questionnaire data
processed_questionnaire_data <- post_8_5_years_data %>%
  group_by(id) %>%
  mutate(closest_exam_age = sapply(Age_yrs, find_closest_value, values = scheduled_exam_years),
         dental_age = Age_yrs - closest_exam_age) %>%
  relocate(c(closest_exam_age, dental_age), .after = Age_yrs) %>%
  group_by(id, closest_exam_age) %>%
  filter(dental_age == min(abs(dental_age))) %>%
  select(id, closest_exam_age, dental_age, Total_mgF, SugarAddedBeverageOzPerDay, BrushingFrequencyPerDay)

# Merge FRI scores with covariates
merged_data <- fri_scores_long %>%
  inner_join(processed_questionnaire_data, by = c("SUBJECT_ID" = "id", "Year" = "closest_exam_age"))

# Process six-weeks-to-8.5-years data for time-stationary covariates
# filtering for age<=5 to be consistent with the justification discussed in Kang et al. (2022)
stationary_covariates <- six_weeks_to_8_5_years_data %>%
  filter(age_yrs <= 5) %>%
  group_by(SUBJECT_ID) %>%
  summarise(Avg_homeppm = mean(homeppm, na.rm = TRUE),
            Prop_IntQ47_DentAppt = mean(IntQ47_DentAppt == 1, na.rm = TRUE),
            Prop_IntQ48_FluorideTreatment = mean(IntQ48_FluorideTreatment == 1, na.rm = TRUE))

# Merge with covariates
final_data_with_NAs <- merged_data %>%
  inner_join(stationary_covariates, by = "SUBJECT_ID")

# Remove NAs
final_cleaned_data <- na.omit(final_data_with_NAs)

# Create Tooth_Zone_Year variable
#final_cleaned_data$Tooth_Zone_Year <- paste(final_cleaned_data$Tooth, final_cleaned_data$Zone, final_cleaned_data$Year, sep = "_")

# Assign index for repeated measures within each cluster
#unique_subject_ids <- unique(final_cleaned_data$SUBJECT_ID)
#final_cleaned_data$rm <- numeric(nrow(final_cleaned_data))

#for (i in seq_along(unique_subject_ids)) {
#  subject_rows <- final_cleaned_data$SUBJECT_ID == unique_subject_ids[i]
#  final_cleaned_data$rm[subject_rows] <- seq_along(final_cleaned_data$rm[subject_rows])
#}

final_cleaned_data = final_cleaned_data[,c(1,5,2:4,6:12)]

colnames(final_cleaned_data)[11] = "Prop_DentAppt"
colnames(final_cleaned_data)[12] = "Prop_FluorideTreatment"

# Introduce Tooth7, Tooth8, Tooth9, and Zone indicators
final_cleaned_data <- final_cleaned_data %>%
  mutate(Tooth8 = as.integer(Tooth == 8),
         Tooth9 = as.integer(Tooth == 9),
         Tooth10 = as.integer(Tooth == 10),
         ZoneM = as.integer(Zone == 'M'),
         ZoneI = as.integer(Zone == 'I'),
         ZoneE = as.integer(Zone == 'O')) #for the maxillary incisors (teeth 7-10) considered, the surface is called the incisal edge; hence the letter E in the Zone indicator.

# Remove the original Tooth and Zone columns
final_cleaned_data$Tooth <- NULL
final_cleaned_data$Zone <- NULL


data_age9 = final_cleaned_data[final_cleaned_data$Year==9,] ; data_age9$Year = NULL
data_age13 = final_cleaned_data[final_cleaned_data$Year==13,]; data_age13$Year = NULL
data_age17 = final_cleaned_data[final_cleaned_data$Year==17,]; data_age17$Year = NULL
data_age23 = final_cleaned_data[final_cleaned_data$Year==23,]; data_age23$Year = NULL

# Save preprocessed data
setwd("path/to/Fluorosis/Results/02_preprocess_data")
write_xlsx(x = data_age9, path = "preprocessed_IFS_data_age9.xlsx")
write_xlsx(x = data_age13, path = "preprocessed_IFS_data_age13.xlsx")
write_xlsx(x = data_age17, path = "preprocessed_IFS_data_age17.xlsx")
write_xlsx(x = data_age23, path = "preprocessed_IFS_data_age23.xlsx")

