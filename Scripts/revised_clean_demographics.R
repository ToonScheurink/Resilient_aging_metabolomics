### Packages
library(dplyr)
library(gtsummary)
library(gt)         
library(broom)
library(forcats)
library(flextable)

### Data
df2 <- readRDS("clean_data_wellcome.rds")

features <-(4:22664)

df1<- df2 %>%
  select(-features)

df1$host_age<-as.numeric(df1$host_age)
df1$host_body_mass_index<-as.numeric(df1$host_body_mass_index)
df1$physresilience_pctl<-as.numeric(df1$physresilience_pctl)
df1$physresilience<-as.numeric(df1$physresilience)
df1$cogresilience<-as.numeric(df1$cogresilience)
df1$cogresilience_pctl<-as.numeric(df1$cogresilience_pctl)

#Setting age target
target_age <- 73.7

#Age harmonization
df <- df1 %>%
  group_by(record_ID) %>%
  slice_min(order_by = abs(host_age - target_age), with_ties = FALSE) %>%
  ungroup()

df <- df %>%
  mutate(merge_id = paste(record_ID, visit, sep = "_"))

#Loading additional demographics
lab <- read.csv("RBS_LAB_BLD.CSV")
lab$merge_id <- paste0(lab$SUBJID, "_V", lab$VISIT)

phys <-read.csv("RBS_PHYSCHAR.CSV")
phys$merge_id <- paste0(phys$SUBJID, "_V", phys$VISIT)


#Merging with df
df <- df %>%
  left_join(lab, by = "merge_id") %>%
  left_join(phys, by = "merge_id")

#Setting visit order
df<-df %>%
  mutate(visit = factor(visit, levels = c("V5", "V8", "V10", "V12")))


#Selecting MAP 
# MAP for measurement 1
df$MAP1 <- (df$SBP1 + 2 * df$DBP1) / 3
# MAP for measurement 2
df$MAP2 <- (df$SBP2 + 2 * df$DBP2) / 3
# Average MAP across both measurements, ignoring missing values
df$MAP_mean <- rowMeans(df[, c("MAP1", "MAP2")], na.rm = TRUE)




### Demographics table generation (Table 1)
tbl <- df %>%
  select(host_age, sex, cogresilience, MAP_mean, WAISTMIN_CM, FATPERC, GLUCOSE, CHOLESTEROL, TRIGLYCERIDES, HDL, LDL, CREATININE, visit) %>%
  tbl_summary(
    statistic = all_continuous() ~ "{mean} ({sd})",
    digits = all_continuous() ~ 1,
    missing = "no",
    label = list(
      host_age ~ "Age",
      sex ~ "Sex",
      MAP_mean ~ "MAP (mm/Hg)",
     WAISTMIN_CM ~ "Waist circumference (cm)",
     FATPERC ~ "Total body fat (%)*",
      cogresilience ~ "Cognitive resilience score",
     GLUCOSE~ "Glucose (mg/dL)*",
     CHOLESTEROL ~ "Total cholesterol (mg/dL)*",
     TRIGLYCERIDES ~ "Triglycerides (mg/dL)*",
     HDL ~ "HDL (mg/dL)*",
     LDL ~ "LDL (mg/dL)*",
     CREATININE ~ "Createnine (mg/dL)*",
      visit ~ "Selected Visit"
    ))

tbl

#as_flex_table(tbl) %>% flextable::save_as_docx(path = "Demographics-aging.docx")


