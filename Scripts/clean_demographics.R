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

#Setting visit order
df<-df %>%
  mutate(visit = factor(visit, levels = c("V5", "V8", "V10", "V12")))

### Demographics table generation (Table 1)
tbl <- df %>%
  select(host_age, sex, cogresilience, visit) %>%
  tbl_summary(
    statistic = all_continuous() ~ "{mean} ({sd})",
    digits = all_continuous() ~ 1,
    missing = "no",
    label = list(
      host_age ~ "Age",
      sex ~ "Sex",
      cogresilience ~ "Cognitive resilience score",
      visit ~ "Selected Visit"
    ))

tbl

#as_flex_table(tbl) %>%
  #flextable::save_as_docx(path = "Demographics-aging.docx")


