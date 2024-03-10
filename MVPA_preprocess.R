rm(list=ls())
cat('\014')

library(data.table)
library(survival)
library(prodlim)
library(ggplot2)
library(reshape2)
library(plyr)
library(Cairo)
library(stringr)
library(Epi)
library(dplyr)

#Data load
data <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/mvpa_VTE5.csv", header = TRUE)
summary(data)
nrow(data) #94482
data$eid <- as.character(data$eid)


#Withdraw
withdrawals <- read.table("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/withdraw100739_15_20231013.txt", header = TRUE)
withdrawals$A1 <- as.character(withdrawals$A1)

data <- data[!(data$eid %in% withdrawals$A1), ]
nrow(data) #94482-1=94481

# Remove blank overall sleep
data$Sleep_overall_average <- as.numeric(data$Sleep_overall_average)
data <- data[!is.na(data$Sleep_overall_average),]
nrow(data) #Still 94481

# Remove blank sleep daily
data$Sleep_day_average <- as.character(data$Sleep_day_average)
data <- data[data$Sleep_day_average != '',]
nrow(data) #94481-8 = 94473


#sb_daily: sleep behavior daily
split_cols <- strsplit(as.character(data$Sleep_day_average), ",")
split_cols <- do.call(rbind, split_cols)
split_cols <- as.data.frame(split_cols)
colnames(split_cols) <- c("sleep1", "sleep2", "sleep3", "sleep4", "sleep5", "sleep6", "sleep7")
data <- cbind(data, split_cols)

#Remove bland daily sb
data$sleep1 <- as.numeric(data$sleep1)
data$sleep2 <- as.numeric(data$sleep2)
data$sleep3 <- as.numeric(data$sleep3)
data$sleep4 <- as.numeric(data$sleep4)
data$sleep5 <- as.numeric(data$sleep5)
data$sleep6 <- as.numeric(data$sleep6)
data$sleep7 <- as.numeric(data$sleep7)

data <- data[!(data$sleep1 <= 0.21 | data$sleep2 <= 0.21 | data$sleep3 <= 0.21 | data$sleep4 <= 0.21 | data$sleep5 <= 0.21 | data$sleep6 <= 0.21 | data$sleep7 <= 0.21), ]
nrow(data) #  90615


# Remove blank overall mvpa
data$MVPA_overall_average <- as.numeric(data$MVPA_overall_average)
data <- data[!is.na(data$MVPA_overall_average),]
nrow(data) #still 90615

# Remove blank mvpa daily
data$MVPA_day_average <- as.character(data$MVPA_day_average)
data <- data[data$MVPA_day_average != '',]
nrow(data) #Still 90615


#mvpa_daily: MVPA daily
split_cols <- strsplit(as.character(data$MVPA_day_average), ",")
split_cols <- do.call(rbind, split_cols)
split_cols <- as.data.frame(split_cols)
colnames(split_cols) <- c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")
data <- cbind(data, split_cols)

#Remove bland daily mvpa
data$mvpa1 <- as.numeric(data$mvpa1)
data$mvpa2 <- as.numeric(data$mvpa2)
data$mvpa3 <- as.numeric(data$mvpa3)
data$mvpa4 <- as.numeric(data$mvpa4)
data$mvpa5 <- as.numeric(data$mvpa5)
data$mvpa6 <- as.numeric(data$mvpa6)
data$mvpa7 <- as.numeric(data$mvpa7)
data <- data[!c(is.na(data$mvpa1) | is.na(data$mvpa2) | is.na(data$mvpa3) | is.na(data$mvpa4) | is.na(data$mvpa5) | is.na(data$mvpa6) | is.na(data$mvpa7)),]
nrow(data) # still 90615



# Age
data$Starttime_wear <- as.Date(data$Starttime_wear)
data$Birth_year <- as.integer(data$Birth_year)
data$Birth_month <- factor(data$Birth_month, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"))
month_mapping <- c(January = 1, February = 2, March = 3, April = 4, May = 5, June = 6, July = 7, August = 8, September = 9, October = 10, November = 11, December = 12)
data$Birth_month <- month_mapping[data$Birth_month]
data$Age <- as.numeric(format(data$Starttime_wear, "%Y")) - data$Birth_year
data$Age[data$Birth_month > as.integer(format(data$Starttime_wear, "%m"))] <- data$Age - 1
summary(data$Age)

data <- data[data$Age >= 60, ]
nrow(data) # 90615-33620 =56995


#PRS numeric
data$Standard_PRS <- as.numeric(data$Standard_PRS)
summary(data$Standard_PRS)
data <- data[!is.na(data$Standard_PRS), ]
nrow(data) # 56995- 1400=55595


#Sex
data$Sex <- ifelse(data$Sex == "Female", 2, 1)
data$Sex <-  as.factor(data$Sex)

#Ethnicity
data$Ethnicity <-  as.factor(data$Ethnicity)

data$Ethnicity <-
  plyr::mapvalues(
    data$Ethnicity,
    c(
      "British" ,
      "Any other white background",
      "Irish",
      "White and Asian",
      "Other ethnic group",
      "Caribbean",
      "Chinese",
      "Indian",
      "Pakistani",
      "White and Black African",
      "Any other mixed background",
      "African",
      "White and Black Caribbean",
      "Prefer not to answer",
      "White",
      "Do not know" ,
      "Any other Black background",
      "Any other Asian background",
      "" ,
      "Bangladeshi",
      "Mixed",
      "Asian or Asian British",
      "Black or Black British"
    ),
    c(
      "White" ,
      "White",
      "White",
      "Mixed_and_other",
      "Mixed_and_other",
      "Black",
      "Mixed_and_other",
      "Asian",
      "Asian",
      "Mixed_and_other",
      "Mixed_and_other",
      "Black",
      "Mixed_and_other",
      NA,
      "White",
      NA ,
      "Black",
      "Asian",
      NA ,
      "Asian",
      "Mixed_and_other",
      "Asian",
      "Black"
    )
  )
summary(data$Ethnicity)
data <- data[!is.na(data$Ethnicity), ]
nrow(data) # 55595-204=55391

#Education
data$Education <-  as.character(data$Education)
unique(data$Education)

data <- data[(data$Education != "Prefer not to answer" & data$Education != ""), ]

data$Education[grepl("(?i)college", data$Education)] <- "college or higher"
data$Education[!grepl("(?i)college", data$Education)] <- "others"
summary(as.factor(data$Education)) 
nrow(data)# still 55391-385 =55006


#Townsend deprivation index
data$Townsend_deprivation_index <- as.numeric(data$Townsend_deprivation_index)
summary(data$Townsend_deprivation_index)
data <- data[!is.na(data$Townsend_deprivation_index), ]
nrow(data) # 55006-57=54949


#Household income
data$Average_totalhousehold_income_beforetax <- as.character(data$Average_totalhousehold_income_beforetax)
data$Average_totalhousehold_income_beforetax <- gsub("Less than 18,000", "<18000", data$Average_totalhousehold_income_beforetax)
data$Average_totalhousehold_income_beforetax <- gsub("18,000 to 30,999", "18000 to 30999", data$Average_totalhousehold_income_beforetax)
data$Average_totalhousehold_income_beforetax <- gsub("31,000 to 51,999", "31000 to 51999", data$Average_totalhousehold_income_beforetax)
data$Average_totalhousehold_income_beforetax <- gsub("52,000 to 100,000|Greater than 100,000", ">52000", data$Average_totalhousehold_income_beforetax)
data$Average_totalhousehold_income_beforetax[(data$Average_totalhousehold_income_beforetax %in% c("Prefer not to answer", "Do not know", ""))] <- NA
summary(as.factor(data$Average_totalhousehold_income_beforetax))
data <- data[!is.na(data$Average_totalhousehold_income_beforetax), ]
nrow(data) #54949-5782=49167

#Employment
sum(grepl("employment", data$Employment_status, ignore.case = TRUE))
data$Employment_status <- ifelse(grepl("(?i)employment", data$Employment_status), "employed", "others")
summary(as.factor(data$Employment_status))
nrow(data) # Still 49167

#Smoke
data$Smoking_status[!(data$Smoking_status %in% c("Current", "Never", "Previous"))] <- NA
data <- data[!is.na(data$Smoking_status), ]

sum(grepl("current", data$Smoking_status, ignore.case = TRUE))
data$Smoking_status <- ifelse(grepl("(?i)current", data$Smoking_status), "current smoking", "others")
summary(as.factor(data$Smoking_status))
nrow(data) # 49167-98=49069


#Alcohol
data$Alcohol_frequency <-
  plyr::mapvalues(
    data$Alcohol_frequency,
    from = c(
      "Never",
      "Three or four times a week",
      "Daily or almost daily",
      "Once or twice a week",
      "One to three times a month",
      "Special occasions only",
      "Prefer not to answer",
      "",
      "Do not know"
    ),
    to = c(
      "Never",
      "3+ times per week",
      "3+ times per week",
      "< 3 times per week",
      "< 3 times per week",
      "< 3 times per week",
      NA,
      NA,
      NA
    )
  )
summary(as.factor(data$Alcohol_frequency))
data <- data[!is.na(data$Alcohol_frequency), ]
nrow(data) #49069-10 = 49059

#Accodring to SB-Dementia


#Fruit: >= 3 servings/day of fruit
unique(data$Fresh_fruit_intake)
summary(as.factor(data$Fresh_fruit_intake))
data$Fresh_fruit_intake <-
  plyr::mapvalues(
    data$Fresh_fruit_intake,
    from = c(
      "Less than one",
      "Do not know",
      "Prefer not to answer"
    ),
    to = c(
      "0.5",
      NA,
      NA
    )
  )

data <- data[!is.na(data$Fresh_fruit_intake), ]
nrow(data) #49059-47=49012

summary(as.factor(data$Dried_fruit_intake))
data$Dried_fruit_intake <-
  plyr::mapvalues(
    data$Dried_fruit_intake,
    from = c(
      "Less than one",
      "Do not know",
      "Prefer not to answer"
    ),
    to = c(
      "0.5",
      NA,
      NA
    )
  )
data <- data[!is.na(data$Dried_fruit_intake), ]
nrow(data) #49012-222=48790

#Fruit score
data$Fresh_fruit_intake <- as.numeric(data$Fresh_fruit_intake)
data$Dried_fruit_intake <- as.numeric(data$Dried_fruit_intake)

data$Fruit_score <- ifelse( (data$Fresh_fruit_intake/1 + data$Dried_fruit_intake/5) >= 3, 1, 0)
summary(as.factor(data$Fruit_score))

#Vegetables: >= 3 servings/day of vegetables
data$Raw_vegetable_intake <-
  plyr::mapvalues(
    data$Raw_vegetable_intake,
    from = c(
      "Less than one",
      "Do not know",
      "Prefer not to answer"
    ),
    to = c(
      "0.5",
      NA,
      NA
    )
  )
summary(as.factor(data$Raw_vegetable_intake))
data <- data[!is.na(data$Raw_vegetable_intake), ]
nrow(data) #48790- 177 = 48613

data$Cooked_vegerable_intake <-
  plyr::mapvalues(
    data$Cooked_vegerable_intake,
    from = c(
      "Less than one",
      "Do not know",
      "Prefer not to answer"
    ),
    to = c(
      "0.5",
      NA,
      NA
    )
  )
summary(as.factor(data$Cooked_vegerable_intake))
data <- data[!is.na(data$Cooked_vegerable_intake), ]
nrow(data) #48613 - 86 = 48527

#Vegetable score
data$Raw_vegetable_intake <- as.numeric(data$Raw_vegetable_intake)
data$Cooked_vegerable_intake <- as.numeric(data$Cooked_vegerable_intake)

data$Vegetable_score <- ifelse( (data$Raw_vegetable_intake/3 + data$Cooked_vegerable_intake/3) >= 3, 1, 0)
summary(as.factor(data$Vegetable_score))


#Fish: >= 2 servings/week of fish
data$Oily_fish_intake <-
  plyr::mapvalues(
    data$Oily_fish_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(
      "0",
      "0.5",
      "1",
      "3",
      NA,
      "5.5",
      NA,
      NA,
      "7"
    )
  )
summary(as.factor(data$Oily_fish_intake))
data <- data[!is.na(data$Oily_fish_intake), ]
nrow(data) #48527 - 63 = 48464

data$Non_oily_fish_intake <-
  plyr::mapvalues(
    data$Non_oily_fish_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(
      "0",
      "0.5",
      "1",
      "3",
      NA,
      "5.5",
      NA,
      NA,
      "7"
    )
  )
summary(as.factor(data$Non_oily_fish_intake))
data <- data[!is.na(data$Non_oily_fish_intake), ]
nrow(data) #48464 - 41 = 48423


#Fish score
data$Oily_fish_intake <- as.numeric(data$Oily_fish_intake)
data$Non_oily_fish_intake <- as.numeric(data$Non_oily_fish_intake)

data$Fish_score <- ifelse( (data$Oily_fish_intake/1 + data$Non_oily_fish_intake/1) >= 2, 1, 0)
summary(as.factor(data$Fish_score))


#Unprocessed meat:<= 2 servings/week of unprocessed red meats
data$Poultry_intake <-
  plyr::mapvalues(
    data$Poultry_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(0, 0.5, 1, 3, NA, 5.5, NA, NA, 7)
  )

summary(as.factor(data$Poultry_intake))
data <- data[!is.na(data$Poultry_intake), ]
nrow(data) #48423 - 7 = 48416

data$Beef_intake <-
  plyr::mapvalues(
    data$Beef_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(0, 0.5, 1, 3, NA, 5.5, NA, NA, 7)
  )
summary(as.factor(data$Beef_intake))
data <- data[!is.na(data$Beef_intake), ]
nrow(data) #48416-32 = 48384

data$Lamb_or_mutton_intake <-
  plyr::mapvalues(
    data$Lamb_or_mutton_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(0, 0.5, 1, 3, NA, 5.5, NA, NA, 7)
  )
summary(as.factor(data$Lamb_or_mutton_intake))
data <- data[!is.na(data$Lamb_or_mutton_intake), ]
nrow(data) #48384-47 = 48337

data$Pork_intake <-
  plyr::mapvalues(
    data$Pork_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(0, 0.5, 1, 3, NA, 5.5, NA, NA, 7)
  )
summary(as.factor(data$Pork_intake))
data <- data[!is.na(data$Pork_intake), ]
nrow(data) #48337 - 34 = 48303

#Unprocessed meat score
data$Poultry_intake <- as.numeric(data$Poultry_intake)
data$Beef_intake <- as.numeric(data$Beef_intake)
data$Pork_intake <- as.numeric(data$Pork_intake)
data$Lamb_or_mutton_intake <- as.numeric(data$Lamb_or_mutton_intake)

data$Unprocessed_meats_score <- ifelse( (data$Poultry_intake/1 + data$Beef_intake/1 + data$Pork_intake/1 + data$Lamb_or_mutton_intake/1 ) <= 2, 1, 0)
summary(as.factor(data$Unprocessed_meats_score))


#Processed meats: <= 1 serving/week of processed meats
data$Processed_meat_intake <-
  plyr::mapvalues(
    data$Processed_meat_intake,
    from = c(
      "Never",
      "Less than once a week",
      "Once a week",
      "2-4 times a week",
      "Do not know",
      "5-6 times a week",
      "",
      "Prefer not to answer",
      "Once or more daily"
    ),
    to =  c(0, 0.5, 1, 3, NA, 5.5, NA, NA, 7)
  )
summary(as.factor(data$Processed_meat_intake))
data <- data[!is.na(data$Processed_meat_intake), ]
nrow(data) #48303 - 9 = 48294

#Processed meat score
data$Processed_meat_intake <- as.numeric(data$Processed_meat_intake)

data$Processed_meats_score <- ifelse( (data$Processed_meat_intake/1 ) <= 1, 1, 0)
summary(as.factor(data$Processed_meats_score))


#Whole grains: >= 3 servings/day of whole grains
unique(data$Bread_intake)
data$Bread_intake <-
  plyr::mapvalues(
    data$Bread_intake,
    from = c(
      "Do not know",
      "Less than one",
      "Prefer not to answer",
      ""
    ),
    to =  c(NA, 0.5, NA, NA)
  )
summary(as.factor(data$Bread_intake))
data <- data[!is.na(data$Bread_intake), ]
nrow(data) #48294 - 141 = 48153

unique(data$Cereal_intake)
data$Cereal_intake <-
  plyr::mapvalues(
    data$Cereal_intake,
    from = c(
      "Less than one",
      "Do not know",
      "Prefer not to answer",
      ""
    ),
    to =  c(0.5, NA, NA, NA)
  )
summary(as.factor(data$Cereal_intake))
data <- data[!is.na(data$Cereal_intake), ]
nrow(data) # 48153-16 = 48137

#Whole grains score
data$Bread_intake <- as.numeric(data$Bread_intake)
data$Cereal_intake <- as.numeric(data$Cereal_intake)

data$i <- 0
data$j <- 0

# 对Bread_type列进行匹配和赋值
data$i[data$Bread_type %in% c("Wholemeal or wholegrain")] <- data$Bread_intake[data$Bread_type %in% c("Wholemeal or wholegrain")]

# 对Cereal_type列进行匹配和赋值
data$j[data$Cereal_type %in% c("Muesli", "Bran cereal (e.g. All Bran, Branflakes)", "Oat cereal (e.g. Ready Brek, porridge)")] <- data$Cereal_intake[data$Cereal_type %in% c("Muesli", "Bran cereal (e.g. All Bran, Branflakes)", "Oat cereal (e.g. Ready Brek, porridge)")]

# 计算score列的值
data$i <- as.numeric(data$i)
data$j <- as.numeric(data$j)
data$Whole_grains_score <- ifelse(data$i + data$j >= 21, 1, 0)
data <- subset(data, select = -c(i, j))
summary(as.factor(data$Whole_grains_score))


#Refined grains:  <= 2 servings/day of refined grains
data$i <- 0
data$j <- 0
data$i[!(data$Bread_type %in% c("Wholemeal or wholegrain"))] <- data$Bread_intake[!(data$Bread_type %in% c("Wholemeal or wholegrain"))]
data$j[!(data$Cereal_type %in% c("Muesli", "Bran cereal (e.g. All Bran, Branflakes)", "Oat cereal (e.g. Ready Brek, porridge)"))] <- data$Cereal_intake[!(data$Cereal_type %in% c("Muesli", "Bran cereal (e.g. All Bran, Branflakes)", "Oat cereal (e.g. Ready Brek, porridge)"))]

data$i <- as.numeric(data$i)
data$j <- as.numeric(data$j)
data$Refined_grains_score <- ifelse(data$i + data$j <= 14, 1, 0)
data <- subset(data, select = -c(i, j))
summary(as.factor(data$Refined_grains_score))


#Vegetable oil score
data$Bread_intake <- as.numeric(data$Bread_intake)

data$i <- 0
data$j <- 0

# 对Bread_type列进行匹配和赋值
data$i[data$Spread_type %in% c("Flora Pro-Active/Benecol")] <- data$Bread_intake[data$Spread_type %in% c("Flora Pro-Active/Benecol")]
data$j[data$Non_butter_spread_types %in% c("Soft (tub) margarine","Olive oil based spread (eg: Bertolli)","Polyunsaturated/sunflower oil based spread (eg: Flora)")] <- data$Bread_intake[data$Non_butter_spread_types %in% c("Soft (tub) margarine","Olive oil based spread (eg: Bertolli)","Polyunsaturated/sunflower oil based spread (eg: Flora)")]


# 计算score列的值
data$i <- as.numeric(data$i)
data$j <- as.numeric(data$j)
data$Vegetabe_oil_score <- ifelse(data$i + data$j >= 28, 1, 0)
data <- subset(data, select = -c(i, j))
summary(as.factor(data$Vegetabe_oil_score))



#Dairy score
summary(as.factor(data$Cheese_intake))
summary(as.factor(data$Milk_type_used))
data$i <- ifelse(data$Cheese_intake %in% c("Once or more daily"), 1, 0)
data$j <- ifelse(data$Milk_type_used %in% c("Full cream", "Other type of milk","Semi-skimmed","Skimmed","Soya"), 1, 0)

# 计算score列的值
data$i <- as.numeric(data$i)
data$j <- as.numeric(data$j)
data$Dairy_score <- ifelse(data$i + data$j >= 2, 1, 0)
data <- subset(data, select = -c(i, j))
summary(as.factor(data$Dairy_score))



#Sugar_sweentened_beverages_score
summary(as.factor(data$Never_eat_eggs.dairy_wheat_sugar))
data$Sugar_sweentened_beverages_score <- ifelse(grepl("sugar|all", data$Never_eat_eggs.dairy_wheat_sugar), 0, 1)
summary(as.factor(data$Sugar_sweentened_beverages_score))


#Healthy_diet_score
data$Sugar_sweentened_beverages_score <-  as.numeric(data$Sugar_sweentened_beverages_score)
data$Fish_score <- as.numeric(data$Fish_score )
data$Dairy_score <- as.numeric(data$Dairy_score)
data$Fruit_score <- as.numeric(data$Fruit_score)
data$Vegetable_score <- as.numeric(data$Vegetable_score)
data$Vegetabe_oil_score <- as.numeric(data$Vegetabe_oil_score)
data$Whole_grains_score <- as.numeric(data$Whole_grains_score)
data$Processed_meats_score <- as.numeric(data$Processed_meats_score)
data$Refined_grains_score <- as.numeric(data$Refined_grains_score)
data$Unprocessed_meats_score <- as.numeric(data$Unprocessed_meats_score)

data$Healthy_diet_score <- data$Sugar_sweentened_beverages_score + data$Fish_score + data$Dairy_score + data$Fruit_score + data$Vegetable_score + data$Vegetabe_oil_score + data$Whole_grains_score + data$Processed_meats_score + data$Refined_grains_score + data$Unprocessed_meats_score
data$Healthy_diet_score <- as.numeric(data$Healthy_diet_score)



#BMI
summary(data$BMI_i0)
data <- data[!is.na(data$BMI_i0), ]
data$BMI_i0 <- round(data$BMI_i0,1)
nrow(data) # 48137 - 80 = 48057

#Overall_health_rating
summary(as.factor(data$Overall_health_rating))
data$Overall_health_rating <-
  plyr::mapvalues(
    data$Overall_health_rating,
    from = c(
      "Do not know",
      "Prefer not to answer"
    ),
    to =  c(NA, NA)
  )
summary(as.factor(data$Overall_health_rating))
data <- data[!is.na(data$Overall_health_rating), ]
nrow(data) # 48057-62= 47995


#Pre-existing chronic conditions

#Diabetes
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_E10 <- as.Date(data$Date_E10)
data$Date_E11 <- as.Date(data$Date_E11)
data$Date_E12 <- as.Date(data$Date_E12)
data$Date_E13 <- as.Date(data$Date_E13)
data$Date_E14 <- as.Date(data$Date_E14)

unique(data$Source_E10)

data$pre_existing_diabetes <- ifelse((data$Source_E10 != "" & data$Date_E10 <= data$Endtime_wear) | (data$Source_E11 != "" & data$Date_E11 <= data$Endtime_wear) | (data$Source_E12 != "" & data$Date_E12 <= data$Endtime_wear) | 
                                       (data$Source_E13 != "" & data$Date_E13 <= data$Endtime_wear) | (data$Source_E14 != "" & data$Date_E14 <= data$Endtime_wear), 1, 0)

summary(as.factor(data$pre_existing_diabetes))

#Hypertension
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_I10 <- as.Date(data$Date_I10)
data$Date_I11 <- as.Date(data$Date_I11)
data$Date_I12 <- as.Date(data$Date_I12)
data$Date_I13 <- as.Date(data$Date_I13)
data$Date_I15 <- as.Date(data$Date_I15)

unique(data$Source_I10)

data$pre_existing_hypertension <- ifelse((data$Source_I10 != "" & data$Date_I10 <= data$Endtime_wear) | (data$Source_I11 != "" & data$Date_I11 <= data$Endtime_wear) | (data$Source_I12 != "" & data$Date_I12 <= data$Endtime_wear) | 
                                       (data$Source_I13 != "" & data$Date_I13 <= data$Endtime_wear) | (data$Source_I15 != "" & data$Date_I15 <= data$Endtime_wear), 1, 0)

summary(as.factor(data$pre_existing_hypertension))

#Dyslipidaemia
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_E78 <- as.Date(data$Date_E78)

unique(data$Source_E78)

data$pre_existing_dyslipidaemia <- ifelse((data$Source_E78 != "" & data$Date_E78 <= data$Endtime_wear), 1, 0)

#Ref to self-reported and defined the one NA to 1 
data$pre_existing_dyslipidaemia[is.na(data$pre_existing_dyslipidaemia)] <- 1

summary(as.factor(data$pre_existing_dyslipidaemia))

#Depression
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_F32 <- as.Date(data$Date_F32)

unique(data$Source_F32)

data$pre_existing_depression <- ifelse((data$Source_F32 != "" & data$Date_F32 <= data$Endtime_wear), 1, 0)

summary(as.factor(data$pre_existing_depression))


#Cancers
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_cancer_i0 <- as.Date(data$Date_cancer_i0)
data$Date_cancer_i1 <- as.Date(data$Date_cancer_i1)
data$Date_cancer_i2 <- as.Date(data$Date_cancer_i2)
data$Date_cancer_i3 <- as.Date(data$Date_cancer_i3)
data$Date_cancer_i4 <- as.Date(data$Date_cancer_i4)
data$Date_cancer_i5 <- as.Date(data$Date_cancer_i5)
data$Date_cancer_i6 <- as.Date(data$Date_cancer_i6)
data$Date_cancer_i7 <- as.Date(data$Date_cancer_i7)
data$Date_cancer_i8 <- as.Date(data$Date_cancer_i8)
data$Date_cancer_i9 <- as.Date(data$Date_cancer_i9)
data$Date_cancer_i10 <- as.Date(data$Date_cancer_i10)
data$Date_cancer_i11 <- as.Date(data$Date_cancer_i11)
data$Date_cancer_i12 <- as.Date(data$Date_cancer_i12)
data$Date_cancer_i13 <- as.Date(data$Date_cancer_i13)
data$Date_cancer_i14 <- as.Date(data$Date_cancer_i14)
data$Date_cancer_i15 <- as.Date(data$Date_cancer_i15)
data$Date_cancer_i16 <- as.Date(data$Date_cancer_i16)
data$Date_cancer_i17 <- as.Date(data$Date_cancer_i17)
data$Date_cancer_i18 <- as.Date(data$Date_cancer_i18)
data$Date_cancer_i19 <- as.Date(data$Date_cancer_i19)
data$Date_cancer_i20 <- as.Date(data$Date_cancer_i20)
data$Date_cancer_i21 <- as.Date(data$Date_cancer_i21)

unique(data$Type_cancer_i1)

data$pre_existing_cancers <- ifelse((data$Type_cancer_i0 != "" & data$Date_cancer_i0 <= data$Endtime_wear) | (data$Type_cancer_i1 != "" & data$Date_cancer_i1 <= data$Endtime_wear) | (data$Type_cancer_i2 != "" & data$Date_cancer_i2 <= data$Endtime_wear) | 
                                    (data$Type_cancer_i3 != "" & data$Date_cancer_i3 <= data$Endtime_wear) | (data$Type_cancer_i4 != "" & data$Date_cancer_i4 <= data$Endtime_wear) | (data$Type_cancer_i5 != "" & data$Date_cancer_i5 <= data$Endtime_wear) | 
                                    (data$Type_cancer_i6 != "" & data$Date_cancer_i6 <= data$Endtime_wear) | (data$Type_cancer_i7 != "" & data$Date_cancer_i7 <= data$Endtime_wear) | (data$Type_cancer_i8 != "" & data$Date_cancer_i8 <= data$Endtime_wear) |
                                    (data$Type_cancer_i9 != "" & data$Date_cancer_i9 <= data$Endtime_wear) | (data$Type_cancer_i10 != "" & data$Date_cancer_i10 <= data$Endtime_wear) | (data$Type_cancer_i11 != "" & data$Date_cancer_i11 <= data$Endtime_wear) |
                                    (data$Type_cancer_i12 != "" & data$Date_cancer_i12 <= data$Endtime_wear) | (data$Type_cancer_i13 != "" & data$Date_cancer_i13 <= data$Endtime_wear) | (data$Type_cancer_i14 != "" & data$Date_cancer_i14 <= data$Endtime_wear) |
                                    (data$Type_cancer_i15 != "" & data$Date_cancer_i15 <= data$Endtime_wear) | (data$Type_cancer_i16 != "" & data$Date_cancer_i16 <= data$Endtime_wear) | (data$Type_cancer_i17 != "" & data$Date_cancer_i17 <= data$Endtime_wear) |
                                    (data$Type_cancer_i18 != "" & data$Date_cancer_i18 <= data$Endtime_wear) | (data$Type_cancer_i19 != "" & data$Date_cancer_i19 <= data$Endtime_wear) | (data$Type_cancer_i20 != "" & data$Date_cancer_i20 <= data$Endtime_wear) |
                                    (data$Type_cancer_i21 != "" & data$Date_cancer_i21 <= data$Endtime_wear), 1, 0)

data$pre_existing_cancers[is.na(data$pre_existing_cancers)] <- 0
summary(as.factor(data$pre_existing_cancers))

#Chronic conditions
data$chronic_conditions <- ifelse( (data$pre_existing_cancers+data$pre_existing_depression+data$pre_existing_diabetes+data$pre_existing_dyslipidaemia+data$pre_existing_hypertension) >=1, 1, 0)
summary(as.factor(data$chronic_conditions))


#VTE
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_I26 <- as.Date(data$Date_I26)
data$Date_I80 <- as.Date(data$Date_I80)
data$Date_I81 <- as.Date(data$Date_I81)
data$Date_I82 <- as.Date(data$Date_I82)

unique(data$Source_I26)

data$pre_existing_VTE <- ifelse((data$Source_I26 != "" & data$Date_I26 <= data$Endtime_wear) | (data$Source_I80 != "" & data$Date_I80 <= data$Endtime_wear) | (data$Source_I81 != "" & data$Date_I81 <= data$Endtime_wear) | 
                                           (data$Source_I82 != "" & data$Date_I82 <= data$Endtime_wear), 1, 0)

summary(as.factor(data$pre_existing_VTE))

#incident outcomes- incident VTE
#data$incident_VTE <- ifelse((data$Source_I26 != "" & data$Date_I26 > data$Endtime_wear) | (data$Source_I80 != "" & data$Date_I80 > data$Endtime_wear) | (data$Source_I81 != "" & data$Date_I81 > data$Endtime_wear) | 
#                              (data$Source_I82 != "" & data$Date_I82 > data$Endtime_wear), 1, 0)

#summary(as.factor(data$incident_VTE))

##try
data$incident_VTE <- ifelse((data$Source_I26 == "Hospital admissions data only" & data$Date_I26 > data$Endtime_wear) |
                              (data$Source_I80 == "Hospital admissions data only" & data$Date_I80 > data$Endtime_wear) |
                              (data$Source_I81 == "Hospital admissions data only" & data$Date_I81 > data$Endtime_wear) |
                              (data$Source_I82 == "Hospital admissions data only" & data$Date_I82 > data$Endtime_wear), 
                            1,
                            ifelse((data$Source_I26 == "Hospital admissions data and other source(s)" & data$Date_I26 > data$Endtime_wear) |
                                     (data$Source_I80 == "Hospital admissions data and other source(s)" & data$Date_I80 > data$Endtime_wear) |
                                     (data$Source_I81 == "Hospital admissions data and other source(s)" & data$Date_I81 > data$Endtime_wear) |
                                     (data$Source_I82 == "Hospital admissions data and other source(s)" & data$Date_I82 > data$Endtime_wear), 
                                   2,
                                   ifelse((data$Source_I26 == "Death register only" & data$Date_I26 > data$Endtime_wear) |
                                            (data$Source_I80 == "Death register only" & data$Date_I80 > data$Endtime_wear) |
                                            (data$Source_I81 == "Death register only" & data$Date_I81 > data$Endtime_wear) |
                                            (data$Source_I82 == "Death register only" & data$Date_I82 > data$Endtime_wear), 
                                          3, 0)))

summary(as.factor(data$incident_VTE))
data$incident_VTE <- as.numeric(data$incident_VTE)
data$incident_VTE <- ifelse(data$incident_VTE >= 1, 1, 0)

#VTE_time
data$VTE_time <- apply(data[, c("Date_I26", "Date_I80", "Date_I81", "Date_I82")], 1, function(x) {
  non_na_values <- x[!is.na(x)]
  if (length(non_na_values) > 0) {
    min(non_na_values)
  } else {
    NA
  }
})

data$VTE_time <- as.Date(data$VTE_time)
sum(!is.na(data$VTE_time))



#accel_to_VTE
data$country <- "England"
data$country[(data$UKB_centers == "Edinburgh") |
              (data$UKB_centers == "Glasgow")] <- "Scotland"
data$country[(data$UKB_centers == "Cardiff") |
              (data$UKB_centers == "Swansea") |
              (data$UKB_centers == "Wrexham")] <- "Wales"

summary(as.factor(data$country))


data$censored <- "2022-10-31"
data$censored[data$country == "Wales"] <- "2022-05-31"
data$censored[data$country == "Scotland"] <- "2022-08-31"

#censorship
data$VTE_censor_date <-
  as.Date(data$censored, "%Y-%m-%d") # censor at this date as after this date while there may be data it is incomplete

#incident event
data$VTE_censor_date[data$incident_VTE == 1] <-
  pmin(data$VTE_censor_date[data$incident_VTE == 1], as.Date(data$VTE_time[data$incident_VTE == 1], "%Y-%m-%d"))

#death
data$Date_death_i0 <- as.Date(data$Date_death_i0)
data$has_died <- ifelse(!is.na(data$Date_death_i0), 1, 0)
summary(as.factor(data$has_died))

data$VTE_censor_date[data$has_died == 1] <-
  pmin(data$VTE_censor_date[data$has_died == 1], as.Date(data$Date_death_i0[data$has_died == 1],  "%Y-%m-%d"))

#lost-to-follow-up
data$Date_lost_to_followup <- as.Date(data$Date_lost_to_followup)
data$has_lost_follow <- ifelse(!is.na(data$Date_lost_to_followup), 1, 0)
summary(as.factor(data$has_lost_follow))

data$VTE_censor_date[data$has_lost_follow == 1] <-
  pmin(data$VTE_censor_date[data$has_lost_follow == 1], as.Date(data$Date_lost_to_followup[data$has_lost_follow == 1],  "%Y-%m-%d"))

rows <- which(data$VTE_censor_date < data$Endtime_wear)

data <- data[-rows, ]
nrow(data)# 47995 - 51 = 47944


#
#exclude the prevalent VTE
data <- data[data$pre_existing_VTE == 0, ]
nrow(data) # 47944 - 1583 = 46361



#time_accel_to_VTE
data$VTE_censor_date <- as.Date(data$VTE_censor_date)
data$Endtime_wear <- as.Date(data$Endtime_wear)

data$time_accel_to_VTE = as.numeric(data$VTE_censor_date - data$Endtime_wear)/365.25

summary(data$time_accel_to_VTE)


names(data)
nrow(data)



#enhanced PRS
data$enhanced_PRS <- as.numeric(data$enhanced_PRS)
summary(data$enhanced_PRS)

#genetic principle components
data$genetic_a1 <- as.numeric(data$genetic_a1)
data$genetic_a2 <- as.numeric(data$genetic_a2)
data$genetic_a3 <- as.numeric(data$genetic_a3)
data$genetic_a4 <- as.numeric(data$genetic_a4)
data$genetic_a5 <- as.numeric(data$genetic_a5)
data$genetic_a6 <- as.numeric(data$genetic_a6)
data$genetic_a7 <- as.numeric(data$genetic_a7)
data$genetic_a8 <- as.numeric(data$genetic_a8)
data$genetic_a9 <- as.numeric(data$genetic_a9)
data$genetic_a10 <- as.numeric(data$genetic_a10)

sum(is.na(data$genetic_a1))
sum(is.na(data$genetic_a2))
sum(is.na(data$genetic_a3))
sum(is.na(data$genetic_a4))
sum(is.na(data$genetic_a5))
sum(is.na(data$genetic_a6))
sum(is.na(data$genetic_a7))
sum(is.na(data$genetic_a8))
sum(is.na(data$genetic_a9))
sum(is.na(data$genetic_a10))

#genotypying array
data$genotyping_batch <- as.character(data$genotyping_batch)
data$genotyping_batch[grepl("(?i)Batch", data$genotyping_batch)] <- "Axiom"
data$genotyping_batch[grepl("(?i)UKBiLEVEAX", data$genotyping_batch)] <- "BiLEVE"
data$genotyping_batch <- as.factor(data$genotyping_batch)
summary(data$genotyping_batch)



#mean(data$time_accel_to_VTE)
#sd(data$time_accel_to_VTE)

#mean(data$Age)
#sd(data$Age)
#summary(data$Sex)

#

#output
columns_icd <- c("eid",'mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7','MVPA_dayhour_average','MVPA_weekdayhour_average','MVPA_weekendhour_average',
                "Age","Sex","Ethnicity","Education","Townsend_deprivation_index","Average_totalhousehold_income_beforetax",
                "Employment_status","BMI_i0","Smoking_status","Alcohol_frequency","Overall_health_rating","chronic_conditions","Healthy_diet_score",
                "Standard_PRS","Sedentary_overall_average","Light_overall_average","MVPA_overall_average","time_accel_to_VTE","incident_VTE",
                "enhanced_PRS",'genotyping_batch',
                'genetic_a1','genetic_a2','genetic_a3','genetic_a4','genetic_a5','genetic_a6','genetic_a7','genetic_a8','genetic_a9','genetic_a10')

data_icd <- data[, columns_icd]

summary(data_icd)


write.csv(data_icd, "E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/data_mvpa.csv", row.names=FALSE)






###################################################################################################################
###################################################################################################################
###################################################################################################################
#Sensitivity analysis:2-years landmark analysis
library(lubridate)

#VTE
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_I26 <- as.Date(data$Date_I26)
data$Date_I80 <- as.Date(data$Date_I80)
data$Date_I81 <- as.Date(data$Date_I81)
data$Date_I82 <- as.Date(data$Date_I82)

unique(data$Source_I26)
data$Endtime_plus2years <- data$Endtime_wear + years(2)
head(data$Endtime_wear)
head(data$Endtime_plus2years)

data$pre_existing_VTE_plus2years <- ifelse((data$Source_I26 != "" & data$Date_I26 <= data$Endtime_plus2years) | 
                                             (data$Source_I80 != "" & data$Date_I80 <= data$Endtime_plus2years) | 
                                             (data$Source_I81 != "" & data$Date_I81 <= data$Endtime_plus2years) | 
                                             (data$Source_I82 != "" & data$Date_I82 <= data$Endtime_plus2years), 1, 0)

#sensitivity outcomes- incident VTE
#data$incident_VTE_landmark <- ifelse((data$Source_I26 != "" & data$Date_I26 > data$Endtime_plus2years) | (data$Source_I80 != "" & data$Date_I80 > data$Endtime_plus2years) | (data$Source_I81 != "" & data$Date_I81 > data$Endtime_plus2years) | 
#                                      (data$Source_I82 != "" & data$Date_I82 > data$Endtime_plus2years), 1, 0)

#summary(as.factor(data$incident_VTE_landmark))

data$incident_VTE_landmark <- ifelse((data$Source_I26 == "Hospital admissions data only" & data$Date_I26 > data$Endtime_plus2years) |
                              (data$Source_I80 == "Hospital admissions data only" & data$Date_I80 > data$Endtime_plus2years) |
                              (data$Source_I81 == "Hospital admissions data only" & data$Date_I81 > data$Endtime_plus2years) |
                              (data$Source_I82 == "Hospital admissions data only" & data$Date_I82 > data$Endtime_plus2years), 
                            1,
                            ifelse((data$Source_I26 == "Hospital admissions data and other source(s)" & data$Date_I26 > data$Endtime_plus2years) |
                                     (data$Source_I80 == "Hospital admissions data and other source(s)" & data$Date_I80 > data$Endtime_plus2years) |
                                     (data$Source_I81 == "Hospital admissions data and other source(s)" & data$Date_I81 > data$Endtime_plus2years) |
                                     (data$Source_I82 == "Hospital admissions data and other source(s)" & data$Date_I82 > data$Endtime_plus2years), 
                                   2,
                                   ifelse((data$Source_I26 == "Death register only" & data$Date_I26 > data$Endtime_plus2years) |
                                            (data$Source_I80 == "Death register only" & data$Date_I80 > data$Endtime_plus2years) |
                                            (data$Source_I81 == "Death register only" & data$Date_I81 > data$Endtime_plus2years) |
                                            (data$Source_I82 == "Death register only" & data$Date_I82 > data$Endtime_plus2years), 
                                          3, 0)))

summary(as.factor(data$incident_VTE_landmark))
data$incident_VTE_landmark <- as.numeric(data$incident_VTE_landmark)
data$incident_VTE_landmark <- ifelse(data$incident_VTE_landmark >= 1, 1, 0)




#VTE_time
data$VTE_time <- apply(data[, c("Date_I26", "Date_I80", "Date_I81", "Date_I82")], 1, function(x) {
  non_na_values <- x[!is.na(x)]
  if (length(non_na_values) > 0) {
    min(non_na_values)
  } else {
    NA
  }
})

data$VTE_time <- as.Date(data$VTE_time)
sum(!is.na(data$VTE_time))



#accel_to_VTE
data$country <- "England"
data$country[(data$UKB_centers == "Edinburgh") |
               (data$UKB_centers == "Glasgow")] <- "Scotland"
data$country[(data$UKB_centers == "Cardiff") |
               (data$UKB_centers == "Swansea") |
               (data$UKB_centers == "Wrexham")] <- "Wales"

summary(as.factor(data$country))


data$censored <- "2022-10-31"
data$censored[data$country == "Wales"] <- "2022-05-31"
data$censored[data$country == "Scotland"] <- "2022-08-31"

#censorship
data$VTE_censor_date <-
  as.Date(data$censored, "%Y-%m-%d") # censor at this date as after this date while there may be data it is incomplete

#incident event
data$VTE_censor_date[data$incident_VTE_landmark == 1] <-
  pmin(data$VTE_censor_date[data$incident_VTE_landmark == 1], as.Date(data$VTE_time[data$incident_VTE_landmark == 1], "%Y-%m-%d"))

#death
data$Date_death_i0 <- as.Date(data$Date_death_i0)
data$has_died <- ifelse(!is.na(data$Date_death_i0), 1, 0)
summary(as.factor(data$has_died))

data$VTE_censor_date[data$has_died == 1] <-
  pmin(data$VTE_censor_date[data$has_died == 1], as.Date(data$Date_death_i0[data$has_died == 1],  "%Y-%m-%d"))

#lost-to-follow-up
data$Date_lost_to_followup <- as.Date(data$Date_lost_to_followup)
data$has_lost_follow <- ifelse(!is.na(data$Date_lost_to_followup), 1, 0)
summary(as.factor(data$has_lost_follow))

data$VTE_censor_date[data$has_lost_follow == 1] <-
  pmin(data$VTE_censor_date[data$has_lost_follow == 1], as.Date(data$Date_lost_to_followup[data$has_lost_follow == 1],  "%Y-%m-%d"))
rows <- which(data$VTE_censor_date < data$Endtime_plus2years)
nrow(data)
data <- data[-rows, ]
nrow(data)# 47995 - 356 = 47639


#
#exclude the prevalent VTE
data <- data[data$pre_existing_VTE_plus2years == 0, ]
nrow(data) # 47639 - 1748 = 45891



#time_accel_to_VTE
data$VTE_censor_date <- as.Date(data$VTE_censor_date)
data$Endtime_plus2years <- as.Date(data$Endtime_plus2years)
data$time_accel_to_VTE = as.numeric(data$VTE_censor_date - data$Endtime_plus2years)/365.25
summary(data$time_accel_to_VTE)

names(data)
nrow(data) #45891
summary(as.factor(data$incident_VTE_landmark))

#genetic principle components
data$genetic_a1 <- as.numeric(data$genetic_a1)
data$genetic_a2 <- as.numeric(data$genetic_a2)
data$genetic_a3 <- as.numeric(data$genetic_a3)
data$genetic_a4 <- as.numeric(data$genetic_a4)
data$genetic_a5 <- as.numeric(data$genetic_a5)
data$genetic_a6 <- as.numeric(data$genetic_a6)
data$genetic_a7 <- as.numeric(data$genetic_a7)
data$genetic_a8 <- as.numeric(data$genetic_a8)
data$genetic_a9 <- as.numeric(data$genetic_a9)
data$genetic_a10 <- as.numeric(data$genetic_a10)

sum(is.na(data$genetic_a1))
sum(is.na(data$genetic_a2))
sum(is.na(data$genetic_a3))
sum(is.na(data$genetic_a4))
sum(is.na(data$genetic_a5))
sum(is.na(data$genetic_a6))
sum(is.na(data$genetic_a7))
sum(is.na(data$genetic_a8))
sum(is.na(data$genetic_a9))
sum(is.na(data$genetic_a10))

#genotypying array
data$genotyping_batch <- as.character(data$genotyping_batch)
data$genotyping_batch[grepl("(?i)Batch", data$genotyping_batch)] <- "Axiom"
data$genotyping_batch[grepl("(?i)UKBiLEVEAX", data$genotyping_batch)] <- "BiLEVE"
data$genotyping_batch <- as.factor(data$genotyping_batch)
summary(data$genotyping_batch)



#data_landmark
columns_landmark <- c("eid",'mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7','MVPA_dayhour_average','MVPA_weekdayhour_average','MVPA_weekendhour_average',
                      "Age","Sex","Ethnicity","Education","Townsend_deprivation_index","Average_totalhousehold_income_beforetax",
                      "Employment_status","BMI_i0","Smoking_status","Alcohol_frequency","Overall_health_rating","chronic_conditions","Healthy_diet_score",
                      "Standard_PRS","Sedentary_overall_average","Light_overall_average","MVPA_overall_average","time_accel_to_VTE","incident_VTE_landmark",
                      'genotyping_batch','genetic_a1','genetic_a2','genetic_a3','genetic_a4','genetic_a5','genetic_a6','genetic_a7','genetic_a8','genetic_a9','genetic_a10')

data_landmark <- data[, columns_landmark]

summary(data_landmark)


write.csv(data_landmark, "E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_landmark.csv", row.names=FALSE)


###################################################################################################################
###################################################################################################################
###################################################################################################################







##################################################################################################################
##################################################################################################################
#Sensitivity analysis-shift worker

data$shift_worked <-
  plyr::mapvalues(
    data$shift_worked,
    from = c(
      "",
      "Do not know",
      "Prefer not to answer",
      "Always",
      'Never/rarely',
      'Sometimes',
      'Usually'
    ),
    to =  c(0,NA,NA,1,0,1,1)
  )
summary(as.factor(data$shift_worked))

data <- data[!is.na(data$shift_worked), ]
nrow(data) # 47966


#VTE
data$Endtime_wear <- as.Date(data$Endtime_wear)
data$Date_I26 <- as.Date(data$Date_I26)
data$Date_I80 <- as.Date(data$Date_I80)
data$Date_I81 <- as.Date(data$Date_I81)
data$Date_I82 <- as.Date(data$Date_I82)

unique(data$Source_I26)

data$pre_existing_VTE <- ifelse((data$Source_I26 != "" & data$Date_I26 <= data$Endtime_wear) | (data$Source_I80 != "" & data$Date_I80 <= data$Endtime_wear) | (data$Source_I81 != "" & data$Date_I81 <= data$Endtime_wear) | 
                                  (data$Source_I82 != "" & data$Date_I82 <= data$Endtime_wear), 1, 0)

summary(as.factor(data$pre_existing_VTE))

#incident outcomes- incident VTE
data$incident_VTE <- ifelse((data$Source_I26 == "Hospital admissions data only" & data$Date_I26 > data$Endtime_wear) |
                              (data$Source_I80 == "Hospital admissions data only" & data$Date_I80 > data$Endtime_wear) |
                              (data$Source_I81 == "Hospital admissions data only" & data$Date_I81 > data$Endtime_wear) |
                              (data$Source_I82 == "Hospital admissions data only" & data$Date_I82 > data$Endtime_wear), 
                            1,
                            ifelse((data$Source_I26 == "Hospital admissions data and other source(s)" & data$Date_I26 > data$Endtime_wear) |
                                     (data$Source_I80 == "Hospital admissions data and other source(s)" & data$Date_I80 > data$Endtime_wear) |
                                     (data$Source_I81 == "Hospital admissions data and other source(s)" & data$Date_I81 > data$Endtime_wear) |
                                     (data$Source_I82 == "Hospital admissions data and other source(s)" & data$Date_I82 > data$Endtime_wear), 
                                   2,
                                   ifelse((data$Source_I26 == "Death register only" & data$Date_I26 > data$Endtime_wear) |
                                            (data$Source_I80 == "Death register only" & data$Date_I80 > data$Endtime_wear) |
                                            (data$Source_I81 == "Death register only" & data$Date_I81 > data$Endtime_wear) |
                                            (data$Source_I82 == "Death register only" & data$Date_I82 > data$Endtime_wear), 
                                          3, 0)))

summary(as.factor(data$incident_VTE))
data$incident_VTE <- as.numeric(data$incident_VTE)
data$incident_VTE <- ifelse(data$incident_VTE >= 1, 1, 0)


#VTE_time
data$VTE_time <- apply(data[, c("Date_I26", "Date_I80", "Date_I81", "Date_I82")], 1, function(x) {
  non_na_values <- x[!is.na(x)]
  if (length(non_na_values) > 0) {
    min(non_na_values)
  } else {
    NA
  }
})

data$VTE_time <- as.Date(data$VTE_time)
sum(!is.na(data$VTE_time))



#accel_to_VTE
data$country <- "England"
data$country[(data$UKB_centers == "Edinburgh") |
               (data$UKB_centers == "Glasgow")] <- "Scotland"
data$country[(data$UKB_centers == "Cardiff") |
               (data$UKB_centers == "Swansea") |
               (data$UKB_centers == "Wrexham")] <- "Wales"

summary(as.factor(data$country))


data$censored <- "2022-10-31"
data$censored[data$country == "Wales"] <- "2022-05-31"
data$censored[data$country == "Scotland"] <- "2022-08-31"

#censorship
data$VTE_censor_date <-
  as.Date(data$censored, "%Y-%m-%d") # censor at this date as after this date while there may be data it is incomplete

#incident event
data$VTE_censor_date[data$incident_VTE == 1] <-
  pmin(data$VTE_censor_date[data$incident_VTE == 1], as.Date(data$VTE_time[data$incident_VTE == 1], "%Y-%m-%d"))

#death
data$Date_death_i0 <- as.Date(data$Date_death_i0)
data$has_died <- ifelse(!is.na(data$Date_death_i0), 1, 0)
summary(as.factor(data$has_died))

data$VTE_censor_date[data$has_died == 1] <-
  pmin(data$VTE_censor_date[data$has_died == 1], as.Date(data$Date_death_i0[data$has_died == 1],  "%Y-%m-%d"))

#lost-to-follow-up
data$Date_lost_to_followup <- as.Date(data$Date_lost_to_followup)
data$has_lost_follow <- ifelse(!is.na(data$Date_lost_to_followup), 1, 0)
summary(as.factor(data$has_lost_follow))

data$VTE_censor_date[data$has_lost_follow == 1] <-
  pmin(data$VTE_censor_date[data$has_lost_follow == 1], as.Date(data$Date_lost_to_followup[data$has_lost_follow == 1],  "%Y-%m-%d"))

rows <- which(data$VTE_censor_date < data$Endtime_wear)

data <- data[-rows, ]
nrow(data)# 47966 - 53 = 47915


#excluede the shifted worker
data <- data[data$shift_worked == 0, ]
nrow(data) # 47915 - 3030 = 44885


#
#exclude the prevalent VTE
data <- data[data$pre_existing_VTE == 0, ]
nrow(data) # 44885 - 1478 = 43407




#time_accel_to_VTE
data$VTE_censor_date <- as.Date(data$VTE_censor_date)
data$Endtime_wear <- as.Date(data$Endtime_wear)

data$time_accel_to_VTE = as.numeric(data$VTE_censor_date - data$Endtime_wear)/365.25

summary(data$time_accel_to_VTE)


names(data)
nrow(data)

#genetic principle components
data$genetic_a1 <- as.numeric(data$genetic_a1)
data$genetic_a2 <- as.numeric(data$genetic_a2)
data$genetic_a3 <- as.numeric(data$genetic_a3)
data$genetic_a4 <- as.numeric(data$genetic_a4)
data$genetic_a5 <- as.numeric(data$genetic_a5)
data$genetic_a6 <- as.numeric(data$genetic_a6)
data$genetic_a7 <- as.numeric(data$genetic_a7)
data$genetic_a8 <- as.numeric(data$genetic_a8)
data$genetic_a9 <- as.numeric(data$genetic_a9)
data$genetic_a10 <- as.numeric(data$genetic_a10)

sum(is.na(data$genetic_a1))
sum(is.na(data$genetic_a2))
sum(is.na(data$genetic_a3))
sum(is.na(data$genetic_a4))
sum(is.na(data$genetic_a5))
sum(is.na(data$genetic_a6))
sum(is.na(data$genetic_a7))
sum(is.na(data$genetic_a8))
sum(is.na(data$genetic_a9))
sum(is.na(data$genetic_a10))

#genotypying array
data$genotyping_batch <- as.character(data$genotyping_batch)
data$genotyping_batch[grepl("(?i)Batch", data$genotyping_batch)] <- "Axiom"
data$genotyping_batch[grepl("(?i)UKBiLEVEAX", data$genotyping_batch)] <- "BiLEVE"
data$genotyping_batch <- as.factor(data$genotyping_batch)
summary(data$genotyping_batch)




#data_shift
columns_shift <-  c("eid",'mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7','MVPA_dayhour_average','MVPA_weekdayhour_average','MVPA_weekendhour_average',
                   "Age","Sex","Ethnicity","Education","Townsend_deprivation_index","Average_totalhousehold_income_beforetax",
                   "Employment_status","BMI_i0","Smoking_status","Alcohol_frequency","Overall_health_rating","chronic_conditions","Healthy_diet_score",
                   "Standard_PRS","Sedentary_overall_average","Light_overall_average","MVPA_overall_average","time_accel_to_VTE","incident_VTE",
                   'genotyping_batch','genetic_a1','genetic_a2','genetic_a3','genetic_a4','genetic_a5','genetic_a6','genetic_a7','genetic_a8','genetic_a9','genetic_a10')



data_shift <- data[, columns_shift]

summary(data_shift)


write.csv(data_shift, "E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_shift.csv", row.names=FALSE)
##################################################################################################################
##################################################################################################################






##################################################################################################################
##################################################################################################################
#Sensitivity analysis-sleep
#
summary(as.numeric(data$selfreported_sleep))
data$selfreported_sleep <- as.numeric(data$selfreported_sleep)

data <- data[!is.na(data$selfreported_sleep), ]
nrow(data) # 46361-52=46309



#data_sleep
columns_sleep <- c("eid",'mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7','MVPA_dayhour_average','MVPA_weekdayhour_average','MVPA_weekendhour_average',
                      'selfreported_sleep',
                      "Age","Sex","Ethnicity","Education","Townsend_deprivation_index","Average_totalhousehold_income_beforetax",
                      "Employment_status","BMI_i0","Smoking_status","Alcohol_frequency","Overall_health_rating","chronic_conditions","Healthy_diet_score",
                      "Standard_PRS","Sedentary_overall_average","Light_overall_average","MVPA_overall_average","time_accel_to_VTE","incident_VTE",
                      'genotyping_batch','genetic_a1','genetic_a2','genetic_a3','genetic_a4','genetic_a5','genetic_a6','genetic_a7','genetic_a8','genetic_a9','genetic_a10')

data_sleep <- data[, columns_sleep]

summary(data_sleep)


write.csv(data_sleep, "E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_sleep.csv", row.names=FALSE)

##################################################################################################################
##################################################################################################################










##################################################################################################################
##################################################################################################################

#MICE：保留covariate的NA，排除prevalent VTE，得到 n = 55011，输出csv文件，进行插补后敏感性分析

#output
data$Education <- as.factor(data$Education)
data$Average_totalhousehold_income_beforetax <- as.factor(data$Average_totalhousehold_income_beforetax)
data$Employment_status <- as.factor(data$Employment_status)
data$Smoking_status <- as.factor(data$Smoking_status)
data$Alcohol_frequency <- as.factor(data$Alcohol_frequency)
data$Overall_health_rating <- as.factor(data$Overall_health_rating)
data$chronic_conditions <- as.factor(data$chronic_conditions)
data$incident_VTE <- as.factor(data$incident_VTE)

columns_mice <- c('mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7','MVPA_dayhour_average',
                 "Age","Sex","Ethnicity","Education","Townsend_deprivation_index","Average_totalhousehold_income_beforetax",
                 "Employment_status","BMI_i0","Smoking_status","Alcohol_frequency","Overall_health_rating","chronic_conditions","Healthy_diet_score",
                 "Standard_PRS","Sedentary_overall_average","Light_overall_average","MVPA_overall_average","time_accel_to_VTE","incident_VTE",
                 'genotyping_batch',
                 'genetic_a1','genetic_a2','genetic_a3','genetic_a4','genetic_a5','genetic_a6','genetic_a7','genetic_a8','genetic_a9','genetic_a10')

data_mice <- data[, columns_mice]

summary(data_mice)


write.csv(data_mice, "E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_mice.csv", row.names=FALSE)










