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

#data_impute load
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_landmark.csv", header = TRUE)
summary(data_impute)

names(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)



##########################################################################################
#绘制RCS曲线
#RCS_fully_adjusted
library(rms)
dd1 <- datadist(data_impute)
options(datadist="dd1")
ref_value <- median(data_impute$MVPA_Meantime, na.rm = TRUE)
dd1$limits$MVPA_Meantime[3] <- ref_value

fit.cox <- cph(Surv(time_accel_to_VTE, incident_VTE_landmark==1)~ rcs(MVPA_Meantime,3)+ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                 Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                 Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + 
                 SB_Meantime+ genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10, data=data_impute, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, MVPA_Meantime, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(MVPA_Meantime,yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="#4875c7") +
  geom_ribbon(data=HR, aes(MVPA_Meantime,ymin = lower, ymax = upper),alpha = 0.1,fill="#2C6F99") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('MVPA, h/d') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) + #坐标轴标题的文字大小
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),  # 设置主标题字体大小
        axis.text.x = element_text(size = 13),  #坐标轴标签的文字大小
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) + #坐标轴标签的文字大小
  coord_cartesian(xlim = c(0, 3.2)) +
  scale_y_continuous(breaks = c(seq(0.5, 1, 0.1), seq(1, 1.2, 0.1)), 
                     labels = function(x) ifelse(x %in% c(0.6, 0.7, 0.8, 0.9, 1.1), '', x))
P2




################################################################################################
#Table 2


# 计算MVPA_meantime变量的第10、50和90百分位数
percentiles <- quantile(data_impute$MVPA_Meantime, probs = c(0.1, 0.5, 0.9))

percentiles
# 使用cut函数将MVPA_meantime变量分为四个等级
data_impute$MVPA_meantime_quartile <- cut(data_impute$MVPA_Meantime, 
                                          breaks = c(-Inf, percentiles, Inf), 
                                          labels = c("Q1", "Q2", "Q3", "Q4"), 
                                          include.lowest = TRUE)

# 查看新创建的MVPA_meantime_quartile变量
head(data_impute$MVPA_meantime_quartile)

#Person-years计算
summary(data_impute$MVPA_meantime_quartile)

Q1_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q1"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q2"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q3"]
sum_Q3 <- sum(Q3_group)
sum_Q3

Q4_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
Q1_cases_VTE <- sum(data_impute$incident_VTE_landmark[data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE_landmark == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_impute$incident_VTE_landmark[data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE_landmark == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_impute$incident_VTE_landmark[data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE_landmark == 1])
Q3_cases_VTE

Q4_cases_VTE <- sum(data_impute$incident_VTE_landmark[data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE_landmark == 1])
Q4_cases_VTE

# 基于Q1作为参考组建立 Cox 比例风险回归模型:sex+age
mod_vte_min <- coxph(Surv(time_accel_to_VTE, incident_VTE_landmark) ~ Age + Sex + 
                       MVPA_meantime_quartile, data = data_impute)

summary(mod_vte_min)



# 基于Q1作为参考组建立 Cox 比例风险回归模型:全变量
mod_vte <- coxph(Surv(time_accel_to_VTE, incident_VTE_landmark) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_meantime_quartile, data = data_impute)

summary(mod_vte)


###############################################################################################


#data_ww load

data_ww <- data_impute
summary(data_ww)

names(data_ww)


#
data_ww$mvpa1_h <- 24*data_ww$mvpa1
data_ww$mvpa2_h <- 24*data_ww$mvpa2
data_ww$mvpa3_h <- 24*data_ww$mvpa3
data_ww$mvpa4_h <- 24*data_ww$mvpa4
data_ww$mvpa5_h <- 24*data_ww$mvpa5
data_ww$mvpa6_h <- 24*data_ww$mvpa6
data_ww$mvpa7_h <- 24*data_ww$mvpa7

data_ww <- subset(data_ww, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_ww)[names(data_ww) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_ww$mvpa1)
names(data_ww)

#
data_ww$mvpa_daily_total <- rowSums(data_ww[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_ww$mvpa_daily_total)

attach(data_ww)
data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total / 2)
  )), 1, 0)
detach(data_ww)



##150min =2.5h
data_ww$who_acc <- ifelse(data_ww$mvpa_daily_total >= 2.5, 1, 0)

data_ww$who_acc_ww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$who_acc_nww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 0, 1, 0)


data_ww$activity_group <- ifelse(data_ww$who_acc == 0, 'Inactive',
                                 ifelse(data_ww$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_ww$activity_group))

##150min =2.5h
data_ww$activity_group <- as.factor(data_ww$activity_group)
data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE_landmark) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group == "Inactive" & data_ww$incident_VTE_landmark == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group == "Active - WW" & data_ww$incident_VTE_landmark == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group == "Active - Regular" & data_ww$incident_VTE_landmark == 1])
activeregular_cases_VTE


#Median:3.84h
data_ww$above_median <- ifelse(data_ww$mvpa_daily_total >= quantile(data_ww$mvpa_daily_total, 0.5), 1, 0)
data_ww$median_ww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$median_nww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 0, 1, 0)
data_ww$activity_group_median <- ifelse(data_ww$above_median == 0, 'Inactive',
                                        ifelse(data_ww$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_ww$activity_group_median))

#Median:3.84h
data_ww$activity_group_median <- as.factor(data_ww$activity_group_median)
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE_landmark) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group_median, data = data_ww)

summary(mod_ww_median)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE_landmark == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE_landmark == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE_landmark[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE_landmark == 1])
activeregular_cases_VTE














##################################################################################################################################################################
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

#data_impute load
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_shift.csv", header = TRUE)
summary(data_impute)

names(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)



##########################################################################################
#绘制RCS曲线
#RCS_fully_adjusted
library(rms)
dd1 <- datadist(data_impute)
options(datadist="dd1")
ref_value <- median(data_impute$MVPA_Meantime, na.rm = TRUE)
dd1$limits$MVPA_Meantime[3] <- ref_value

fit.cox <- cph(Surv(time_accel_to_VTE, incident_VTE==1)~ rcs(MVPA_Meantime,3)+ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                 Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                 Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + 
                 SB_Meantime +genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +, data=data_impute, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, MVPA_Meantime, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(MVPA_Meantime,yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="#4875c7") +
  geom_ribbon(data=HR, aes(MVPA_Meantime,ymin = lower, ymax = upper),alpha = 0.1,fill="#2C6F99") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('MVPA, h/d') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) + #坐标轴标题的文字大小
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),  # 设置主标题字体大小
        axis.text.x = element_text(size = 13),  #坐标轴标签的文字大小
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) + #坐标轴标签的文字大小
  coord_cartesian(xlim = c(0, 3.2)) +
  scale_y_continuous(breaks = c(seq(0.5, 1, 0.1), seq(1, 1.2, 0.1)), 
                     labels = function(x) ifelse(x %in% c(0.6, 0.7, 0.8, 0.9, 1.1), '', x))
P2




################################################################################################
#Table 2


# 计算MVPA_meantime变量的第10、50和90百分位数
percentiles <- quantile(data_impute$MVPA_Meantime, probs = c(0.1, 0.5, 0.9))

percentiles
# 使用cut函数将MVPA_meantime变量分为四个等级
data_impute$MVPA_meantime_quartile <- cut(data_impute$MVPA_Meantime, 
                                          breaks = c(-Inf, percentiles, Inf), 
                                          labels = c("Q1", "Q2", "Q3", "Q4"), 
                                          include.lowest = TRUE)

# 查看新创建的MVPA_meantime_quartile变量
head(data_impute$MVPA_meantime_quartile)

#Person-years计算
summary(data_impute$MVPA_meantime_quartile)

Q1_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q1"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q2"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q3"]
sum_Q3 <- sum(Q3_group)
sum_Q3

Q4_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
Q1_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1])
Q3_cases_VTE

Q4_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1])
Q4_cases_VTE

# 基于Q1作为参考组建立 Cox 比例风险回归模型:sex+age
mod_vte_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + 
                       MVPA_meantime_quartile, data = data_impute)

summary(mod_vte_min)



# 基于Q1作为参考组建立 Cox 比例风险回归模型:全变量
mod_vte <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_meantime_quartile, data = data_impute)

summary(mod_vte)



###############################################################################################

#data_ww load

data_ww <- data_impute
summary(data_ww)

names(data_ww)


#
data_ww$mvpa1_h <- 24*data_ww$mvpa1
data_ww$mvpa2_h <- 24*data_ww$mvpa2
data_ww$mvpa3_h <- 24*data_ww$mvpa3
data_ww$mvpa4_h <- 24*data_ww$mvpa4
data_ww$mvpa5_h <- 24*data_ww$mvpa5
data_ww$mvpa6_h <- 24*data_ww$mvpa6
data_ww$mvpa7_h <- 24*data_ww$mvpa7

data_ww <- subset(data_ww, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_ww)[names(data_ww) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_ww$mvpa1)
names(data_ww)

#
data_ww$mvpa_daily_total <- rowSums(data_ww[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_ww$mvpa_daily_total)

attach(data_ww)
data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total / 2)
  )), 1, 0)
detach(data_ww)



##150min =2.5h
data_ww$who_acc <- ifelse(data_ww$mvpa_daily_total >= 2.5, 1, 0)

data_ww$who_acc_ww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$who_acc_nww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 0, 1, 0)


data_ww$activity_group <- ifelse(data_ww$who_acc == 0, 'Inactive',
                                 ifelse(data_ww$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_ww$activity_group))

##150min =2.5h
data_ww$activity_group <- as.factor(data_ww$activity_group)
data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE


#Median:3.84h
data_ww$above_median <- ifelse(data_ww$mvpa_daily_total >= quantile(data_ww$mvpa_daily_total, 0.5), 1, 0)
data_ww$median_ww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$median_nww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 0, 1, 0)
data_ww$activity_group_median <- ifelse(data_ww$above_median == 0, 'Inactive',
                                        ifelse(data_ww$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_ww$activity_group_median))

#Median:3.84h
data_ww$activity_group_median <- as.factor(data_ww$activity_group_median)
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         activity_group_median, data = data_ww)

summary(mod_ww_median)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE


















###################################################################################################################################################################
#self-reported sleep
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

#data_impute load
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_sleep.csv", header = TRUE)
summary(data_impute)

names(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)



##########################################################################################
#绘制RCS曲线
#RCS_fully_adjusted
library(rms)
dd1 <- datadist(data_impute)
options(datadist="dd1")
ref_value <- median(data_impute$MVPA_Meantime, na.rm = TRUE)
dd1$limits$MVPA_Meantime[3] <- ref_value

fit.cox <- cph(Surv(time_accel_to_VTE, incident_VTE==1)~ rcs(MVPA_Meantime,3)+ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                 Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                 Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + 
                 SB_Meantime + selfreported_sleep + 
                 genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10, data=data_impute, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, MVPA_Meantime, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(MVPA_Meantime,yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="#4875c7") +
  geom_ribbon(data=HR, aes(MVPA_Meantime,ymin = lower, ymax = upper),alpha = 0.1,fill="#2C6F99") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('MVPA, h/d') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) + #坐标轴标题的文字大小
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),  # 设置主标题字体大小
        axis.text.x = element_text(size = 13),  #坐标轴标签的文字大小
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) + #坐标轴标签的文字大小
  coord_cartesian(xlim = c(0, 3.2)) +
  scale_y_continuous(breaks = c(seq(0.5, 1, 0.1), seq(1, 1.2, 0.1)), 
                     labels = function(x) ifelse(x %in% c(0.6, 0.7, 0.8, 0.9, 1.1), '', x))
P2




################################################################################################
#Table 2


# 计算MVPA_meantime变量的第10、50和90百分位数
percentiles <- quantile(data_impute$MVPA_Meantime, probs = c(0.1, 0.5, 0.9))

percentiles
# 使用cut函数将MVPA_meantime变量分为四个等级
data_impute$MVPA_meantime_quartile <- cut(data_impute$MVPA_Meantime, 
                                          breaks = c(-Inf, percentiles, Inf), 
                                          labels = c("Q1", "Q2", "Q3", "Q4"), 
                                          include.lowest = TRUE)

# 查看新创建的MVPA_meantime_quartile变量
head(data_impute$MVPA_meantime_quartile)

#Person-years计算
summary(data_impute$MVPA_meantime_quartile)

Q1_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q1"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q2"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q3"]
sum_Q3 <- sum(Q3_group)
sum_Q3

Q4_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
Q1_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1])
Q3_cases_VTE

Q4_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1])
Q4_cases_VTE

# 基于Q1作为参考组建立 Cox 比例风险回归模型:sex+age
mod_vte_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + 
                       MVPA_meantime_quartile, data = data_impute)

summary(mod_vte_min)



# 基于Q1作为参考组建立 Cox 比例风险回归模型:全变量
mod_vte <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                   selfreported_sleep +
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_meantime_quartile, data = data_impute)

summary(mod_vte)


###############################################################################################
###############################################################################################
###############################################################################################

#data_ww load

data_ww <- data_impute
summary(data_ww)

names(data_ww)


#
data_ww$mvpa1_h <- 24*data_ww$mvpa1
data_ww$mvpa2_h <- 24*data_ww$mvpa2
data_ww$mvpa3_h <- 24*data_ww$mvpa3
data_ww$mvpa4_h <- 24*data_ww$mvpa4
data_ww$mvpa5_h <- 24*data_ww$mvpa5
data_ww$mvpa6_h <- 24*data_ww$mvpa6
data_ww$mvpa7_h <- 24*data_ww$mvpa7

data_ww <- subset(data_ww, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_ww)[names(data_ww) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_ww$mvpa1)
names(data_ww)

#
data_ww$mvpa_daily_total <- rowSums(data_ww[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_ww$mvpa_daily_total)

attach(data_ww)
data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total / 2)
  )), 1, 0)
detach(data_ww)



##150min =2.5h
data_ww$who_acc <- ifelse(data_ww$mvpa_daily_total >= 2.5, 1, 0)

data_ww$who_acc_ww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$who_acc_nww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 0, 1, 0)


data_ww$activity_group <- ifelse(data_ww$who_acc == 0, 'Inactive',
                                 ifelse(data_ww$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_ww$activity_group))

##150min =2.5h
data_ww$activity_group <- as.factor(data_ww$activity_group)
data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                  selfreported_sleep +
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE


#Median:3.84h
data_ww$above_median <- ifelse(data_ww$mvpa_daily_total >= quantile(data_ww$mvpa_daily_total, 0.5), 1, 0)
data_ww$median_ww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$median_nww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 0, 1, 0)
data_ww$activity_group_median <- ifelse(data_ww$above_median == 0, 'Inactive',
                                        ifelse(data_ww$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_ww$activity_group_median))

#Median:3.84h
data_ww$activity_group_median <- as.factor(data_ww$activity_group_median)
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                         selfreported_sleep +
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         activity_group_median, data = data_ww)

summary(mod_ww_median)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE












######################################################################################################
#Enhanced PRS
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/data_mvpa.csv", header = TRUE)
summary(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)


summary(data_impute$enhanced_PRS)


#
data_EnhancedPRS <- data_impute[!is.na(data_impute$enhanced_PRS), ]
nrow(data_EnhancedPRS)


###############################################################################
#绘制RCS曲线
#RCS_fully_adjusted
library(rms)
dd1 <- datadist(data_EnhancedPRS)
options(datadist="dd1")
ref_value <- median(data_EnhancedPRS$MVPA_Meantime, na.rm = TRUE)
dd1$limits$MVPA_Meantime[3] <- ref_value

fit.cox <- cph(Surv(time_accel_to_VTE, incident_VTE==1)~ rcs(MVPA_Meantime,3)+ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                 Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                 Overall_health_rating + chronic_conditions + Healthy_diet_score + enhanced_PRS + LightB_Meantime + 
                 SB_Meantime +
                 genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 , data=data_EnhancedPRS, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, MVPA_Meantime, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(MVPA_Meantime,yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="#4875c7") +
  geom_ribbon(data=HR, aes(MVPA_Meantime,ymin = lower, ymax = upper),alpha = 0.1,fill="#2C6F99") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('MVPA, h/d') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) + #坐标轴标题的文字大小
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),  # 设置主标题字体大小
        axis.text.x = element_text(size = 13),  #坐标轴标签的文字大小
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) + #坐标轴标签的文字大小
  coord_cartesian(xlim = c(0, 3.2)) +
  scale_y_continuous(breaks = c(seq(0.5, 1, 0.1), seq(1, 1.2, 0.1)), 
                     labels = function(x) ifelse(x %in% c(0.6, 0.7, 0.8, 0.9, 1.1), '', x))
P2


################################################################################
#Table 2
# 计算MVPA_meantime变量的第10、50和90百分位数
percentiles <- quantile(data_EnhancedPRS$MVPA_Meantime, probs = c(0.1, 0.5, 0.9))

percentiles
# 使用cut函数将MVPA_meantime变量分为四个等级
data_EnhancedPRS$MVPA_meantime_quartile <- cut(data_EnhancedPRS$MVPA_Meantime, 
                                          breaks = c(-Inf, percentiles, Inf), 
                                          labels = c("Q1", "Q2", "Q3", "Q4"), 
                                          include.lowest = TRUE)

# 查看新创建的MVPA_meantime_quartile变量
head(data_EnhancedPRS$MVPA_meantime_quartile)

#Person-years计算
summary(data_EnhancedPRS$MVPA_meantime_quartile)

Q1_group <- data_EnhancedPRS$time_accel_to_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q1"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_EnhancedPRS$time_accel_to_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q2"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_EnhancedPRS$time_accel_to_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q3"]
sum_Q3 <- sum(Q3_group)
sum_Q3

Q4_group <- data_EnhancedPRS$time_accel_to_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
Q1_cases_VTE <- sum(data_EnhancedPRS$incident_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q1" & data_EnhancedPRS$incident_VTE == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_EnhancedPRS$incident_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q2" & data_EnhancedPRS$incident_VTE == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_EnhancedPRS$incident_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q3" & data_EnhancedPRS$incident_VTE == 1])
Q3_cases_VTE

Q4_cases_VTE <- sum(data_EnhancedPRS$incident_VTE[data_EnhancedPRS$MVPA_meantime_quartile == "Q4" & data_EnhancedPRS$incident_VTE == 1])
Q4_cases_VTE

# 基于Q1作为参考组建立 Cox 比例风险回归模型:sex+age
mod_vte_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + 
                       MVPA_meantime_quartile, data = data_EnhancedPRS)

summary(mod_vte_min)



# 基于Q1作为参考组建立 Cox 比例风险回归模型:全变量
mod_vte <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score + enhanced_PRS + LightB_Meantime + SB_Meantime + 
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_meantime_quartile, data = data_EnhancedPRS)

summary(mod_vte)

#######################################################################################
#
data_ww <- data_EnhancedPRS
summary(data_ww)

names(data_ww)


#
data_ww$mvpa1_h <- 24*data_ww$mvpa1
data_ww$mvpa2_h <- 24*data_ww$mvpa2
data_ww$mvpa3_h <- 24*data_ww$mvpa3
data_ww$mvpa4_h <- 24*data_ww$mvpa4
data_ww$mvpa5_h <- 24*data_ww$mvpa5
data_ww$mvpa6_h <- 24*data_ww$mvpa6
data_ww$mvpa7_h <- 24*data_ww$mvpa7

data_ww <- subset(data_ww, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_ww)[names(data_ww) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_ww$mvpa1)
names(data_ww)

#
data_ww$mvpa_daily_total <- rowSums(data_ww[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_ww$mvpa_daily_total)

attach(data_ww)
data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total / 2)
  )), 1, 0)
detach(data_ww)



##150min =2.5h
data_ww$who_acc <- ifelse(data_ww$mvpa_daily_total >= 2.5, 1, 0)

data_ww$who_acc_ww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$who_acc_nww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 0, 1, 0)


data_ww$activity_group <- ifelse(data_ww$who_acc == 0, 'Inactive',
                                 ifelse(data_ww$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_ww$activity_group))

##150min =2.5h
data_ww$activity_group <- as.factor(data_ww$activity_group)
data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + enhanced_PRS + LightB_Meantime + SB_Meantime +
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE


#Median:3.84h
data_ww$above_median <- ifelse(data_ww$mvpa_daily_total >= quantile(data_ww$mvpa_daily_total, 0.5), 1, 0)
data_ww$median_ww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$median_nww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 0, 1, 0)
data_ww$activity_group_median <- ifelse(data_ww$above_median == 0, 'Inactive',
                                        ifelse(data_ww$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_ww$activity_group_median))

#Median:3.84h
data_ww$activity_group_median <- as.factor(data_ww$activity_group_median)
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + enhanced_PRS + LightB_Meantime + SB_Meantime + 
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         activity_group_median, data = data_ww)

summary(mod_ww_median)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE










#######################################################################################################################################
#######################################################################################################################################
##MICE:插补后完整数据集
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

#data_impute load
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/sensitivity_data/data_mice.csv", header = TRUE)
summary(data_impute)

names(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)


data_impute$mvpa1_h <- 24*data_impute$mvpa1
data_impute$mvpa2_h <- 24*data_impute$mvpa2
data_impute$mvpa3_h <- 24*data_impute$mvpa3
data_impute$mvpa4_h <- 24*data_impute$mvpa4
data_impute$mvpa5_h <- 24*data_impute$mvpa5
data_impute$mvpa6_h <- 24*data_impute$mvpa6
data_impute$mvpa7_h <- 24*data_impute$mvpa7

data_impute <- subset(data_impute, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_impute)[names(data_impute) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_impute$mvpa1)
names(data_impute)

#
data_impute$mvpa_daily_total <- rowSums(data_impute[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_impute$mvpa_daily_total)

attach(data_impute)
data_impute$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total / 2) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total / 2) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total / 2)
  )), 1, 0)
detach(data_impute)



##150min =2.5h
data_impute$who_acc <- ifelse(data_impute$mvpa_daily_total >= 2.5, 1, 0)

data_impute$who_acc_ww <- ifelse(data_impute$who_acc == 1 & data_impute$ww_pattern == 1, 1, 0)
data_impute$who_acc_nww <- ifelse(data_impute$who_acc == 1 & data_impute$ww_pattern == 0, 1, 0)


data_impute$activity_group <- ifelse(data_impute$who_acc == 0, 'Inactive',
                                 ifelse(data_impute$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_impute$activity_group))

##150min =2.5h
data_impute$activity_group <- as.factor(data_impute$activity_group)
data_impute$activity_group <- relevel(data_impute$activity_group, ref = "Inactive")


#Median:3.84h
data_impute$above_median <- ifelse(data_impute$mvpa_daily_total >= quantile(data_impute$mvpa_daily_total, 0.5), 1, 0)
data_impute$median_ww <- ifelse(data_impute$above_median == 1 & data_impute$ww_pattern == 1, 1, 0)
data_impute$median_nww <- ifelse(data_impute$above_median == 1 & data_impute$ww_pattern == 0, 1, 0)
data_impute$activity_group_median <- ifelse(data_impute$above_median == 0, 'Inactive',
                                        ifelse(data_impute$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_impute$activity_group_median))

#Median:3.84h
data_impute$activity_group_median <- as.factor(data_impute$activity_group_median)
data_impute$activity_group_median <- relevel(data_impute$activity_group_median, ref = "Inactive")


# 计算MVPA_meantime变量的第10、50和90百分位数
percentiles <- quantile(data_impute$MVPA_Meantime, probs = c(0.1, 0.5, 0.9))

percentiles
# 使用cut函数将MVPA_meantime变量分为四个等级
data_impute$MVPA_meantime_quartile <- cut(data_impute$MVPA_Meantime, 
                                          breaks = c(-Inf, percentiles, Inf), 
                                          labels = c("Q1", "Q2", "Q3", "Q4"), 
                                          include.lowest = TRUE)

# 查看新创建的MVPA_meantime_quartile变量
head(data_impute$MVPA_meantime_quartile)

#Person-years计算
summary(data_impute$MVPA_meantime_quartile)

Q1_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q1"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q2"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q3"]
sum_Q3 <- sum(Q3_group)
sum_Q3

Q4_group <- data_impute$time_accel_to_VTE[data_impute$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
Q1_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1])
Q3_cases_VTE

Q4_cases_VTE <- sum(data_impute$incident_VTE[data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1])
Q4_cases_VTE



#准备插补

data_ww <- subset(data_impute, select = c(Age, Sex, Ethnicity, Education, Townsend_deprivation_index, 
                                          Average_totalhousehold_income_beforetax, Employment_status,
                                          BMI_i0, Smoking_status, Alcohol_frequency, Overall_health_rating,
                                          chronic_conditions, Healthy_diet_score, Standard_PRS, LightB_Meantime,
                                          SB_Meantime, genotyping_batch, genetic_a1, genetic_a2, genetic_a3, genetic_a4, 
                                          genetic_a5, genetic_a6, genetic_a7, genetic_a8, genetic_a9, genetic_a10, 
                                          activity_group,activity_group_median,
                                          MVPA_meantime_quartile,MVPA_Meantime,
                                          time_accel_to_VTE, incident_VTE)) 

data_ww$Sex <- as.factor(data_ww$Sex)
data_ww$Ethnicity <- as.factor(data_ww$Ethnicity)
data_ww$Education <- as.factor(data_ww$Education)
data_ww$Average_totalhousehold_income_beforetax <- as.factor(data_ww$Average_totalhousehold_income_beforetax)
data_ww$Employment_status <- as.factor(data_ww$Employment_status)
data_ww$Smoking_status <- as.factor(data_ww$Smoking_status)
data_ww$Alcohol_frequency <- as.factor(data_ww$Alcohol_frequency)
data_ww$Overall_health_rating <- as.factor(data_ww$Overall_health_rating)
data_ww$chronic_conditions <- as.factor(data_ww$chronic_conditions)
data_ww$incident_VTE <- as.numeric(data_ww$incident_VTE)
data_ww$genotyping_batch <- as.factor(data_ww$genotyping_batch)
data_ww$MVPA_meantime_quartile <- as.factor(data_ww$MVPA_meantime_quartile)

summary(data_ww)

data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

missing_percentages <- colMeans(is.na(data_ww)) * 100
missing_percentages


library(mice)

imputed_data <- mice(data_ww, m = 10, seed = 123,
                     maxit = 5, print = TRUE)

#Table 1
#min
fit_min <- with(imputed_data, 
                    coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex+
                            MVPA_meantime_quartile))


est_min <-pool(fit_min)

summary(est_min) 

#
summary_data <- data.frame(
  term = c("MVPA_meantime_quartileQ2", "MVPA_meantime_quartileQ3", "MVPA_meantime_quartileQ4"),
  estimate = c(-0.58497570, -0.81580527, -1.02099381),
  std.error = c(0.086866685, 0.090436473, 0.135416105)
)

summary_data$HR <- exp(summary_data$estimate)
summary_data$lower_CI <- exp(summary_data$estimate - 1.96 * summary_data$std.error)
summary_data$upper_CI <- exp(summary_data$estimate + 1.96 * summary_data$std.error)

print(summary_data[, c("term", "HR", "lower_CI", "upper_CI")])



#fully
fit_fully <- with(imputed_data, 
                  coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                        Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                        Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime +
                        genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                        genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                        MVPA_meantime_quartile))


est_fully <-pool(fit_fully)

summary(est_fully) 

#
summary_data <- data.frame(
  term = c("MVPA_meantime_quartileQ2", "MVPA_meantime_quartileQ3", "MVPA_meantime_quartileQ4"),
  estimate = c(-0.3249662678, -0.4059147682, -0.5746548277),
  std.error = c(0.09057656, 0.09951302, 0.1473081)
)

summary_data$HR <- exp(summary_data$estimate)
summary_data$lower_CI <- exp(summary_data$estimate - 1.96 * summary_data$std.error)
summary_data$upper_CI <- exp(summary_data$estimate + 1.96 * summary_data$std.error)

print(summary_data[, c("term", "HR", "lower_CI", "upper_CI")])



###
###
###
#Table 2
#2.5h
fit_group <- with(imputed_data, 
             coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
             Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
             Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime +
             genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
             genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
             activity_group))


est_group <-pool(fit_group)

summary(est_group) 

#
summary_data <- data.frame(
  term = c("Active - Regular", "Active - WW"),
  estimate = c(-0.2074684074, -0.2104382715),
  std.error = c(0.09068372, 0.07275768)
)

summary_data$HR <- exp(summary_data$estimate)
summary_data$lower_CI <- exp(summary_data$estimate - 1.96 * summary_data$std.error)
summary_data$upper_CI <- exp(summary_data$estimate + 1.96 * summary_data$std.error)

print(summary_data[, c("term", "HR", "lower_CI", "upper_CI")])




#median
fit_group_median <- with(imputed_data, 
                    coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                    Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                    Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime +
                    genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                    genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                    activity_group_median))

est_group_median <- pool(fit_group_median)

summary(est_group_median) 

#
summary_data <- data.frame(
  term = c("Active - Regular", "Active - WW"),
  estimate = c(-0.1831612488, -0.1531615326),
  std.error = c(0.09155874, 0.07666074)
)

summary_data$HR <- exp(summary_data$estimate)
summary_data$lower_CI <- exp(summary_data$estimate - 1.96 * summary_data$std.error)
summary_data$upper_CI <- exp(summary_data$estimate + 1.96 * summary_data$std.error)

print(summary_data[, c("term", "HR", "lower_CI", "upper_CI")])








############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
#MVPA cutoffs

#
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/data_mvpa.csv", header = TRUE)
summary(data_impute)

data_impute$SB_Meantime <- 24*data_impute$Sedentary_overall_average
summary(data_impute$SB_Meantime)

data_impute$LightB_Meantime <- 24*data_impute$Light_overall_average
summary(data_impute$LightB_Meantime)

data_impute$MVPA_Meantime <- 24*data_impute$MVPA_overall_average
summary(data_impute$MVPA_Meantime)

data_impute$Sex <- as.factor(data_impute$Sex)
data_impute$chronic_conditions <- as.factor(data_impute$chronic_conditions)



#
data_ww <- data_impute
summary(data_ww)

names(data_ww)


#
data_ww$mvpa1_h <- 24*data_ww$mvpa1
data_ww$mvpa2_h <- 24*data_ww$mvpa2
data_ww$mvpa3_h <- 24*data_ww$mvpa3
data_ww$mvpa4_h <- 24*data_ww$mvpa4
data_ww$mvpa5_h <- 24*data_ww$mvpa5
data_ww$mvpa6_h <- 24*data_ww$mvpa6
data_ww$mvpa7_h <- 24*data_ww$mvpa7

data_ww <- subset(data_ww, select = -c(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7)) 
names(data_ww)[names(data_ww) %in% c("mvpa1_h", "mvpa2_h", "mvpa3_h", "mvpa4_h", "mvpa5_h", "mvpa6_h", "mvpa7_h")] <- paste0("mvpa", 1:7) 

summary(data_ww$mvpa1)
names(data_ww)

#
data_ww$mvpa_daily_total <- rowSums(data_ww[, c("mvpa1", "mvpa2", "mvpa3", "mvpa4", "mvpa5", "mvpa6", "mvpa7")], na.rm = TRUE)
summary(data_ww$mvpa_daily_total)


#55%
attach(data_ww)
data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total * 0.55) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total * 0.55) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total * 0.55) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total * 0.55) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total * 0.55) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total * 0.55) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total * 0.55) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total * 0.55) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total * 0.55) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total * 0.55) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total * 0.55) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total * 0.55) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total * 0.55) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total * 0.55) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total * 0.55) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total * 0.55) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total * 0.55) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total * 0.55) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total * 0.55) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total * 0.55) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total * 0.55)
  )), 1, 0)

data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total * 0.60) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total * 0.60) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total * 0.60) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total * 0.60) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total * 0.60) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total * 0.60) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total * 0.60) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total * 0.60) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total * 0.60) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total * 0.60) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total * 0.60) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total * 0.60) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total * 0.60) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total * 0.60) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total * 0.60) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total * 0.60) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total * 0.60) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total * 0.60) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total * 0.60) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total * 0.60) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total * 0.60)
  )), 1, 0)

data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total * 0.65) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total * 0.65) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total * 0.65) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total * 0.65) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total * 0.65) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total * 0.65) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total * 0.65) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total * 0.65) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total * 0.65) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total * 0.65) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total * 0.65) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total * 0.65) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total * 0.65) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total * 0.65) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total * 0.65) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total * 0.65) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total * 0.65) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total * 0.65) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total * 0.65) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total * 0.65) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total * 0.65)
  )), 1, 0)

data_ww$ww_pattern <- ifelse(
  ((
    (mvpa1 + mvpa2) >= (mvpa_daily_total * 0.70) |
      (mvpa1 + mvpa3) >= (mvpa_daily_total * 0.70) |
      (mvpa1 + mvpa4) >= (mvpa_daily_total * 0.70) |
      (mvpa1 + mvpa5) >= (mvpa_daily_total * 0.70) |
      (mvpa1 + mvpa6) >= (mvpa_daily_total * 0.70) |
      (mvpa1 + mvpa7) >= (mvpa_daily_total * 0.70) |
      (mvpa2 + mvpa3) >= (mvpa_daily_total * 0.70) |
      (mvpa2 + mvpa4) >= (mvpa_daily_total * 0.70) |
      (mvpa2 + mvpa5) >= (mvpa_daily_total * 0.70) |
      (mvpa2 + mvpa6) >= (mvpa_daily_total * 0.70) |
      (mvpa2 + mvpa7) >= (mvpa_daily_total * 0.70) |
      (mvpa3 + mvpa4) >= (mvpa_daily_total * 0.70) |
      (mvpa3 + mvpa5) >= (mvpa_daily_total * 0.70) |
      (mvpa3 + mvpa6) >= (mvpa_daily_total * 0.70) |
      (mvpa3 + mvpa7) >= (mvpa_daily_total * 0.70) |
      (mvpa4 + mvpa5) >= (mvpa_daily_total * 0.70) |
      (mvpa4 + mvpa6) >= (mvpa_daily_total * 0.70) |
      (mvpa4 + mvpa7) >= (mvpa_daily_total * 0.70) |
      (mvpa5 + mvpa6) >= (mvpa_daily_total * 0.70) |
      (mvpa5 + mvpa7) >= (mvpa_daily_total * 0.70) |
      (mvpa6 + mvpa7) >= (mvpa_daily_total * 0.70)
  )), 1, 0)

detach(data_ww)



##150min =2.5h
data_ww$who_acc <- ifelse(data_ww$mvpa_daily_total >= 2.5, 1, 0)

data_ww$who_acc_ww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$who_acc_nww <- ifelse(data_ww$who_acc == 1 & data_ww$ww_pattern == 0, 1, 0)


data_ww$activity_group <- ifelse(data_ww$who_acc == 0, 'Inactive',
                                 ifelse(data_ww$who_acc_ww == 1, 'Active - WW', 'Active - Regular'))

summary(as.factor(data_ww$activity_group))

##150min =2.5h
data_ww$activity_group <- as.factor(data_ww$activity_group)
data_ww$activity_group <- relevel(data_ww$activity_group, ref = "Inactive")

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime +
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE


#Median:3.84h
data_ww$above_median <- ifelse(data_ww$mvpa_daily_total >= quantile(data_ww$mvpa_daily_total, 0.5), 1, 0)
data_ww$median_ww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 1, 1, 0)
data_ww$median_nww <- ifelse(data_ww$above_median == 1 & data_ww$ww_pattern == 0, 1, 0)
data_ww$activity_group_median <- ifelse(data_ww$above_median == 0, 'Inactive',
                                        ifelse(data_ww$median_ww == 1, 'Active - WW', 'Active - Regular'))
summary(as.factor(data_ww$activity_group_median))

#Median:3.84h
data_ww$activity_group_median <- as.factor(data_ww$activity_group_median)
data_ww$activity_group_median <- relevel(data_ww$activity_group_median, ref = "Inactive")

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         activity_group_median, data = data_ww)

summary(mod_ww_median)

#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE



