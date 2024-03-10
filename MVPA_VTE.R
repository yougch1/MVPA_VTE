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
library("survminer")

#data_mvpa load
data_impute <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/data_mvpa.csv", header = TRUE)
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


#mean(data_impute$time_accel_to_VTE)
#sd(data_impute$time_accel_to_VTE)
#sum(data_impute$time_accel_to_VTE)
#mean(data_impute$Age)
#sd(data_impute$Age)
#summary(data_impute$Sex)
#summary(as.factor(data_impute$Ethnicity))


#########################################################################################
# 绘制直方图
vte_positive <- data_impute$MVPA_Meantime[data_impute$incident_VTE == 1]

P1 <- ggplot(data_impute, aes(x=MVPA_Meantime)) +
  geom_histogram(binwidth=0.2, fill = "#4C6F80" ,color="black") +
  geom_segment(data=data.frame(x=vte_positive), aes(x=x, xend=x, y=-Inf, yend=-150), color="black") +  # 使用 geom_segment
  labs(title="Distribution of MVPA time and VTE cases", x="MVPA, h/d", y="No. of participants") +
  coord_cartesian(xlim = c(0, 3.2)) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text = element_text(size=13), #坐标轴上的文本的字体大小
        axis.title = element_text(size=14), # 坐标轴标题的字体大小
        plot.title = element_text(size = 16),  # 设置标题字体大小
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"))
P1



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
                 SB_Meantime + genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 +genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10, data=data_impute, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, MVPA_Meantime, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(MVPA_Meantime,yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="#4C6F80") +
  geom_ribbon(data=HR, aes(MVPA_Meantime,ymin = lower, ymax = upper),alpha = 0.1,fill="#4C6F99") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('MVPA, h/d') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) +  # 设置坐标轴刻度线长度为0.2
  coord_cartesian(xlim = c(0, 3.2)) +
  scale_y_continuous(breaks = c(seq(0.5, 1, 0.1), seq(1, 1.2, 0.1)), 
                     labels = function(x) ifelse(x %in% c(0.6, 0.7, 0.8, 0.9, 1.1), '', x))

P2


################################################################################################
#基线表
library(gtsummary)

# 创建基线表格
table1 <- data_impute %>%
  select(Age, Sex, Ethnicity, Education, Townsend_deprivation_index,
         Average_totalhousehold_income_beforetax, Employment_status, BMI_i0,
         Smoking_status, Alcohol_frequency, Healthy_diet_score, chronic_conditions,
         Overall_health_rating, SB_Meantime, LightB_Meantime, Standard_PRS,
         MVPA_Meantime, incident_VTE) %>%
  tbl_summary(by = incident_VTE, type = all_continuous() ~ "continuous2", missing = "ifany") %>%
  add_p()


table1

# 将tbl_summary对象转换为tibble
table1_tibble <- as_tibble(table1)

# 将结果写入CSV文件
write.csv(table1_tibble, "C:/Users/yougch/Desktop/baseline/table1.csv")




##sTable基线
library(gtsummary)

# 创建基线表格
table1 <- data_impute %>%
  select(Age, Sex, Ethnicity, Education, Townsend_deprivation_index,
         Average_totalhousehold_income_beforetax, Employment_status, BMI_i0,
         Smoking_status, Alcohol_frequency, Healthy_diet_score, chronic_conditions,
         Overall_health_rating, SB_Meantime, LightB_Meantime, Standard_PRS,
         MVPA_Meantime, MVPA_meantime_quartile) %>%
  tbl_summary(by = MVPA_meantime_quartile, type = all_continuous() ~ "continuous2", missing = "ifany") %>%
  add_p()

table1

# 将tbl_summary对象转换为tibble
table1_tibble <- as_tibble(table1)

# 将结果写入CSV文件
write.csv(table1_tibble, "C:/Users/yougch/Desktop/stable/table1.csv")



#WW_sTable_guideline 基线
library(gtsummary)

# 创建基线表格
table1 <- data_ww %>%
  select(Age, Sex, Ethnicity, Education, Townsend_deprivation_index,
         Average_totalhousehold_income_beforetax, Employment_status, BMI_i0,
         Smoking_status, Alcohol_frequency, Healthy_diet_score, chronic_conditions,
         Overall_health_rating, SB_Meantime, LightB_Meantime, Standard_PRS,
         MVPA_Meantime, activity_group) %>%
  tbl_summary(by = activity_group, type = all_continuous() ~ "continuous2", missing = "ifany") %>%
  add_p()

table1

# 将tbl_summary对象转换为tibble
table1_tibble <- as_tibble(table1)

# 将结果写入CSV文件
write.csv(table1_tibble, "C:/Users/yougch/Desktop/WW_stable/table1.csv")


#WW_sTable_median 基线
library(gtsummary)

# 创建基线表格
table1 <- data_ww %>%
  select(Age, Sex, Ethnicity, Education, Townsend_deprivation_index,
         Average_totalhousehold_income_beforetax, Employment_status, BMI_i0,
         Smoking_status, Alcohol_frequency, Healthy_diet_score, chronic_conditions,
         Overall_health_rating, SB_Meantime, LightB_Meantime, Standard_PRS,
         MVPA_Meantime, activity_group_median) %>%
  tbl_summary(by = activity_group_median, type = all_continuous() ~ "continuous2", missing = "ifany") %>%
  add_p()

table1

# 将tbl_summary对象转换为tibble
table1_tibble <- as_tibble(table1)

# 将结果写入CSV文件
write.csv(table1_tibble, "C:/Users/yougch/Desktop/WW_median/table1.csv")


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

test_ph <- cox.zph(mod_vte_min)
test_ph


# 基于Q1作为参考组建立 Cox 比例风险回归模型:全变量
mod_vte <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_meantime_quartile, data = data_impute)

summary(mod_vte)


test_ph <- cox.zph(mod_vte)
test_ph

#ggcoxzph(test_ph,resid=T,se=T,point.col = "red",point.size = 1) 


#PRS与MVPA相乘交互作用的计算
#The interaction between PRS categories and MVPA categories were exaimed by liklihood ratio test

mod_vte0 <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                    Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                    Overall_health_rating + chronic_conditions + Healthy_diet_score + LightB_Meantime + SB_Meantime +
                    genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                    genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                    PRS_tertiles + MVPA_meantime_quartile, data = data_impute)

mod_vte1 <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                    Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                    Overall_health_rating + chronic_conditions + Healthy_diet_score + LightB_Meantime + SB_Meantime +
                    genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                    genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                    PRS_tertiles*MVPA_meantime_quartile, data = data_impute)

summary(mod_vte0)
summary(mod_vte1)
anova(mod_vte0)
anova(mod_vte1)

library(lmtest)
lrtest(mod_vte1, mod_vte0)


test_ph <- cox.zph(mod_vte1)
test_ph





##

#Association between PRS and incident VTE_Supplementary Table

# 使用quantile函数进行三等分割
tertiles <- quantile(data_impute$Standard_PRS, probs = c(1/3, 2/3))
tertiles

# 将数据分割并依次命名为Q1, Q2, Q3
data_impute$PRS_tertiles <- cut(data_impute$Standard_PRS, breaks = c(-Inf, tertiles, Inf), labels = c("Q1", "Q2", "Q3"), include.lowest = TRUE)

# 打印结果
head(data_impute$PRS_tertiles)


#Minimally 
mod_prs_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + 
                       PRS_tertiles, data = data_impute)

summary(mod_prs_min)

test_ph <- cox.zph(mod_prs_min)
test_ph

# Full variables
mod_prs_full <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                   Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                   Overall_health_rating + chronic_conditions + Healthy_diet_score  + LightB_Meantime + SB_Meantime + 
                   genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                   genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                   MVPA_Meantime + PRS_tertiles, data = data_impute)

summary(mod_prs_full)

test_ph <- cox.zph(mod_prs_full)
test_ph

#
summary(data_impute$PRS_tertiles)

Q1_cases_VTE <- sum(data_impute$incident_VTE[data_impute$PRS_tertiles == "Q1" & data_impute$incident_VTE == 1])
Q1_cases_VTE

Q2_cases_VTE <- sum(data_impute$incident_VTE[data_impute$PRS_tertiles == "Q2" & data_impute$incident_VTE == 1])
Q2_cases_VTE

Q3_cases_VTE <- sum(data_impute$incident_VTE[data_impute$PRS_tertiles == "Q3" & data_impute$incident_VTE == 1])
Q3_cases_VTE

#
#Table 3
data_highPRS <- data_impute[data_impute$PRS_tertiles == "Q3", ]

mod_highPRS <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                       Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                       Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime +SB_Meantime + 
                       genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                       genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                       MVPA_meantime_quartile , data = data_highPRS)

summary(mod_highPRS)

test_ph <- cox.zph(mod_highPRS)
test_ph

#0/1
cases_highPRS <- sum(data_impute$PRS_tertiles == "Q3" & data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1)
cases_highPRS

cases_highPRS <- sum(data_impute$PRS_tertiles == "Q3" & data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1)
cases_highPRS

cases_highPRS <- sum(data_impute$PRS_tertiles == "Q3" & data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1)
cases_highPRS

cases_highPRS <- sum(data_impute$PRS_tertiles == "Q3" & data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1)
cases_highPRS



#
data_mediumPRS <- data_impute[data_impute$PRS_tertiles == "Q2", ]

mod_mediumPRS <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS  + LightB_Meantime +SB_Meantime +
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         MVPA_meantime_quartile , data = data_mediumPRS)

summary(mod_mediumPRS)

test_ph <- cox.zph(mod_mediumPRS)
test_ph

#
cases_mediumPRS <- sum(data_impute$PRS_tertiles == "Q2" & data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1)
cases_mediumPRS

cases_mediumPRS <- sum(data_impute$PRS_tertiles == "Q2" & data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1)
cases_mediumPRS

cases_mediumPRS <- sum(data_impute$PRS_tertiles == "Q2" & data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1)
cases_mediumPRS

cases_mediumPRS <- sum(data_impute$PRS_tertiles == "Q2" & data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1)
cases_mediumPRS



#
data_lowPRS <- data_impute[data_impute$PRS_tertiles == "Q1", ]

mod_lowPRS <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                         Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                         Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                         genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                         genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                         MVPA_meantime_quartile , data = data_lowPRS)

summary(mod_lowPRS)

test_ph <- cox.zph(mod_lowPRS)
test_ph

cases_lowPRS <- sum(data_impute$PRS_tertiles == "Q1" & data_impute$MVPA_meantime_quartile == "Q1" & data_impute$incident_VTE == 1)
cases_lowPRS

cases_lowPRS <- sum(data_impute$PRS_tertiles == "Q1" & data_impute$MVPA_meantime_quartile == "Q2" & data_impute$incident_VTE == 1)
cases_lowPRS

cases_lowPRS <- sum(data_impute$PRS_tertiles == "Q1" & data_impute$MVPA_meantime_quartile == "Q3" & data_impute$incident_VTE == 1)
cases_lowPRS

cases_lowPRS <- sum(data_impute$PRS_tertiles == "Q1" & data_impute$MVPA_meantime_quartile == "Q4" & data_impute$incident_VTE == 1)
cases_lowPRS




####Figure 3

library(dplyr)

data_impute <- data_impute %>%
  mutate(JointPRS = case_when(
    PRS_tertiles == "Q3" & MVPA_meantime_quartile == "Q1" ~ "J1",
    PRS_tertiles == "Q3" & MVPA_meantime_quartile == "Q2" ~ "J2",
    PRS_tertiles == "Q3" & MVPA_meantime_quartile == "Q3" ~ "J3",
    PRS_tertiles == "Q3" & MVPA_meantime_quartile == "Q4" ~ "J4",
    PRS_tertiles == "Q2" & MVPA_meantime_quartile == "Q1" ~ "J5",
    PRS_tertiles == "Q2" & MVPA_meantime_quartile == "Q2" ~ "J6",
    PRS_tertiles == "Q2" & MVPA_meantime_quartile == "Q3" ~ "J7",
    PRS_tertiles == "Q2" & MVPA_meantime_quartile == "Q4" ~ "J8",
    PRS_tertiles == "Q1" & MVPA_meantime_quartile == "Q1" ~ "J9",
    PRS_tertiles == "Q1" & MVPA_meantime_quartile == "Q2" ~ "J10",
    PRS_tertiles == "Q1" & MVPA_meantime_quartile == "Q3" ~ "J11",
    PRS_tertiles == "Q1" & MVPA_meantime_quartile == "Q4" ~ "J12",
    TRUE ~ NA_character_
  ))

mod_jointPRS <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                      Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                      Overall_health_rating + chronic_conditions + Healthy_diet_score  + LightB_Meantime + SB_Meantime + 
                      genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                      genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                      JointPRS , data = data_impute)

summary(mod_jointPRS)

test_ph <- cox.zph(mod_jointPRS)
test_ph


#
library(ggplot2)
library(dplyr)

# 将数据整理成适合绘图的格式
figure3 <- data.frame(
  JointPRS = factor(c("J1","J2", "J3", "J4", "J5", "J6", "J7", "J8", "J9", "J10", "J11", "J12"), levels=c("J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8", "J9", "J10", "J11", "J12")),
  HR = c(1, 0.7999, 0.8029, 0.7208, 0.7123, 0.6305, 0.4980, 0.4419, 0.6563,0.4357, 0.3791, 0.1687),
  lower = c(NA, 0.59277, 0.58619, 0.46359, 0.48426, 0.46179, 0.35525, 0.26171, 0.43897,0.31219, 0.26540, 0.07966),
  upper = c(NA, 1.0795, 1.0996, 1.1208, 1.0477, 0.8609, 0.6980, 0.7462, 0.9811,0.6080, 0.5414, 0.3572)
)


# 创建图形
ggplot(figure3, aes(x = JointPRS, y = HR, ymin = lower, ymax = upper)) +
  geom_point(size = 2.3, shape=16) +  
  geom_errorbar(width = 0.07) +
  labs(title = "MVPA, h/d",
       x = "",
       y = "Adjusted HR (95% CI) for VTE risk")+
  theme_minimal() +
  theme(
    axis.text = element_text(size=13), 
    axis.title = element_text(size=14), 
    plot.title = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_line(size = 0.2), 
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.y = element_line(size = 0.2)  
  ) +
  scale_y_continuous(breaks = c(seq(0, 0.5, 0.25), seq(0.5, 1, 0.1), 1.1), 
                     labels = function(x) {
                       ifelse(x %in% c(0.25, 0.5, 1.0, 1.1), 
                              ifelse(x %% 1 == 0, sprintf("%.2f", x), sprintf("%.2f", x)), 
                              "")
                     })+
  geom_hline(yintercept = 1, linetype="dashed", color = "black")




###############################################################################################
###############################################################################################
###############################################################################################

#data_ww load
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


data_ww <- read.csv("E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/data_mvpa.csv", header = TRUE)
summary(data_ww)

names(data_ww)

data_ww$SB_Meantime <- 24*data_ww$Sedentary_overall_average
summary(data_ww$SB_Meantime)

data_ww$LightB_Meantime <- 24*data_ww$Light_overall_average
summary(data_ww$LightB_Meantime)

data_ww$MVPA_Meantime <- 24*data_ww$MVPA_overall_average
summary(data_ww$MVPA_Meantime)

data_ww$Sex <- as.factor(data_ww$Sex)
data_ww$chronic_conditions <- as.factor(data_ww$chronic_conditions)

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


mod_ww_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex +
                  activity_group, data = data_ww)

summary(mod_ww_min)

mod_ww <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime + 
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group, data = data_ww)

summary(mod_ww)

test_ph <- cox.zph(mod_ww)
test_ph

#
Q1_group <- data_ww$time_accel_to_VTE[data_ww$activity_group == "Inactive"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_ww$time_accel_to_VTE[data_ww$activity_group == "Active - Regular"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_ww$time_accel_to_VTE[data_ww$activity_group == "Active - WW"]
sum_Q3 <- sum(Q3_group)
sum_Q3



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

#
mod_ww_median_min <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex +
                         activity_group_median, data = data_ww)

summary(mod_ww_median_min)

mod_ww_median <- coxph(Surv(time_accel_to_VTE, incident_VTE) ~ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                  Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                  Overall_health_rating + chronic_conditions + Healthy_diet_score + Standard_PRS + LightB_Meantime + SB_Meantime +
                  genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 + genetic_a4 + genetic_a5 + genetic_a6 +
                  genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10 +
                  activity_group_median, data = data_ww)

summary(mod_ww_median)

test_ph <- cox.zph(mod_ww_median)
test_ph

#
Q1_group <- data_ww$time_accel_to_VTE[data_ww$activity_group_median == "Inactive"]
sum_Q1 <- sum(Q1_group)
sum_Q1

Q2_group <- data_ww$time_accel_to_VTE[data_ww$activity_group_median == "Active - Regular"]
sum_Q2 <- sum(Q2_group)
sum_Q2

Q3_group <- data_ww$time_accel_to_VTE[data_ww$activity_group_median == "Active - WW"]
sum_Q3 <- sum(Q3_group)
sum_Q3


#
inactive_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Inactive" & data_ww$incident_VTE == 1])
inactive_cases_VTE

activeww_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - WW" & data_ww$incident_VTE == 1])
activeww_cases_VTE

activeregular_cases_VTE <- sum(data_ww$incident_VTE[data_ww$activity_group_median == "Active - Regular" & data_ww$incident_VTE == 1])
activeregular_cases_VTE



########################################################################################密度图1

###150min =2.5h,ww
data_activeww <- subset(data_ww, activity_group == "Active - WW")
# 选择mvpa1-7列并进行大小比较
nrow(data_activeww)

library(dplyr)
#
data_activeww <- data_activeww %>%
  mutate(top2_mvpa = (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                        pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                             pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                             pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                             pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                             pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))) / 2)


data_activeww <- data_activeww %>%
  mutate(remain5_mvpa = (mvpa_daily_total- (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                        pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                             pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                             pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                             pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                             pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))))/ 5)
head(data_activeww)

summary(data_activeww$top2_mvpa)


#
ggplot(data_activeww, aes(x = top2_mvpa, fill = "Top 2 d")) + 
  geom_density(alpha = 0.5) +
  geom_density(aes(x = remain5_mvpa, fill = "Remaining 5 d"), alpha = 0.5) +
  labs(title = "Weekend warrior activity",
       x = "MVPA, h/d",
       y = "Density") +
  scale_fill_manual(values = c("#f18e0c","#3C6F80"), guide = "none") +
  coord_cartesian(xlim = c(0,4)) +
  scale_x_continuous(breaks = c(0, 1, 2, 4)) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        plot.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"))





###150min =2.5h,regular
data_activeregular <- subset(data_ww, activity_group == "Active - Regular")
# 选择mvpa1-7列并进行大小比较
nrow(data_activeregular)

library(dplyr)
#
data_activeregular <- data_activeregular %>%
  mutate(top2_mvpa = (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                        pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                             pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                             pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                             pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                             pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))) / 2)


data_activeregular <- data_activeregular %>%
  mutate(remain5_mvpa = (mvpa_daily_total- (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                                              pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                                                   pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                                                   pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                                                   pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                                                   pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))))/ 5)
head(data_activeregular)

summary(data_activeregular$top2_mvpa)



#
ggplot(data_activeregular, aes(x = top2_mvpa, fill = "Top 2 d")) + 
  geom_density(alpha = 0.5) +
  geom_density(aes(x = remain5_mvpa, fill = "Remaining 5 d"), alpha = 0.5) +
  labs(title = "Regular activity",
       x = "MVPA, h/d",
       y = "Density") +
  scale_fill_manual(values = c("#f18e0c","#3C6F80")) +
  coord_cartesian(xlim = c(0,4), ylim = c(0, 2)) +
  scale_x_continuous(breaks = c(0, 1, 2, 4)) +
  theme_minimal() +
  theme(legend.position = c(0.90, 0.90),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank())

########################################################################################密度图2



###
###Median:3.84h, ww
data_activeww <- subset(data_ww, activity_group_median == "Active - WW")
# 选择mvpa1-7列并进行大小比较
nrow(data_activeww)

library(dplyr)
#
data_activeww <- data_activeww %>%
  mutate(top2_mvpa = (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                        pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                             pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                             pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                             pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                             pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))) / 2)


data_activeww <- data_activeww %>%
  mutate(remain5_mvpa = (mvpa_daily_total- (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                                              pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                                                   pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                                                   pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                                                   pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                                                   pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))))/ 5)
head(data_activeww)

summary(data_activeww$top2_mvpa)


#
ggplot(data_activeww, aes(x = top2_mvpa, fill = "Top 2 d")) + 
  geom_density(alpha = 0.5) +
  geom_density(aes(x = remain5_mvpa, fill = "Remaining 5 d"), alpha = 0.5) +
  labs(title = "Weekend warrior activity",
       x = "MVPA, h/d",
       y = "Density") +
  scale_fill_manual(values = c('#D93D25', "#1077A9"), guide = "none") +
  coord_cartesian(xlim = c(0,4)) +
  scale_x_continuous(breaks = c(0, 1, 2, 4)) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text = element_text(size=13),
        axis.title = element_text(size=14),
        plot.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"))



###Median:3.84h,regular
data_activeregular <- subset(data_ww, activity_group_median == "Active - Regular")
# 选择mvpa1-7列并进行大小比较
nrow(data_activeregular)

library(dplyr)
#
data_activeregular <- data_activeregular %>%
  mutate(top2_mvpa = (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                        pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                             pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                             pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                             pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                             pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))) / 2)


data_activeregular <- data_activeregular %>%
  mutate(remain5_mvpa = (mvpa_daily_total- (pmax(mvpa1, mvpa2, mvpa3, mvpa4, mvpa5, mvpa6, mvpa7) +
                                              pmax(pmin(mvpa1, mvpa2), pmin(pmax(mvpa1, mvpa2), mvpa3), 
                                                   pmin(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), 
                                                   pmin(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), 
                                                   pmin(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), 
                                                   pmin(pmax(pmax(pmax(pmax(pmax(mvpa1, mvpa2), mvpa3), mvpa4), mvpa5), mvpa6), mvpa7))))/ 5)
head(data_activeregular)

summary(data_activeregular$top2_mvpa)



#
ggplot(data_activeregular, aes(x = top2_mvpa, fill = "Top 2 d")) + 
  geom_density(alpha = 0.5) +
  geom_density(aes(x = remain5_mvpa, fill = "Remaining 5 d"), alpha = 0.5) +
  labs(title = "Regular activity",
       x = "MVPA, h/d",
       y = "Density") +
  scale_fill_manual(values = c('#D93D25', "#1077A9")) +
  coord_cartesian(xlim = c(0,4), ylim = c(0, 2)) +
  scale_x_continuous(breaks = c(0, 1, 2, 4)) +
  theme_minimal() +
  theme(legend.position = c(0.90, 0.90),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.title = element_blank())




#Forest plot 图太丑，先不做了，做Table 4
library(forestplot)
library(ggplot2)
rs_forest <- read.csv('E:/CVDs_Team/UKB/Acceleratory_analysis/20231212_MVPA/Data/plot/Forestplot_WW.csv',header = FALSE)
rs_forest
rs_forest$V6 <- as.numeric(rs_forest$V6)
rs_forest$V7 <- as.numeric(rs_forest$V7)
rs_forest$V8 <- as.numeric(rs_forest$V8)


# 绘制森林图

forestplot(labeltext = as.matrix(rs_forest[,c(1:5)]),
           mean = rs_forest$V6,    
           lower = rs_forest$V7,      
           upper = rs_forest$V8, 
           is.summary = c(TRUE, F, FALSE, FALSE, FALSE, F, FALSE, FALSE, FALSE),
           align='l',
           zero = 1,
           xlab = "Hazard ratio (95% CI)",
           col = fpColors(box = "#8c232a", line = "black", summary = "royalblue"),
           graphwidth = unit(.20,"npc"),
           boxsize = 0.1,
           colgap=unit(9, 'mm'),
           #lineheight=unit(2, 'mm'), 
           #line.margin=unit(2, 'mm'),
           txt_gp=fpTxtGp(label=gpar(cex=1.25,fontface = "plain"), ticks=gpar(cex=1.15), xlab=gpar(cex = 1.15), title=gpar(cex = 1.15)),
           graph.pos=6,
           xticks=c(0.5,0.6,0.7,0.8,0.9,1.0,1.5)
)



###绘制cumulative risk曲线,guideline-based
data_ww$activity_group1 <- data_ww$activity_group
prodlim_vte <- prodlim(Hist(time_accel_to_VTE,incident_VTE)~activity_group1,data=data_ww)

#CairoPDF(file='C:/Users/yougch/Desktop/11/km_guideline_vte.pdf',height=3,width=3.5,
#         pointsize=4)

#par(oma=c(3,1,1,1),mar=c(3,1,1,1))

plot(prodlim_vte,"cuminc",ylim=c(0,0.04),xlim=c(0,7),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=seq(0,7,1),axis1.labels=as.character(0:7),
     atrisk.times=seq(0,7,1),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=1.4,axis2.cex.axis=1.4,axis1.padj=0.3,
     legend.cex=1.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),xlab='',ylab='')
     #atrisk.title=("                     "),atrisk.pos=-2,#atrisk.line=c(1.2,2.8,4.4),
     #atrisk.cex=1.8,atrisk.interspace=1.8,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=2.5,at=0.02,cex=1.6)
mtext("Years",side=1, line=-4.5,cex=1.6)
#mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)

#dev.off()


###Median-based
data_ww$activity_group2 <- data_ww$activity_group_median
prodlim_vte_median <- prodlim(Hist(time_accel_to_VTE,incident_VTE)~activity_group2,data=data_ww)
plot(prodlim_vte_median,"cuminc",ylim=c(0,0.04),xlim=c(0,7),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=seq(0,7,1),axis1.labels=as.character(0:7),
     atrisk.times=seq(0,7,1),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=1.4,axis2.cex.axis=1.4,axis1.padj=0.3,
     legend.cex=1.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),xlab='',ylab='')
#atrisk.title=("                     "),atrisk.pos=-2,#atrisk.line=c(1.2,2.8,4.4),
#atrisk.cex=1.8,atrisk.interspace=1.8,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=2.5,at=0.02,cex=1.6)
mtext("Years",side=1, line=-4.5,cex=1.6)

data_ww <- data_ww[, !names(data_ww) %in% c("activity_group1", "activity_group2")]














###########################################################################################################################################################
###########################################################################################################################################################
#sFigure 1

# 绘制直方图


P1 <- ggplot(data_impute, aes(x=Standard_PRS)) +
  geom_histogram(binwidth=0.2, fill = "grey" ,color="black") +
  labs(title="Distribution of standard PRS", x="Standardized PRS", y="No. of participants") +
  coord_cartesian(xlim = c(-2.6, 2.6)) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text = element_text(size=13), #坐标轴上的文本的字体大小
        axis.title = element_text(size=14), # 坐标轴标题的字体大小
        plot.title = element_text(size = 16),  # 设置标题字体大小
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"))+
        scale_x_continuous(breaks = c(seq(-2, 2, 1),-2.5,2.5), 
        labels = c(seq(-2, 2, 1),"-2.5","2.5"))
P1



##########################################################################################
#绘制RCS曲线
#RCS_fully_adjusted
library(rms)
dd1 <- datadist(data_impute)
options(datadist="dd1")
ref_value <- median(data_impute$Standard_PRS, na.rm = TRUE)
dd1$limits$Standard_PRS[3] <- ref_value

fit.cox <- cph(Surv(time_accel_to_VTE, incident_VTE==1)~ rcs(Standard_PRS,3)+ Age + Sex + Ethnicity + Education + Townsend_deprivation_index +
                 Average_totalhousehold_income_beforetax + Employment_status + BMI_i0 + Smoking_status + Alcohol_frequency + 
                 Overall_health_rating + chronic_conditions + Healthy_diet_score + LightB_Meantime + MVPA_Meantime +
                 SB_Meantime + genotyping_batch + genetic_a1 +genetic_a2 +genetic_a3 +genetic_a4 + genetic_a5 + genetic_a6 +
                 genetic_a7 + genetic_a8 + genetic_a9 + genetic_a10, data=data_impute, x=TRUE, y=TRUE)

HR <- Predict(fit.cox, Standard_PRS, fun = exp, ref.zero = TRUE)
anova(fit.cox)

P2 <- ggplot() +
  geom_line(data=HR, aes(Standard_PRS, yhat),linetype="solid",linewidth=1,alpha = 0.7,colour="red") +
  geom_ribbon(data=HR, aes(Standard_PRS, ymin = lower, ymax = upper),alpha = 0.1,fill="dimgrey") +
  theme_classic() +
  geom_hline(yintercept=1, linetype=2,size=0.5) +
  xlab('Standardized PRS') +
  theme(axis.title.x = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ylab('HR (95% CI) for VTE risk') +
  theme(axis.title.y = element_text(size=14,vjust=0.5,hjust=0.5)) +
  ggtitle('Fully adjusted association model') +
  theme(plot.title = element_text(size = 16, vjust = 0.5, hjust = 0),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.ticks.length = unit(0.2, "cm")) +  # 设置坐标轴刻度线长度为0.2
  coord_cartesian(xlim = c(-2.6, 2.6),ylim=c(0.24,3)) +
  scale_y_continuous(breaks = c(0.25,0.5,1,2,3), 
                     labels = c("0.25","0.5","1","2","3"))+
  scale_x_continuous(breaks = c(seq(-2.0, 2.0, 1.0),-2.5,2.5), 
                     labels = c(seq(-2.0, 2.0, 1.0),"-2.5","2.5"))

P2






###########################################################################################################################################################
###########################################################################################################################################################
#sFigure 2


high_group <- data_highPRS$time_accel_to_VTE[data_highPRS$MVPA_meantime_quartile == "Q4"]
sum_Q4 <- sum(Q4_group)
sum_Q4

#
J1_cases_VTE <- sum(data_highPRS$MVPA_meantime_quartile == "Q1" & data_highPRS$incident_VTE == 1)
J1_cases_VTE

J2_cases_VTE <- sum(data_highPRS$MVPA_meantime_quartile == "Q2" & data_highPRS$incident_VTE == 1)
J2_cases_VTE

J3_cases_VTE <- sum(data_highPRS$MVPA_meantime_quartile == "Q3" & data_highPRS$incident_VTE == 1)
J3_cases_VTE

J4_cases_VTE <- sum(data_highPRS$MVPA_meantime_quartile == "Q4" & data_highPRS$incident_VTE == 1)
J4_cases_VTE


J1_person_time <- data_highPRS$time_accel_to_VTE[data_highPRS$MVPA_meantime_quartile == "Q1"]
sum_J1 <- sum(J1_person_time)
sum_J1

J2_person_time <- data_highPRS$time_accel_to_VTE[data_highPRS$MVPA_meantime_quartile == "Q2"]
sum_J2 <- sum(J2_person_time)
sum_J2

J3_person_time <- data_highPRS$time_accel_to_VTE[data_highPRS$MVPA_meantime_quartile == "Q3"]
sum_J3 <- sum(J3_person_time)
sum_J3

J4_person_time <- data_highPRS$time_accel_to_VTE[data_highPRS$MVPA_meantime_quartile == "Q4"]
sum_J4 <- sum(J4_person_time)
sum_J4





J12_cases_VTE <- sum(data_lowPRS$MVPA_meantime_quartile == "Q4" & data_lowPRS$incident_VTE == 1)
J12_cases_VTE


J12_person_time <- data_lowPRS$time_accel_to_VTE[data_lowPRS$MVPA_meantime_quartile == "Q4"]
sum_J12 <- sum(J12_person_time)
sum_J12




