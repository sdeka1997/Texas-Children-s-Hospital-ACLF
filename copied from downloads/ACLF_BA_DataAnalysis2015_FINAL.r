# Investigating acute-on-chronic liver failure on Biliary Atresia pediatric
# patients on the liver transplant list

# Charles Ho
# *in conjunction with STAT 450 project*
#
# April 2015
#
# Summary of Data Analysis
# 
# 1. Preliminary data setup and calculations
# 2. Significance testing
# 3. Building a risk assessment model for ACLF development
#
####### 1 - PRELIMINARY DATA SETUP  #######
# reading data
setwd("C:/Users/CharlesHo/Desktop/Spring 2015/STAT 450")
library(BMA)
library(pROC)

data <- read.csv("ba_data_2015_v4.csv")
ekg_data <- read.csv("ekg_data_2015_modified_pretrans.csv")
echo_data <- read.csv("echo_data_v1.csv")

# computing peld and  adjusted balf score
data$listing_peld <- 10*(0.480*log(data$total_bili)+1.857*log(data$inr)
                         -0.687*log(data$albumin)+0.463*data$age.factor)
data[data$listing_peld <0,]$listing_peld <- 0   

data$listing_balf <- (7.196+1.438*log(data$total_bili)+0.434*log(data$ggt)
                      -3.491*log(data$albumin)-0.670*log(data$listing_age/365))
data$listing_balf[107] <- 0   

# data$adj_balf <- (data$listing_balf-4)*3   #compute "adjusted" balf score

# Combine BA and EKG data
ba_ekg_data <- merge(ekg_data,data, by = "mrn")
ba_ekg_data$delta <- ba_ekg_data$qtc - ba_ekg_data$qtc_norm

# Combine BA and Echo data
ba_echo_data <- merge(data,echo_data, by = "mrn")

####### 2 - HYPOTHESIS TESTING OF SIGNIFICANT DIFFERENCES #######
analysis1 <- function(var, data1, data2){
  a <- subset(data1, select = var)
  a1 <- a[,1]
  #print(median(a1, na.rm = T))
  #print(mean(a1, na.rm = T))
  print(summary(a1))
  print(sd(a1, na.rm = T))
  
  
  b <- subset(data2, select = var)
  b1 <- b[,1]
  #print(median(b1, na.rm = T))
  #print(mean(b1, na.rm = T))
  print(summary(b1))
  print(sd(b1, na.rm = T))
  
  print(t.test(a1,b1))
  print(wilcox.test(a1,b1))
  
  print(shapiro.test(a1))
  print(shapiro.test(b1))
  
}

# Scenario 1: BA ACLF vs. non-ACLF groups
aclf <- subset(data, aclf_status == "ACLF")
non_aclf <- subset(data, aclf_status == "non-ACLF")

# Scenario 2: BA ACLF vs. non-ACLF (max wait time) groups
# stratified analysis to exclude non-ACLF patients with short waitlist times
# max(aclf$wait)    # longest waiting time for an ACLF patient was 168 days
# sum(non_aclf$wait > 167)   #25 non-ACLF patients remains
non_aclf <- subset(non_aclf, wait > 167)

# Scenario 3: EKG ACLF vs. non-ACLF groups
ekg_aclf <- subset(ekg_data, aclf_status == "ACLF")
ekg_nonaclf <- subset(ekg_data, aclf_status == "non-ACLF")

ekg_aclf <- subset(ba_ekg_data, aclf_status.x == "ACLF")
ekg_nonaclf <- subset(ba_ekg_data, aclf_status.x == "non-ACLF")

# Scenario 4: Echo ACLF vs. non-ACLF groups
aclf <- subset(ba_echo_data, aclf_status == "ACLF")
non_aclf <- subset(ba_echo_data, aclf_status == "non-ACLF")


# Propotions testing

# analyze qtc and gf
summary(ekg_data$qtc)
summary(ekg_data$qtc > 445)
summary(ekg_data$qtc >= 445 & ekg_data$aclf_status == "ACLF")
summary(ekg_data$qtc < 445 & ekg_data$aclf_status == "ACLF")

data$gf <- 0
data$gf[which(data$z_wt.age < -2)] <- 1
data$gf[which(data$z_len.age < -2)] <- 1
data$gf <- as.factor(data$gf)
data$nogf <- 0
data$nogf[which(data$z_wt.age > 0)] <- 1
data$nogf[which(data$z_len.age > 0)] <- 1
data$nogf <- as.factor(data$nogf)
summary(data$gf)
summary(data$gf == 1 & data$aclf_status == 1)
summary(data$gf == 0 & data$aclf_status == 1)

sum(aclf$kasai_status == 1)   # 11 ACLF had kasai
sum(aclf$kasai_status == 0)   # 6 ACLF had no kasai
sum(non_aclf$kasai_status == 1, na.rm = TRUE)  # 50 non-ACLF had kasai
sum(non_aclf$kasai_status == 0, na.rm = TRUE)  # 16 non-ACLF had no kasai

prop.test(c(11,50), c(17,66))
fisher.test(cbind(c(11,50),c(6,16)))        # p = 0.37 for Kasai and ACLF

sum(aclf$kasai_age > 42, na.rm = TRUE)  # 9 ACLF kasais were late
sum(aclf$kasai_age < 42, na.rm = TRUE)  # 1 ACLF kasais were not late
sum(non_aclf$kasai_age > 42, na.rm = TRUE) # 32 non-ACLF kasais were late
sum(non_aclf$kasai_age < 42, na.rm = TRUE) # 6 non-ACLF kasais were not late

prop.test(c(9,32), c(10,38))
fisher.test(cbind(c(9,32),c(1,6)))     # p = 1


# looking at only those 1 yr or less for selected echo parameters
young <- subset(ba_echo_data, listing_age < 366)
young$ea_ratio <- young$tv_e/young$tv_a
aclf <- subset(young, aclf_status == "ACLF")
non_aclf <- subset(young, aclf_status == "non-ACLF")

aclf <- subset(ba_echo_data, aclf_status == "ACLF")
non_aclf <- subset(ba_echo_data, aclf_status == "non-ACLF")
analysis1("e", aclf, non_aclf)

# odds ratios calculations for selected parameters
levels(data$aclf_status) = c("ACLF" = 1, "non-ACLF" = 0)   
data$aclf_status <- factor(data$aclf_status, levels = c("0", "1"))
glm1 <- glm(aclf_status ~ total_bili, data = data, family = "binomial")
glm2 <- glm(aclf_status ~ inr, data = data, family = "binomial")
glm3 <- glm(aclf_status ~ platelets, data = data, family = "binomial")
glm4 <- glm(aclf_status ~ ggt, data = data, family = "binomial")
glm5 <- glm(aclf_status ~ z_bmi.age, data = data, family = "binomial")
levels(ba_echo_data$aclf_status) = c("ACLF" = 1, "non-ACLF" = 0)   
ba_echo_data$aclf_status <- factor(ba_echo_data$aclf_status, levels = c("0", "1"))
glm6 <- glm(aclf_status ~ mmlv_mass_z, data = ba_echo_data, family = "binomial")
glm7 <- glm(aclf_status ~ lv_midwall_sys_z, data = ba_echo_data, family = "binomial")

roc_tb <- roc(data$aclf_status, data$total_bili)
coords(roc_tb, "best", best.method = "youden")
summary(aclf$total_bili < 7.9)
summary(non_aclf$total_bili < 7.9)
table_tb <- cbind(c(18,2), c(37,57))
fisher.test(table_tb)

roc_inr <- roc(data$aclf_status, data$inr)
coords(roc_inr, "best", best.method = "youden")
summary(aclf$inr < 1.55)
summary(non_aclf$inr < 1.55)
table_inr <- cbind(c(14,6), c(16,78))
fisher.test(table_inr)

roc_bmi <- roc(data$aclf_status, data$z_bmi.age)
coords(roc_bmi, "best", best.method = "youden")
summary(aclf$z_bmi.age < -0.81)
summary(non_aclf$z_bmi.age < -0.81)
table_bmi <- cbind(c(9,11),c(16,69))
fisher.test(table_bmi)

roc_platelets <- roc(data$aclf_status, data$platelets)
coords(roc_platelets, "best", best.method = "youden")
summary(aclf$platelets < 162.5)
summary(non_aclf$platelets < 162.5)
table_platelets <- cbind(c(14,5),c(48,38))
fisher.test(table_platelets)


roc_ggt <- roc(data$aclf_status, data$ggt)
coords(roc_ggt, "best", best.method = "youden")
summary(aclf$ggt < 379.5)
summary(non_aclf$ggt < 379.5)
table_ggt <- cbind(c(17,3),c(52,36))
fisher.test(table_ggt)

roc_mass <- roc(ba_echo_data$aclf_status, ba_echo_data$mmlv_mass_z)
coords(roc_mass, "best", best.method = "youden")
summary(aclf$mmlv_mass_z < 1.90)
summary(non_aclf$mmlv_mass_z < 1.90)
table_mass <- cbind(c(7,11), c(14,61))
fisher.test(table_mass)

roc_sys <- roc(ba_echo_data$aclf_status, ba_echo_data$lv_midwall_sys_z)
coords(roc_sys, "best", best.method = "youden")
summary(aclf$lv_midwall_sys_z < 0.605)
summary(non_aclf$lv_midwall_sys_z < 0.605)
table_sys <- cbind(c(11,2), c(15,37))
fisher.test(table_sys)

####### 3 - MODEL BUILDING  #######

### Part A: Model on only BA+EKG data #####
levels(data$aclf_status) = c("ACLF" = 1, "non-ACLF" = 0)   
data$aclf_status <- factor(data$aclf_status, levels = c("0", "1"))

# data_rmna is for BA data only
# subset data to exclude missing values for GGT, platemets z-scores, sodium 
data_rmna <- data[!is.na(data$ggt),]
data_rmna <- data_rmna[!is.na(data_rmna$platelets),]  
data_rmna <- data_rmna[!is.na(data_rmna$z_wt.len),]
data_rmna <- data_rmna[!is.na(data_rmna$sodium),]
#data_rmna <- subset(data_rmna, listing_age < 1000)
#data_rmna <- data[!is.na(data$ggt),]
#data_rmna <- data_rmna[!is.na(data_rmna$z_wt.len),]
data_rmna$gf <- 0
data_rmna$gf[which(data_rmna$z_wt.age < -2)] <- 1
data_rmna$gf[which(data_rmna$z_len.age < -2)] <- 1

# data_rmna_2 is for BA + EKG combined dataset
data_rmna_2 <- ba_ekg_data[!is.na(ba_ekg_data$ggt),]     # 98 -> 94
data_rmna_2 <- data_rmna_2[!is.na(data_rmna_2$platelets),]  # 94 --> 91
data_rmna_2 <- data_rmna_2[!is.na(data_rmna_2$z_wt.len),]   #91 --> 87
data_rmna_2 <- data_rmna_2[!is.na(data_rmna_2$sodium),]   #87 --> 76
data_rmna_2$gf <- 0
data_rmna_2$gf[which(data_rmna_2$z_wt.age < -2)] <- 1
data_rmna_2$gf[which(data_rmna_2$z_len.age < -2)] <- 1
levels(data_rmna_2$aclf_status.x) = c("ACLF" = 1, "non-ACLF" = 0)   
data_rmna_2$aclf_status.x <- factor(data_rmna_2$aclf_status.x, levels = c("0", "1"))
data_rmna_2 <- data_rmna_2[-43,]

# FINAL MODEL
# use data_rmna_2 as dataset since it includes both EKG and BA Data
data_rmna_2$aclf_status <- data_rmna_2$aclf_status.x
y <- as.numeric(data_rmna_2$aclf_status)- 1
x <- data.frame(data_rmna_2[,c(11:15,26,27,31,42,44:47,66:71,76)]) 

estimates <- rep(0.5,20)
estimates[c(11,13)] <- 1
model_b1 <- bic.glm(x,y,strict = FALSE, OR = 20, glm.family = "binomial",
                    factor.type = TRUE, prior.param = estimates)
#from this: select (1) inr (2) TB (3) listing_age (4) z_len.age (5) GF

form <- formula(aclf_status ~ (inr + total_bili + listing_age + z_len.age + gf)^2)
estimates <- rep(0.5,15)
estimates[c(1,2)] <- 1

model_b2 <- bic.glm(form, data = data_rmna_2, glm.family = "binomial",
                    OR = 20, strict = FALSE, prior.param = estimates)

#from this: select (1) inr (2) TB (3) TB:gf (4) listing_age (5) inr:listing_age

logmodel_b2 <- glm(aclf_status.x ~ total_bili + inr + total_bili:gf +
                     listing_age, family = "binomial",
                   data = data_rmna_2)
roc_lb2 <- roc(data_rmna_2$aclf_status.x, predict(logmodel_b2))
roc_peld_temp <- roc(data_rmna_2$aclf_status.x, data_rmna_2$listing_peld)

# remove influential points
data_rmna_2_new <- data_rmna_2[-c(13,44,47),]
logmodel_b2 <- glm(aclf_status.x ~ total_bili + inr + total_bili:gf +
                     listing_age, family = "binomial",
                   data = data_rmna_2_new)
roc_lb2 <- roc(data_rmna_2_new$aclf_status.x, predict(logmodel_b2))
roc_peld_temp <- roc(data_rmna_2_new$aclf_status.x, data_rmna_2_new$listing_peld)
roc_balf_temp <- roc(data_rmna_2_new$aclf_status.x, data_rmna_2_new$listing_balf)
coords(roc_lb2, "best", best.method = "closest.topleft")
coords(roc_peld_temp, "best", best.method = "closest.topleft")
coords(roc_balf_temp, "best", best.method = "youden")

plot(roc_lb2) #cex.axis = 1.25, cex.lab = 1.5)
plot(roc_peld_temp, add = TRUE, lty = 2)
plot(roc_balf_temp, add = TRUE, lty = 3)

roc.test(roc_lb2, roc_peld_temp)


data_rmna_2_new$propscore <- (0.11046*data_rmna_2_new$total_bili +
                              2.68606*data_rmna_2_new$inr -
                              0.03337*data_rmna_2_new$listing_age + 
                              0.23804*(data_rmna_2_new$total_bili*data_rmna_2_new$gf))
data_rmna_2_new$propscore2 <- 54 + data_rmna_2_new$propscore
rocblah <- roc(data_rmna_2_new$aclf_status, data_rmna_2_new$propscore2)

# validate on larger data

data_rmna$propscore <- (0.11046*data_rmna$total_bili +
                                2.68606*data_rmna$inr -
                                0.03337*data_rmna$listing_age + 
                                0.23804*(data_rmna$total_bili*data_rmna$gf))
data_rmna$propscore2 <- 54 + data_rmna$propscore
rocblah2 <- roc(data_rmna$aclf_status, data_rmna$propscore2)


#library(ResourceSelection)
#h1 <- hoslem.test(logmodel_b2$y, fitted(logmodel_b2), g = 10)  # 0.80
#library(rms)
#data_rmna_2_new$aclf_statusz <- as.n
#m1 <- lrm(y ~ total_bili + inr + total_bili:gf +
#            listing_age, data = data_rmna_2_new)


## Pt. B Try to incorporate Echo results
# data_rmna_2 is for BA + EKG combined dataset

data_rmna_3 <- ba_echo_data[!is.na(ba_echo_data$lv_midwall_sys_z),]  #98 --> 92
data_rmna_3 <- data_rmna_3[!is.na(data_rmna_3$platelets),]
levels(data_rmna_3$aclf_status) = c("ACLF" = 1, "non-ACLF" = 0)   
data_rmna_3$aclf_status <- factor(data_rmna_3$aclf_status, levels = c("0", "1"))


data_rmna_3 <- ba_echo_data[!is.na(ba_echo_data$platelets),]  #98 --> 92
data_rmna_3 <- data_rmna_3[!is.na(data_rmna_3$z_wt.len),]  #92 --> 84
data_rmna_3 <- data_rmna_3[!is.na(data_rmna_3$lv_midwall_sys_z),] #84 --> 54
data_rmna_3$gf <- 0
data_rmna_3$gf[which(data_rmna_3$z_wt.age < -2)] <- 1
data_rmna_3$gf[which(data_rmna_3$z_len.age < -2)] <- 1
levels(data_rmna_3$aclf_status) = c("ACLF" = 1, "non-ACLF" = 0)   
data_rmna_3$aclf_status <- factor(data_rmna_3$aclf_status, levels = c("0", "1"))

y <- as.numeric(data_rmna_3$aclf_status)- 1
x <- data.frame(data_rmna_3[,c(7,8,26,27,28,47:51,52,97)]) 

estimates <- rep(0.5,12)
estimates[c(3,5)] <- 1
model_b1 <- bic.glm(x,y,strict = FALSE, OR = 20, glm.family = "binomial",
                    factor.type = TRUE, prior.param = estimates)
#from this: select (1) inr (2) TB (3) listing_age (4) z_len.age (5) GF


noint_model <- glm(aclf_status ~ inr + total_bili + lv_midwall_sys_z + platelets,  
                   data = data_rmna_3, family = "binomial")
roc_temp_obj <- roc(data_rmna_3$aclf_status, predict(noint_model))  #ROC: 0.91
roc_peld_comp <- roc(data_rmna_3$aclf_status, data_rmna_3$listing_peld) #ROC: 0.86


noint_model <- glm(aclf_status ~ inr + total_bili + listing_age, 
                   data = data_rmna, family = "binomial")
roc_temp_obj <- roc(data_rmna$aclf_status, predict(noint_model))  #ROC: 0.91
print(roc_temp_obj)
roc_peld_comp <- roc(data_rmna$aclf_status, data_rmna$listing_peld) #ROC: 0.86
print(roc_peld_comp)

# remove influential points

data_rmna_noip <- data_rmna[-c(15,20,28,65,66),]
noint_model <- glm(aclf_status ~ inr + total_bili + z_len.age + z_wt.age, 
                   data = data_rmna_noip, family = "binomial")
roc_temp_obj <- roc(data_rmna_noip$aclf_status, predict(noint_model))  #ROC: 0.91
print(roc_temp_obj)
roc_peld_comp <- roc(data_rmna_noip$aclf_status, data_rmna_noip$listing_peld) #ROC: 0.86
print(roc_peld_comp)



form <- formula(aclf_status ~ (inr + total_bili + listing_age + z_len.age + gf)^2)
estimates <- rep(0.5,15)
estimates[c(1,2)] <- 1

model_b2 <- bic.glm(form, data = data_rmna_2, glm.family = "binomial",
                    OR = 20, strict = FALSE, prior.param = estimates)

#from this: select (1) inr (2) TB (3) TB:gf (4) listing_age (5) inr:listing_age

logmodel_b2 <- glm(aclf_status.x ~ total_bili + inr + total_bili:gf +
                     listing_age, family = "binomial",
                   data = data_rmna_2)
roc_lb2 <- roc(data_rmna_2$aclf_status.x, predict(logmodel_b2))
roc_peld_temp <- roc(data_rmna_2$aclf_status.x, data_rmna_2$listing_peld)




####### 4 - MODEL VALIDATION  #######

# 4-fold cross validation

fourfoldCV <- function(fml){
  
  rand_nums_aclf <- sample(1:17,17)
  rand_nums_non_aclf <- sample(1:68,68)
  
  aclf_fold1 <- aclf[rand_nums_aclf[1:5],]
  aclf_fold2 <- aclf[rand_nums_aclf[6:9],]
  aclf_fold3 <- aclf[rand_nums_aclf[10:13],]
  aclf_fold4 <- aclf[rand_nums_aclf[14:17],]
  
  non_aclf_fold1 <- non_aclf[rand_nums_non_aclf[1:17],]
  non_aclf_fold2 <- non_aclf[rand_nums_non_aclf[18:34],]
  non_aclf_fold3 <- non_aclf[rand_nums_non_aclf[35:51],]
  non_aclf_fold4 <- non_aclf[rand_nums_non_aclf[52:68],]
  
  fold1 <- rbind(aclf_fold1, non_aclf_fold1)
  fold2 <- rbind(aclf_fold2, non_aclf_fold2)
  fold3 <- rbind(aclf_fold3, non_aclf_fold3)
  fold4 <- rbind(aclf_fold4, non_aclf_fold4)
  train1 <- rbind(fold1, fold2, fold3)
  test1 <- fold4
  
  round1_bma_1 <- glm(fml, data = train1, family = "binomial")
  round1_bma_1_predictions <- predict(round1_bma_1, test1)
  round1_roc_1 <- roc(test1$aclf_status, round1_bma_1_predictions)   #AUC: 0.85
  round1_roc_peld <- roc(test1$aclf_status, test1$listing_peld)   #AUC: 0.77
  
  train2 <- rbind(fold1, fold2, fold4)
  test2 <- fold3
  
  round2_bma_1 <- glm(fml, data = train2, family = "binomial")
  round2_bma_1_predictions <- predict(round2_bma_1, test2)
  round2_roc_1 <- roc(test2$aclf_status, round2_bma_1_predictions)   #AUC: 0.78
  round2_roc_peld <- roc(test2$aclf_status, test2$listing_peld)   #AUC: 0.74
  
  train3 <- rbind(fold1, fold3, fold4)
  test3 <- fold2
  
  round3_bma_1 <- glm(fml,data = train3, family = "binomial")
  round3_bma_1_predictions <- predict(round3_bma_1, test3)
  round3_roc_1 <- roc(test3$aclf_status, round3_bma_1_predictions)   #AUC: 0.72
  round3_roc_peld <- roc(test3$aclf_status, test3$listing_peld)   #AUC: 0.75
  
  train4 <- rbind(fold2, fold3, fold4)
  test4 <- fold1
  
  round4_bma_1 <- glm(fml,data = train4, family = "binomial")
  round4_bma_1_predictions <- predict(round4_bma_1, test4)
  round4_roc_1 <- roc(test4$aclf_status, round4_bma_1_predictions)   #AUC: 0.72
  round4_roc_peld <- roc(test4$aclf_status, test4$listing_peld)   #AUC: 0.75
  
  peld_auc <- c(as.numeric(round1_roc_peld$auc),as.numeric(round2_roc_peld$auc),
                as.numeric(round3_roc_peld$auc),as.numeric(round4_roc_peld$auc))
  
  model_auc <- c(as.numeric(round1_roc_1$auc),as.numeric(round2_roc_1$auc),
                 as.numeric(round3_roc_1$auc),as.numeric(round4_roc_1$auc))
  
  print(peld_auc)
  print(mean(peld_auc))
  print(model_auc)
  print(mean(model_auc))
}

# Final model 1
fml <- formula(aclf_status ~ total_bili + inr + listing_age + total_bili:gf)
aclf <- subset(data_rmna_2_new, aclf_status == 1)
non_aclf <- subset(data_rmna_2_new, aclf_status == 0)
fourfoldCV_v2(fml)

fourfoldCV_v2 <- function(fml){
  
  rand_nums_aclf <- sample(1:12,12)
  rand_nums_non_aclf <- sample(1:60,60)
  
  aclf_fold1 <- aclf[rand_nums_aclf[1:3],]
  aclf_fold2 <- aclf[rand_nums_aclf[4:6],]
  aclf_fold3 <- aclf[rand_nums_aclf[7:9],]
  aclf_fold4 <- aclf[rand_nums_aclf[10:12],]
  
  non_aclf_fold1 <- non_aclf[rand_nums_non_aclf[1:15],]
  non_aclf_fold2 <- non_aclf[rand_nums_non_aclf[16:30],]
  non_aclf_fold3 <- non_aclf[rand_nums_non_aclf[31:45],]
  non_aclf_fold4 <- non_aclf[rand_nums_non_aclf[46:60],]
  
  fold1 <- rbind(aclf_fold1, non_aclf_fold1)
  fold2 <- rbind(aclf_fold2, non_aclf_fold2)
  fold3 <- rbind(aclf_fold3, non_aclf_fold3)
  fold4 <- rbind(aclf_fold4, non_aclf_fold4)
  train1 <- rbind(fold1, fold2, fold3)
  test1 <- fold4
  
  round1_bma_1 <- glm(fml, data = train1, family = "binomial")
  round1_bma_1_predictions <- predict(round1_bma_1, test1)
  round1_roc_1 <- roc(test1$aclf_status, round1_bma_1_predictions)   #AUC: 0.85
  round1_roc_peld <- roc(test1$aclf_status, test1$listing_peld)   #AUC: 0.77
  
  train2 <- rbind(fold1, fold2, fold4)
  test2 <- fold3
  
  round2_bma_1 <- glm(fml, data = train2, family = "binomial")
  round2_bma_1_predictions <- predict(round2_bma_1, test2)
  round2_roc_1 <- roc(test2$aclf_status, round2_bma_1_predictions)   #AUC: 0.78
  round2_roc_peld <- roc(test2$aclf_status, test2$listing_peld)   #AUC: 0.74
  
  train3 <- rbind(fold1, fold3, fold4)
  test3 <- fold2
  
  round3_bma_1 <- glm(fml,data = train3, family = "binomial")
  round3_bma_1_predictions <- predict(round3_bma_1, test3)
  round3_roc_1 <- roc(test3$aclf_status, round3_bma_1_predictions)   #AUC: 0.72
  round3_roc_peld <- roc(test3$aclf_status, test3$listing_peld)   #AUC: 0.75
  
  train4 <- rbind(fold2, fold3, fold4)
  test4 <- fold1
  
  round4_bma_1 <- glm(fml,data = train4, family = "binomial")
  round4_bma_1_predictions <- predict(round4_bma_1, test4)
  round4_roc_1 <- roc(test4$aclf_status, round4_bma_1_predictions)   #AUC: 0.72
  round4_roc_peld <- roc(test4$aclf_status, test4$listing_peld)   #AUC: 0.75
  
  peld_auc <- c(as.numeric(round1_roc_peld$auc),as.numeric(round2_roc_peld$auc),
                as.numeric(round3_roc_peld$auc),as.numeric(round4_roc_peld$auc))
  
  model_auc <- c(as.numeric(round1_roc_1$auc),as.numeric(round2_roc_1$auc),
                 as.numeric(round3_roc_1$auc),as.numeric(round4_roc_1$auc))
  
  print(peld_auc)
  print(mean(peld_auc))
  print(model_auc)
  print(mean(model_auc))
}

####### 5 - OTHER MISC. ANALYSES #######
# Kasai Analysis
kasai_data <- subset(data, kasai_status == 1)
nokasai_data <- subset(data, kasai_status == 0)
kasai_knowndate_data <- subset(kasai_data, is.na(kasai_data$kasai_age) == FALSE)
kasai_early <- subset(kasai_knowndate_data, kasai_age < 43)
kasai_late <- subset(kasai_knowndate_data, kasai_age > 42)
kasai_early_aclf <- subset(kasai_early, aclf_status == "ACLF")
kasai_early_nonaclf <- subset(kasai_early, aclf_status == "non-ACLF")
kasai_late_aclf <- subset(kasai_late, aclf_status == "ACLF")
kasai_late_nonaclf <- subset(kasai_late, aclf_status == "non-ACLF")
kasai_unknowndate <- subset(kasai_data , is.na(kasai_data$kasai_age) == TRUE)
kasai_unknowndate_aclf <- subset(kasai_unknowndate, aclf_status == "ACLF")
kasai_unknowndate_nonaclf <- subset(kasai_unknowndate, aclf_status == "non-ACLF")
nokasai_aclf <- subset(nokasai_data, aclf_status == "ACLF")
nokasai_nonaclf <- subset(nokasai_data, aclf_status == "non-ACLF")

# normality assessment among echo parameters
qqnorm(aclf$lv_midwall_sys_z, main = "Q-Q Plot for 
       LV Midwall Sys Z (ACLF)")
qqline(aclf$lv_midwall_sys_z)
shapiro.test(aclf$lv_midwall_sys_z)  # 0.26

qqnorm(non_aclf$lv_midwall_sys_z, main = "Q-Q Plot for
       LV Midwall Sys Z (Non-ACLF)")
qqline(non_aclf$lv_midwall_sys_z)
shapiro.test(non_aclf$lv_midwall_sys_z)



# APPENDIX 

# Run hypothesis tests on BA data 
analysis1('listing_age', aclf, non_aclf)
analysis1('treatment_qtydays', aclf, non_aclf)
analysis1('wait', aclf, non_aclf)
analysis1('alk_phos', aclf, non_aclf)
analysis1('creat', aclf, non_aclf)
analysis1('sodium', aclf, non_aclf)
analysis1('total_bili', aclf, non_aclf)
analysis1('albumin', aclf, non_aclf)
analysis1('inr', aclf, non_aclf)
analysis1('length', aclf, non_aclf)
analysis1('weight', aclf, non_aclf)
analysis1('bmi', aclf, non_aclf)
analysis1('z_wt.len', aclf, non_aclf)
analysis1('z_wt.age', aclf, non_aclf)
analysis1('z_len.age', aclf, non_aclf)
analysis1('z_bmi.age', aclf, non_aclf)
analysis1('ggt', aclf, non_aclf)
analysis1('platelets', aclf, non_aclf)
analysis1('conj_bili', aclf, non_aclf)
analysis1('listing_peld', aclf, non_aclf)
analysis1('listing_balf', aclf, non_aclf)
analysis1('adj_balf', aclf, non_aclf)
analysis1('kasai_age', aclf, non_aclf)
analysis1('wait', aclf, non_aclf)

# Run hypothesis tests on EKG data
analysis1("vent_rate", ekg_aclf, ekg_nonaclf)
analysis1("pr_int", ekg_aclf, ekg_nonaclf)
analysis1("qrs", ekg_aclf, ekg_nonaclf)
analysis1("qt", ekg_aclf, ekg_nonaclf)
analysis1("qtc", ekg_aclf, ekg_nonaclf)
analysis1("p.axis", ekg_aclf, ekg_nonaclf)
analysis1("r.axis", ekg_aclf, ekg_nonaclf)
analysis1("t.axis", ekg_aclf, ekg_nonaclf)

# Run hypothesis tests on echo data
analysis1("lvmi", aclf, non_aclf)
analysis1("mmlv_mass_height", aclf, non_aclf)
analysis1("mmlv_mass_2.7", aclf, non_aclf)
analysis1("rwt", aclf, non_aclf)
analysis1("e", aclf, non_aclf)
analysis1("a", aclf, non_aclf)
analysis1("s", aclf, non_aclf)
analysis1("mv_e", aclf, non_aclf)
analysis1("mv_a", aclf, non_aclf)
analysis1("mv_s", aclf, non_aclf)
analysis1("tv_e", aclf, non_aclf)
analysis1("tv_a", aclf, non_aclf)
analysis1("tv_s", aclf, non_aclf)
analysis1("bsa", aclf, non_aclf)
analysis1("mmlv_frac_short", aclf, non_aclf)
analysis1("mmlv_frac_short_z", aclf, non_aclf)
analysis1("mmlv_mass", aclf, non_aclf)
analysis1("mmlv_mass_z", aclf, non_aclf)
analysis1("endo_fs", aclf, non_aclf)
analysis1("endo_fs_z", aclf, non_aclf)
analysis1("midwall_fs", aclf, non_aclf)
analysis1("midwall_fs_z", aclf, non_aclf)
analysis1("lv_post", aclf, non_aclf)
analysis1("lv_sys", aclf, non_aclf)
analysis1("lv_sys_z", aclf, non_aclf)
analysis1("lv_dias", aclf, non_aclf)
analysis1("lv_dias_z", aclf, non_aclf)
analysis1("lv_midwall_sys", aclf, non_aclf)
analysis1("lv_midwall_sys_z", aclf, non_aclf)
analysis1("lv_midwall_dias", aclf, non_aclf)
analysis1("lv_midwall_dias_z", aclf, non_aclf)

# Stepwise regression method
model1 <- glm(aclf_status ~ total_bili + inr, data = data_rmna, 
              family = "binomial")
model2 <- glm(aclf_status ~ total_bili + inr + ggt, data = data_rmna,
              family = "binomial")
fullmodel <- glm(aclf_status ~  ggt  + albumin
                 + inr + total_bili + z_wt.len + z_wt.age + z_len.age + 
                   z_bmi.age + listing_age + wait, data = data_rmna,
                 family = "binomial")
back_model <- step(fullmodel)
model3 <- glm(aclf_status ~ total_bili + inr + ggt + z_wt.len + z_len.age +
                z_wt.age, data = data_rmna, family = "binomial")
back_model_2 <- step(back_model, ~.^2)
model5 <- glm(aclf_status ~ total_bili + inr + ggt + z_wt.len + z_len.age +
                inr:z_len.age + ggt:total_bili + ggt:z_wt.len, data = data_rmna,
              family = "binomial")


# 3 Fold cross validation function
threefoldCV <- function(fold1,fold2,fold3,fml){
  train1 <- rbind(fold1, fold2)
  test1 <- fold3
  
  round1_bma_1 <- glm(fml, data = train1, family = "binomial")
  round1_bma_1_predictions <- predict(round1_bma_1, test1)
  round1_roc_1 <- roc(test1$aclf_status, round1_bma_1_predictions)   #AUC: 0.85
  round1_roc_peld <- roc(test1$aclf_status, test1$listing_peld)   #AUC: 0.77
  
  train2 <- rbind(fold1, fold3)
  test2 <- fold2
  
  round2_bma_1 <- glm(fml, data = train2, family = "binomial")
  round2_bma_1_predictions <- predict(round2_bma_1, test2)
  round2_roc_1 <- roc(test2$aclf_status, round2_bma_1_predictions)   #AUC: 0.78
  round2_roc_peld <- roc(test2$aclf_status, test2$listing_peld)   #AUC: 0.74
  
  train3 <- rbind(fold2, fold3)
  test3 <- fold1
  
  round3_bma_1 <- glm(fml,data = train3, family = "binomial")
  round3_bma_1_predictions <- predict(round3_bma_1, test3)
  round3_roc_1 <- roc(test3$aclf_status, round3_bma_1_predictions)   #AUC: 0.72
  round3_roc_peld <- roc(test3$aclf_status, test3$listing_peld)   #AUC: 0.75
  
  print(round1_roc_1)
  print(round2_roc_1)
  print(round3_roc_1)
  print(round1_roc_peld)
  print(round2_roc_peld)
  print(round3_roc_peld)
  
}

# Tables for paper
plot(ba_ekg_data$listing_peld,ba_ekg_data$qtc, ylim = c(420,480))
plot(ba_echo_data$listing_peld,ba_echo_data$lv_midwall_sys_z)
lm <- lm(lv_midwall_sys_z ~ listing_peld, data = ba_echo_data)
abline(0.39128,0.01399)


# nonsurviving ACLF
#add age at kasai


# Table 2 Work
aclf <- subset(data, aclf_status == "ACLF")
non_aclf <- subset(data, aclf_status == "non-ACLF")

table1 <- cbind(c(10,10),c(31,54))
chisq.test(table1)
table2 <- cbind(c(3,17),c(29,56))
chisq.test(table2)
table3 <- cbind(c(7,13),c(25,60))
chisq.test(table3)
table<- cbind(c(10,3,7),c(31,29,25))
chisq.test(table)

gf_cohort <- subset(data, gf == 1)
middle_cohort <- subset(data, gf == 0 & nogf == 0)
nogf_cohort <- subset(data, nogf == 1)

# Table 3 Work

summary(as.factor(aclf$kasai_status))
summary(as.factor(non_aclf$kasai_status))
chisq.test(cbind(c(6,14),c(24,70)))
chisq.test(cbind(c(5,6),c(24,70)))

summary(aclf$kasai_age > 42)
summary(non_aclf$kasai_age > 42)
chisq.test(cbind(c(11,2),c(42,8)))
fisher.test(cbind(c(11,2),c(42,8)))

summary(aclf$kasai_age < 42 & aclf$kasai_age > 28)
summary(non_aclf$kasai_age < 42 & non_aclf$kasai_age > 28)
fisher.test(cbind(c(1,12),c(4,46)))

# Table 4

ekg_fem <- subset(ba_ekg_data, sex == "F")
ekg_fem_01month <- subset(ekg_fem, ecg_age < 30)
ekg_fem_13month <- subset(ekg_fem, ecg_age > 30 & ecg_age < 90)   #3
ekg_fem_36month <- subset(ekg_fem, ecg_age >= 90 & ecg_age < 180)   #22
mean(ekg_fem_36month$qtc, na.rm = TRUE)
sd(ekg_fem_36month$qtc, na.rm = TRUE)
t.test(ekg_fem_36month$qtc, mu = 418)
wilcox.test(ekg_fem_36month$qtc, mu = 418)
ekg_fem_612month <- subset(ekg_fem, ecg_age > 180 & ecg_age < 360)  #23
mean(ekg_fem_612month$qtc)
sd(ekg_fem_612month$qtc)
t.test(ekg_fem_612month$qtc, mu = 414)
wilcox.test(ekg_fem_612month$qtc, mu = 414)
ekg_fem_13year <- subset(ekg_fem, ecg_age > 360 & ecg_age < 1080)  #10
ekg_fem_35year <- subset(ekg_fem, ecg_age > 1080 & ecg_age < 1800)  # 2
ekg_fem_58year <- subset(ekg_fem, ecg_age > 1800 & ecg_age < 2880) #1
ekg_fem_812year <- subset(ekg_fem, ecg_age > 2880 & ecg_age < 4320) # 1

ekg_male <- subset(ba_ekg_data, sex == "M")
ekg_male_01month <- subset(ekg_male, ecg_age < 30)
ekg_male_13month <- subset(ekg_male, ecg_age> 30 & ecg_age < 90)   # 1
ekg_male_36month <- subset(ekg_male, ecg_age >= 90 & ecg_age < 180)   # 11
mean(ekg_male_36month$qtc)
sd(ekg_male_36month$qtc)
t.test(ekg_male_36month$qtc, mu = 422)
wilcox.test(ekg_male_36month$qtc, mu = 422)
ekg_male_612month <- subset(ekg_male, ecg_age >= 180 & ecg_age < 360)  #12
mean(ekg_male_612month$qtc)
sd(ekg_male_612month$qtc)
t.test(ekg_male_612month$qtc, mu = 411)
wilcox.test(ekg_male_612month$qtc, mu = 411)
ekg_male_13year <- subset(ekg_male, ecg_age > 360 & ecg_age < 1080)  #7
ekg_male_35year <- subset(ekg_male, ecg_age > 1080 & ecg_age < 1800)  #3
ekg_male_58year <- subset(ekg_male, ecg_age > 1800 & ecg_age < 2880) #1
ekg_male_812year <- subset(ekg_male, ecg_age > 2880 & ecg_age < 4320) #0
ekg_male_1216year <- subset(ekg_male, ecg_age > 4320)  #1


# QTC delta analysis
aclf <- subset(ba_ekg_data, aclf_status.x == "ACLF")
non_aclf <- subset(ba_ekg_data, aclf_status.x == "non-ACLF")
analysis1('delta', aclf, non_aclf)

aclf_new <- subset(aclf, sex == "M" & ecg_age > 360 & ecg_age < 1080)
non_aclf_new <- subset(non_aclf, sex == "M" & ecg_age > 360 & ecg_age < 1080)


# Scatterplots
plot(ba_ekg_data$listing_peld, ba_ekg_data$qtc, xlab = "Listing PELD",
     ylab = "QTC")
abline(lm(qtc~listing_peld, data = ba_ekg_data))

plot(ba_echo_data$listing_peld,ba_echo_data$lv_midwall_sys_z, 
     xlab = "Listing PELD", ylab = "Left ventricular midwall systolic z-score")
abline(lm(lv_midwall_sys_z ~ listing_peld, data = ba_echo_data))
summary(lm(lv_midwall_sys_z ~ listing_peld, data = ba_echo_data))

plot(ba_echo_data$listing_peld, ba_echo_data$mmlv_mass_z,
     xlab = "Listing PELD", ylab = "mm LV mass z-score")
abline(lm(mmlv_mass_z~listing_peld, data = ba_echo_data))
summary(lm(mmlv_mass_z ~ listing_peld, data= ba_echo_data))

plot(ba_echo_data$listing_peld, ba_echo_data$e, xlab = "Listing PELD",
     ylab = "E'")
abline(lm(e~listing_peld, data = ba_echo_data))
summary(lm(e ~ listing_peld, data= ba_echo_data))


# Work for May 2015

# use a new dataset

ba_echo_data_mod <- read.csv("combo_echo_data_v2.csv")
# summary(ba_echo_data_mod$mmlv_mass_2.7_delta)    #35 NA's --> 63 values
aclf <- subset(ba_echo_data_mod, aclf_status == "ACLF")
non_aclf <- subset(ba_echo_data_mod, aclf_status == "non-ACLF")

analysis1('mmlv_mass_2.7_delta', aclf, non_aclf)


aclf_young <- subset(aclf, listing_age < 365)
non_aclf_young <- subset(non_aclf, listing_age < 365)
analysis1('e', aclf_young, non_aclf_young)
analysis1('a', aclf_young, non_aclf_young)
analysis1('s', aclf_young, non_aclf_young)
analysis1('mv_e', aclf_young, non_aclf_young)
analysis1('mv_a', aclf_young, non_aclf_young)
analysis1('mv_s', aclf_young, non_aclf_young)
analysis1('tv_e', aclf_young, non_aclf_young)
analysis1('tv_a', aclf_young, non_aclf_young)
analysis1('tv_s', aclf_young, non_aclf_young)


summary(ba_echo_data_mod$mmlv_mass_2.7_delta, na.rm = TRUE)
sd(ba_echo_data_mod$mmlv_mass_2.7_delta, na.rm = TRUE)
t.test(ba_echo_data_mod$mmlv_mass_2.7_delta)
wilcox.test(ba_echo_data_mod$mmlv_mass_2.7_delta)


summary(ba_echo_data_mod$mmlv_frac_short_z)
sd(ba_echo_data_mod$mmlv_frac_short_z, na.rm = TRUE)
t.test(ba_echo_data_mod$mmlv_frac_short_z)
wilcox.test(ba_echo_data_mod$mmlv_frac_short_z)

summary(ba_echo_data_mod$mmlv_mass_z)
sd(ba_echo_data_mod$mmlv_mass_z, na.rm = TRUE)
t.test(ba_echo_data_mod$mmlv_mass_z)
wilcox.test(ba_echo_data_mod$mmlv_mass_z)

summary(ba_echo_data_mod$endo_fs_z)
sd(ba_echo_data_mod$endo_fs_z, na.rm = TRUE)
t.test(ba_echo_data_mod$endo_fs_z)
wilcox.test(ba_echo_data_mod$endo_fs_z)

summary(ba_echo_data_mod$midwall_fs_z)
sd(ba_echo_data_mod$midwall_fs_z, na.rm = TRUE)
t.test(ba_echo_data_mod$midwall_fs_z)
wilcox.test(ba_echo_data_mod$midwall_fs_z)

summary(ba_echo_data_mod$lv_sys_z)
sd(ba_echo_data_mod$lv_sys_z, na.rm = TRUE)
t.test(ba_echo_data_mod$lv_sys_z)
wilcox.test(ba_echo_data_mod$lv_sys_z)

summary(ba_echo_data_mod$lv_dias_z)
sd(ba_echo_data_mod$lv_dias_z, na.rm = TRUE)
t.test(ba_echo_data_mod$lv_dias_z)
wilcox.test(ba_echo_data_mod$lv_dias_z)

summary(ba_echo_data_mod$lv_midwall_sys_z)
sd(ba_echo_data_mod$lv_midwall_sys_z, na.rm = TRUE)
t.test(ba_echo_data_mod$lv_midwall_sys_z)
wilcox.test(ba_echo_data_mod$lv_midwall_sys_z)

summary(ba_echo_data$lv_midwall_dias_z)
sd(ba_echo_data$lv_midwall_dias_z, na.rm = TRUE)
t.test(ba_echo_data$lv_midwall_dias_z)
wilcox.test(ba_echo_data$lv_midwall_dias_z)



