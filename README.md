# Condom-Scale-Data
---
title: "Condom_GRPA_merged_data"
output: html_document
---
Libraring packages
```{r}
library(prettyR)
library(psych)
library(car)
```



Cleaning data
```{r}
#gpra baseline data, renaming ID 
#grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98, 99, 97))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB)
dim(grpa_short)
names(grpa_short)


#condom scale baseline data, renaming ID
#condom=read.csv("Condom Scale - Baseline.csv", header=TRUE, na.strings=c(98, 99, 97))
dim(condom)
names(condom)[1]= "PARTID"
names(condom)


#merging both datasets 
describe.factor(grpa_short$SEX_PR)
condom_grpa=merge(condom,grpa_short, by="PARTID", all.x=TRUE)
dim(condom_grpa)
names(condom_grpa)
head(condom_grpa)
```

Redcap data
```{r}
#Loading redcap data
#grpa_redcap = read.csv("CCPE_RedCap_GRPA_Condom_Data.csv", header = TRUE, na.strings=c(98, 99, 97))
names(grpa_redcap)[1:211]<- toupper(names(grpa_redcap)[1:211])

#creating new data frame with appropriate variables 
grpa_redcap_short=data.frame(PARTID=grpa_redcap$PARTID, R_BLACK_N=grpa_redcap$R_BLACK_N, R_WHITE_N=grpa_redcap$R_WHITE_N, GENDER=grpa_redcap$GENDER, SEX_PR=grpa_redcap$SEX_PR, YOB=grpa_redcap$YOB, PreventPregnancy=grpa_redcap$RAPREVENTPREGNANCY, EffectivePreventSTD=grpa_redcap$RAEFFECTIVEPREVENTSTD, EffectivePreventHIV=grpa_redcap$RAEFFECTIVEPREVENTHIV, Comfortable=grpa_redcap$RACOMFORTABLE, ConvenientToUse=grpa_redcap$RACONVENIENTTOUSE, SexualPleasure=grpa_redcap$RASEXUALPLEASURE, Obtain=grpa_redcap$RAOBTAIN, Freinds=grpa_redcap$RAFRIENDS, SexualPartner=grpa_redcap$RASEXUALPARTNER, ExcitingDull=grpa_redcap$RAEXCITINGDULL, Embarrassing=grpa_redcap$RAEMBARRASSING, DiscussWithPartner=grpa_redcap$RADISCUSSWITHPARTNER, EasyHardToUse=grpa_redcap$RAEASYHARDTOUSE, NeatMessy=grpa_redcap$RANEATMESSY)
dim(grpa_redcap_short)
head(grpa_redcap_short)
names(grpa_redcap_short)

#binding grpa_redcap shortened dataset with the merged condom_grpa dataset 
binded_data=rbind(grpa_redcap_short, condom_grpa)
head(binded_data)
dim(binded_data)
names(binded_data)

write.csv(binded_data, file="binded_data.csv")

colnames(binded_data)[colnames(binded_data)=="YOB"] <- "Age"
binded_data$Age=2019-binded_data$Age
names(binded_data)
head(binded_data)

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender

```
Getting larger data set with possible predictors
```{r}
#gpra baseline data, renaming ID 
#grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98,-88,-99, 99, 97))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB, EDLEVEL_N=grpa$EDLEVEL_N, WRGSEX_UNP_A=grpa$WRGSEX_UNP_A, RSKANYSEX_UNP=grpa$RSKANYSEX_UNP, RSKSEX_ALCDRG=grpa$RSKSEX_ALCDRG, CNTRL_REFUSEMOOD=grpa$CNTRL_REFUSEMOOD, CNTRL_WAITCNDM=grpa$CNTRL_WAITCNDM, CNTRL_TREAT=grpa$CNTRL_TREAT, CNTRL_SEXPRAC=grpa$CNTRL_SEXPRAC, CNTRL_ASKCNDM=grpa$CNTRL_ASKCNDM, CNTRL_REFUSECNDM=grpa$CNTRL_REFUSECNDM, KNOW_HIV=grpa$KNOW_HIV, MENTLH30D=grpa$MENTLH30D, SEX_ANY30D=grpa$SEX_ANY30D , LASTSEX_UNP=grpa$LASTSEX_UNP, SEX_MNY_3MOS=grpa$SEX_MNY_3MOS, WHENSEXALCDRG=grpa$WHENSEXALCDRG, ANYABUSE_3M=grpa$ANYABUSE_3M, SEXUNWANT_12M=grpa$SEXUNWANT_12M, LIVE_N=grpa$LIVE_N, RSKCIG=grpa$RSKCIG, RSKMJ=grpa$RSKMJ, RSKALC=grpa$RSKALC, GET_MEDHLP=grpa$GET_MEDHLP, HIV_RESULTS_N=grpa$HIV_RESULTS_N, HOMETYPE_N=grpa$HOMETYPE_N, HC_HAVE_N=grpa$HC_HAVE_N, PEERBINGE_A=grpa$PEERBINGE_A, WRGBINGE_A=grpa$WRGBINGE_A, KNOW_SA=grpa$KNOW_SA, BINGE530D=grpa$BINGE530D, ALC30D=grpa$ALC30D, MJ30D=grpa$MJ30D, CUT_ALC=grpa$CUT_ALC, ANNOY_ALC=grpa$ANNOY_ALC, GUILT_ALC=grpa$GUILT_ALC, EMO_AFT=grpa$EMO_AFT, HIV_DRGS_N=grpa$HIV_DRGS_N, HIV_CURE_N=grpa$HIV_CURE_N )
names(grpa_short)

#condom scale baseline data, renaming ID
#condom=read.csv("Condom Scale - Baseline.csv", header=TRUE, na.strings=c(98,-88,-99, 98, 97))
names(condom)[1]= "PARTID"

#merging both datasets 
condom_grpa=merge(condom,grpa_short, by="PARTID", all.x=TRUE)
```

Redcap data
```{r}
#Loading redcap data
#grpa_redcap = read.csv("CCPE_RedCap_GRPA_Condom_Data.csv", header = TRUE, na.strings=c(98,-88,-99))
names(grpa_redcap)[1:211]<- toupper(names(grpa_redcap)[1:211])

#creating new data frame with appropriate variables 
grpa_redcap_short=data.frame(PARTID=grpa_redcap$PARTID, R_BLACK_N=grpa_redcap$R_BLACK_N, R_WHITE_N=grpa_redcap$R_WHITE_N, GENDER=grpa_redcap$GENDER, SEX_PR=grpa_redcap$SEX_PR, YOB=grpa_redcap$YOB, EDLEVEL_N=grpa_redcap$EDLEVEL_N, WRGSEX_UNP_A=grpa_redcap$WRGSEX_UNP_A, RSKANYSEX_UNP=grpa_redcap$RSKANYSEX_UNP, RSKSEX_ALCDRG=grpa_redcap$RSKSEX_ALCDRG, CNTRL_REFUSEMOOD=grpa_redcap$CNTRL_REFUSEMOOD, CNTRL_WAITCNDM=grpa_redcap$CNTRL_WAITCNDM, CNTRL_TREAT=grpa_redcap$CNTRL_TREAT, CNTRL_SEXPRAC=grpa_redcap$CNTRL_SEXPRAC, CNTRL_ASKCNDM=grpa_redcap$CNTRL_ASKCNDM, CNTRL_REFUSECNDM=grpa_redcap$CNTRL_REFUSECNDM, KNOW_HIV=grpa_redcap$KNOW_HIV, MENTLH30D=grpa_redcap$MENTLH30D, SEX_ANY30D=grpa_redcap$SEX_ANY30D , LASTSEX_UNP=grpa_redcap$LASTSEX_UNP, SEX_MNY_3MOS=grpa_redcap$SEX_MNY_3MOS, WHENSEXALCDRG=grpa_redcap$WHENSEXALCDRG, ANYABUSE_3M=grpa_redcap$ANYABUSE_3M, SEXUNWANT_12M=grpa_redcap$SEXUNWANT_12M, LIVE_N=grpa_redcap$LIVE_N, RSKCIG=grpa_redcap$RSKCIG, RSKMJ=grpa_redcap$RSKMJ, RSKALC=grpa_redcap$RSKALC, GET_MEDHLP=grpa_redcap$GET_MEDHLP, HIV_RESULTS_N=grpa_redcap$HIV_RESULTS_N, HOMETYPE_N=grpa_redcap$HOMETYPE_N, HC_HAVE_N=grpa_redcap$HC_HAVE_N, PEERBINGE_A=grpa_redcap$PEERBINGE_A, WRGBINGE_A=grpa_redcap$WRGBINGE_A, KNOW_SA=grpa_redcap$KNOW_SA, ALC30D=grpa_redcap$ALC30D, BINGE530D=grpa_redcap$BINGE530D, MJ30D=grpa_redcap$MJ30D, CUT_ALC=grpa_redcap$CUT_ALC, ANNOY_ALC=grpa_redcap$ANNOY_ALC, GUILT_ALC=grpa_redcap$GUILT_ALC, EMO_AFT=grpa_redcap$EMO_AFT, HIV_DRGS_N=grpa_redcap$HIV_DRGS_N, HIV_CURE_N=grpa_redcap$HIV_CURE_N, PreventPregnancy=grpa_redcap$RAPREVENTPREGNANCY, EffectivePreventSTD=grpa_redcap$RAEFFECTIVEPREVENTSTD, EffectivePreventHIV=grpa_redcap$RAEFFECTIVEPREVENTHIV, Comfortable=grpa_redcap$RACOMFORTABLE, ConvenientToUse=grpa_redcap$RACONVENIENTTOUSE, SexualPleasure=grpa_redcap$RASEXUALPLEASURE, Obtain=grpa_redcap$RAOBTAIN, Freinds=grpa_redcap$RAFRIENDS, SexualPartner=grpa_redcap$RASEXUALPARTNER, ExcitingDull=grpa_redcap$RAEXCITINGDULL, Embarrassing=grpa_redcap$RAEMBARRASSING, DiscussWithPartner=grpa_redcap$RADISCUSSWITHPARTNER, EasyHardToUse=grpa_redcap$RAEASYHARDTOUSE, NeatMessy=grpa_redcap$RANEATMESSY)

names(grpa_redcap_short)

#binding grpa_redcap shortened dataset with the merged condom_grpa dataset 
binded_data=rbind(grpa_redcap_short, condom_grpa)

colnames(binded_data)[colnames(binded_data)=="YOB"] <- "Age"
binded_data$Age=2019-binded_data$Age

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender
```

```{r}
#collapse education variable 
#1= minority/what most interested in--> reference group is always 0 ; 1=greater than highschool education, 0=high school education or less 


### This doesn't work, because the NAs turn into zeros
Education=recode(binded_data$EDLEVEL_N, "c('1', '2,', '3') ='0'; 98='98'; else='1'")

## Try this
binded_data$EDLEVEL_N_bin = ifelse(binded_data$EDLEVEL_N <= 3, 0, 1)


#collapse sexunwant variable
#0=never have had unwanted sex in last 12 months 1=ever had unwanted sex in the last 12 months (2,3,4)
Sexunwant=recode(binded_data$SEXUNWANT_12M, "c('2','3','4')='1'; 98='98'; else='0'")
binded_data$Sexunwant = ifelse(binded_data$SEXUNWANT_12M == 1, 0, 0)
describe.factor(binded_data$Sexunwant)
## This variable is not useful, because there is no variation
#collapse anyabuse variable 
#0= never have been abused emotionally, physically, or sexually within past 3 months with someone with whom youve been in a relationship in, 1=ever been abused emotionally, physically, or sexually within the past 3 months 
Anyabuse=recode(binded_data$ANYABUSE_3M, "c('2','3','4', '5')='1'; 98='98'; else='0'")
describe.factor(binded_data$ANYABUSE_3M)
binded_data$Anyabuse = ifelse(binded_data$ANYABUSE_3M == 1, 0, 1)
describe.factor(binded_data$Anyabuse)
```

```{r}
#total score for perception of risk variables 
Risk_Total= data.frame(RSKCIG=binded_data$RSKCIG, RSKMJ=binded_data$RSKMJ, RSKALC=binded_data$RSKALC)

### This is summing the everything.  We want a sum across
Perception_Risk_Total=sum(Risk_Total, na.rm=TRUE)

## Not right, because it is turning into 0's into 
Perception_Risk_Total= rowSums(Risk_Total, na.rm=TRUE)
setwd("C:/Users/Matthew.Hanauer/Desktop")
write.csv(Perception_Risk_Total, "Perception_Risk_Total.csv", row.names = FALSE)
Perception_Risk_Total = read.csv("Perception_Risk_Total.csv", header = TRUE, na.strings= c(0))
describe.factor(Perception_Risk_Total)


Perception_Risk_Total
range(Perception_Risk_Total)
describe.factor(Perception_Risk_Total)
#Total score for relationship with (main) partner confidence 
Confidence_Total=data.frame(CNTRL_REFUSEMOOD=binded_data$CNTRL_REFUSEMOOD, CNTRL_WAITCNDM=binded_data$CNTRL_WAITCNDM, CNTRL_TREAT=binded_data$CNTRL_TREAT, CNTRL_SEXPRAC=binded_data$CNTRL_SEXPRAC, CNTRL_ASKCNDM=binded_data$CNTRL_ASKCNDM, CNTRL_REFUSECNDM=binded_data$CNTRL_REFUSECNDM)
Partner_Confidence_Total=sum(Confidence_Total, na.rm=TRUE)
Partner_Confidence_Total

#binding total score variables and collapsed variables with main binded dataset 
binded_data2=cbind(binded_data, Education, Sexunwant, Perception_Risk_Total, Partner_Confidence_Total, Anyabuse)
names(binded_data2)
  
```

```{r}
#invalid data cleaning
binded_data_error = subset(binded_data, EffectivePreventSTD > 7) 
binded_data_error
#participant ID #1767 invalid response for EffectivePrevent STD --> answer of 33 
#find out where their data is living/where the hard copy of the data is 
binded_data_error1=subset(binded_data, ConvenientToUse > 7) 
binded_data_error1
#participant ID #1671 invalid response for ConvenientToUse --> answer of 22 
```

```{r}
#deletes any row that has any missing data 
binded_data_complete= na.omit(binded_data)
dim(binded_data_complete)
1-(dim(binded_data_complete)[1]/dim(binded_data)[1]) 
```

```{r}
install.packages("prettyR")

#Apply function; percentages for all MISSING variables 
apply(binded_data, 2, function(col)sum(is.na(col))/length(col))
#percentages for all missing data for each participant
apply(binded_data, 1, function(col)sum(is.na(col))/length(col))

#count of complete/incomplete cases 
sum(complete.cases(binded_data))
sum(!complete.cases(binded_data))

#Number of missing per row 
rowSums(is.na(binded_data))

#Number of missing per column/variable
colSums(is.na(binded_data))

#considers the pattern in which the users have missed a particular question
missing_data=lapply(binded_data[,-1], function(x){binded_data[!complete.cases(x), 'PARTID']})
missing_data

#all rows in a dataset with atleast one NA
unique (unlist(lapply(binded_data, function (x) which (is.na(x)))))
head(binded_data)
apply(binded_data2, 2, function(x){describe.factor(x)})
```

Create factors for variables
```{r}
#counts and percentages for all demographic variables 
binded_data_cat = binded_data[,2:5]
binded_data_cat = apply(binded_data_cat, 2, function(x){describe.factor(x)})
binded_data_cat
```

```{r}
#mean and standard deviation for all variables  
library(prettyR)
describe(binded_data)
```

```{r}
#creating new data set for condom scale items without demographic variables 
Condom_Scale_TotalScore=data.frame(PreventPregnancy=binded_data$PreventPregnancy, EffectivePreventSTD=binded_data$EffectivePreventSTD, EffectivePreventHIV=binded_data$EffectivePreventHIV, Comfortable=binded_data$Comfortable, ConvenientToUse=binded_data$ConvenientToUse, SexualPleasure=binded_data$SexualPleasure, Obtain=binded_data$Obtain, Freinds=binded_data$Freinds, SexualPartner=binded_data$SexualPartner, ExcitingDull=binded_data$ExcitingDull, Embarrassing=binded_data$Embarrassing,  DiscussWithPartner=binded_data$DiscussWithPartner,  EasyHardToUse=binded_data$EasyHardToUse, NeatMessy=binded_data$NeatMessy)

#creating three total variables and attaching them to the dataset 
Perceived_effectiveness= data.frame(PreventPregnancy=binded_data$PreventPregnancy, EffectivePreventSTD=binded_data$EffectivePreventSTD, EffectivePreventHIV=binded_data$EffectivePreventHIV)

Affective=data.frame(Comfortable=binded_data$Comfortable, SexualPleasure=binded_data$SexualPleasure, Freinds=binded_data$Freinds, SexualPartner=binded_data$SexualPartner, ExcitingDull=binded_data$ExcitingDull, NeatMessy=binded_data$NeatMessy,ConvenientToUse=binded_data$ConvenientToUse)

Manageability=data.frame(EasyHardToUse=binded_data$EasyHardToUse, Embarrassing=binded_data$Embarrassing, DiscussWithPartner=binded_data$DiscussWithPartner)


#Total score for condom variables 
Total_Score_CondomScale=apply(Condom_Scale_TotalScore, 1, sum)
mean(Total_Score_CondomScale, na.rm=TRUE)

PE_Total_Score=sum(Perceived_effectiveness, na.rm=TRUE)
PE_Total_Score

Affective_Total_Score=sum(Affective, na.rm=TRUE)
Affective_Total_Score

Manageability_Total_Score=sum(Manageability, na.rm=TRUE)
Manageability_Total_Score


#Putting Total_score variable back into the original dataset 
binded_data1 =  cbind(binded_data2, Total_Score_CondomScale,PE_Total_Score,Affective_Total_Score,Manageability_Total_Score)
names(binded_data1)

#age ranges 
range(binded_data1$Age, na.rm = TRUE)
```
Reiability

```{r}
head(binded_data[45:58])

omegaItems = na.omit((binded_data[45:58]))
dim(omegaItems)
head(binded_data)
omegaResults = omega(omegaItems)
summary(omegaResults)


```
Splitting Data

Create a binded data items for the efa to make it easier.  We are keeping the missing data for EFA and CFA, so the N should not change throughout this process.
```{r}
binded_data_items = binded_data[45:58]
binded_data_items$ID = 1:dim(binded_data)[1]

inTrain = createDataPartition(binded_data_items$ID, p = .50, list = FALSE)
efaData = binded_data_items[inTrain,]
cfaData = binded_data_items[-inTrain,]
efaData$ID = NULL
```
EFA
```{r}
efa1 = fa(r = efaData, nfactors = 1, fm = "gls", cor = "poly")

efa2 = fa(r = efaData, nfactors = 2, fm = "gls", cor = "poly")
fa.diagram(efa2)

efa3 = fa(r = efaData, nfactors = 3, fm = "gls", cor = "poly")
fa.diagram(efa3)
summary(efa3)


efa4 = fa(r = efaData, nfactors = 4, fm = "gls", cor = "poly")
fa.diagram(efa4)
summary(efa4)

anova(efa1, efa2)
anova(efa1, efa3)
anova(efa2, efa3)
anova(efa3, efa4)

```
EFA with just the three foru items
```{r}
efaData_pe = data.frame(efaData[,1:3],EasyHardToUse = efaData$EasyHardToUse)
efa1_pe = fa(r = efaData_pe, nfactors = 1, fm = "gls", cor = "poly")
fa.diagram(efa1_pe)
summary(efa1_pe)

vss(efaData_pe, n = 3, rotate = "oblimin", fm = "mle", cor = "poly")

efaData_pe_complete = na.omit(efaData_pe)

paran(efaData_pe_complete, centile = 95, iterations = 1000, graph = TRUE, cfa = TRUE)

```


VSS
```{r}
vss(efaData, n = 3, rotate = "oblimin", fm = "mle", cor = "poly")
vss(efaData, n = 4, rotate = "oblimin", fm = "mle", cor = "poly")
```
paran
```{r}
efaData_complete = na.omit(efaData)

paran(efaData_complete, centile = 95, iterations = 1000, graph = TRUE, cfa = TRUE)

```




CFA
```{r}
library(lavaan)


## Original model presented
model1 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull + NeatMessy
          M =~ ConvenientToUse + Obtain  + Embarrassing + DiscussWithPartner + EasyHardToUse'

fit1 = cfa(model1, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit1, fit.measures = TRUE, standardized = TRUE)

head(binded_data)

## Try a one construct model
model2 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV + Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull + NeatMessy + ConvenientToUse + EasyHardToUse + Embarrassing + DiscussWithPartner + Obtain'

fit2 = cfa(model2, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit2, fit.measures = TRUE, standardized = TRUE)

### CFA 

## Try exact model from EFA
model3 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ ConvenientToUse +Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull +     NeatMessy
          M =~ Obtain  + Embarrassing + DiscussWithPartner + EasyHardToUse'

fit3 = cfa(model3, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit3, fit.measures = TRUE, standardized = TRUE)

modificationindices(fit3)


### Try making changes based on mod indicies
model4 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable + Freinds + SexualPartner + ExcitingDull   
          M =~ Obtain  + Embarrassing + DiscussWithPartner + EasyHardToUse + NeatMessy + SexualPleasure + ConvenientToUse'

fit4 = cfa(model4, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit4, fit.measures = TRUE, standardized = TRUE)

modificationindices(fit3)


### Ok now try double loading, but what does that mean???? So I guess it goes into the total score twice for the different constructs???
model5 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ ConvenientToUse +Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull +     NeatMessy
          M =~ Obtain  + Embarrassing + DiscussWithPartner + EasyHardToUse + NeatMessy + SexualPleasure + ConvenientToUse'

fit5 = cfa(model5, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit5, fit.measures = TRUE, standardized = TRUE)

modificationindices(fit3)


#### Try four factor model based on EFA
model6 = '
          F1_m =~ DiscussWithPartner + Embarrassing + EasyHardToUse + NeatMessy + Obtain
          F2 =~ ExcitingDull + SexualPleasure + SexualPartner + Freinds
          PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          F4 =~ ConvenientToUse + Comfortable'

fit6 = cfa(model6, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit6, fit.measures = TRUE, standardized = TRUE)



### Try a model with just the effective questions.  Those seem to be the ones that cling together the most.
model7 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV  + EasyHardToUse'

fit7 = cfa(model7, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit7, fit.measures = TRUE, standardized = TRUE)

### Try with the second factor of easy to use
model8 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV  + EasyHardToUse
          Ease =~ ConvenientToUse + Comfortable'

fit8 = cfa(model8, estimator  = "MLR", missing = "ML", data = binded_data_items)
summary(fit8, fit.measures = TRUE, standardized = TRUE)

```
Final Models

I don't think it is happening.  Even the good model is not that good and I cannot justify the easy hard to use inclusion.  
```{r}
model5 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull +ConvenientToUse
M =~ EasyHardToUse + Embarrassing + DiscussWithPartner'



fit5 = cfa(model5, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit5, fit.measures = TRUE, standardized = TRUE)

### Try a model with just the effective questions.  Those seem to be the ones that cling together the most.
model7 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV  + EasyHardToUse'

fit7 = cfa(model7, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit7, fit.measures = TRUE, standardized = TRUE)

```
########################################################
Trying psychomterics with confidence variables from GPRA
########################################################
```{r}
head(Confidence_Total)

Confidence_Total$ID = 1:dim(Confidence_Total)[1]

inTrain = createDataPartition(Confidence_Total$ID, p = .50, list = FALSE)
efaData = Confidence_Total[inTrain,]
cfaData = Confidence_Total[-inTrain,]
efaData$ID = NULL


efa1 = fa(r = efaData, nfactors = 1, fm = "gls", cor = "poly")
fa.diagram(efa1)
summary(efa1)

efa2= fa(r = efaData, nfactors = 2, fm = "gls", cor = "poly")
fa.diagram(efa2)
summary(efa2)

anova(efa1, efa2)

efa3= fa(r = efaData, nfactors = 3, fm = "gls", cor = "poly")
fa.diagram(efa3)
summary(efa3)

anova(efa2, efa3)

vss(efaData, n = 3, rotate = "oblimin", fm = "mle", cor = "poly")

efaData_complete = na.omit(efaData)

paran(efaData_complete, centile = 95, iterations = 1000, graph = TRUE, cfa = TRUE)


#### Try CFA with one-factor
model_con1 = 'PE =~ CNTRL_REFUSEMOOD + CNTRL_WAITCNDM + CNTRL_TREAT  + CNTRL_SEXPRAC + CNTRL_ASKCNDM + CNTRL_REFUSECNDM'

fit_con1 = cfa(model_con1, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit_con1, fit.measures = TRUE, standardized = TRUE)


model_con2 = 
            'F1 =~  CNTRL_WAITCNDM  + CNTRL_ASKCNDM + CNTRL_REFUSECNDM
             F2 =~ CNTRL_REFUSEMOOD + CNTRL_TREAT + CNTRL_SEXPRAC'

fit_con2 = cfa(model_con2, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit_con2, fit.measures = TRUE, standardized = TRUE)
```



