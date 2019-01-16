# Condom-Scale-Data
---
title: "Condom_GRPA_merged_data"
output: html_document
---

Cleaning data
Cleaning data
```{r}
#gpra baseline data, renaming ID 
grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98,-88,-99))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB)

#condom scale baseline data, renaming ID
condom=read.csv("Condom Scale - Baseline.csv", header=TRUE, na.strings=c(98,-88,-99))
names(condom)[1]= "PARTID"

#merging both datasets 
condom_grpa=merge(condom,grpa_short, by="PARTID", all.x=TRUE)
```

Redcap data
```{r}
#Loading redcap data
grpa_redcap = read.csv("CCPE_RedCap_GRPA_Condom_Data.csv", header = TRUE, na.strings=c(98,-88,-99))
names(grpa_redcap)[1:211]<- toupper(names(grpa_redcap)[1:211])

#creating new data frame with appropriate variables 
grpa_redcap_short=data.frame(PARTID=grpa_redcap$PARTID, R_BLACK_N=grpa_redcap$R_BLACK_N, R_WHITE_N=grpa_redcap$R_WHITE_N, GENDER=grpa_redcap$GENDER, SEX_PR=grpa_redcap$SEX_PR, YOB=grpa_redcap$YOB, PreventPregnancy=grpa_redcap$RAPREVENTPREGNANCY, EffectivePreventSTD=grpa_redcap$RAEFFECTIVEPREVENTSTD, EffectivePreventHIV=grpa_redcap$RAEFFECTIVEPREVENTHIV, Comfortable=grpa_redcap$RACOMFORTABLE,b ConvenientToUse=grpa_redcap$RACONVENIENTTOUSE, SexualPleasure=grpa_redcap$RASEXUALPLEASURE, Obtain=grpa_redcap$RAOBTAIN, Freinds=grpa_redcap$RAFRIENDS, SexualPartner=grpa_redcap$RASEXUALPARTNER, ExcitingDull=grpa_redcap$RAEXCITINGDULL, Embarrassing=grpa_redcap$RAEMBARRASSING, DiscussWithPartner=grpa_redcap$RADISCUSSWITHPARTNER, EasyHardToUse=grpa_redcap$RAEASYHARDTOUSE, NeatMessy=grpa_redcap$RANEATMESSY)

#binding grpa_redcap shortened dataset with the merged condom_grpa dataset 
binded_data=rbind(grpa_redcap_short, condom_grpa)

colnames(binded_data)[colnames(binded_data)=="YOB"] <- "Age"
binded_data$Age=2019-binded_data$Age

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender
```

```{r}
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
#Apply function; percentages for all MISSING variables 
apply(binded_data, 2, function(col)sum(is.na(col))/length(col))
```

Create factors for variables
```{r}
#counts and percentages for all variables 
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

#Total score for condom variables 
Total_score=apply(Condom_Scale_TotalScore, 1, sum)
Total_score

#Putting Total_score variable back into the original dataset 
binded_data =  cbind(binded_data, Total_score)
names(binded_data)

#age ranges 
range(binded_data$Age, na.rm = TRUE)
```
