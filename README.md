# Condom-Scale-Data
---
title: "Condom_GRPA_merged_data"
output: html_document
---

Cleaning data
```{r}
#gpra baseline data, renaming ID 
grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB)
dim(grpa_short)
names(grpa_short)


#condom scale baseline data, renaming ID
condom=read.csv("Condom Scale - Baseline.csv", header=TRUE, na.strings=c(98))
dim(condom)
names(condom)[1]= "PARTID"
names(condom)


#merging both datasets 
install.packages("prettyR")
library(prettyR)
describe.factor(grpa_short$SEX_PR)
condom_grpa=merge(condom,grpa_short, by="PARTID", all.x=TRUE)
dim(condom_grpa)
names(condom_grpa)
head(condom_grpa)
```

Redcap data
```{r}
#Loading redcap data
grpa_redcap = read.csv("CCPE_RedCap_GRPA_Condom_Data.csv", header = TRUE, na.strings=c(98))
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
Age
names(binded_data)
head(binded_data)

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender

library(prettyR)
describe.factor
```
Next steps
Make YOB of their age
Check for errors: can use describe, summary
```{r}
library(psych)
```



