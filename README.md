# Condom-Scale-Data
---
title: "Condom_GRPA_merged_data"
output: html_document
---

```{r}
#gpra baseline data, renaming ID 
grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98,-88,-99))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB)
dim(grpa_short)
names(grpa_short)


#condom scale baseline data, renaming ID
condom=read.csv("Condom Scale - Baseline.csv", header=TRUE, na.strings=c(98,-88,-99))
dim(condom)
names(condom)[1]= "PARTID"
names(condom)


#merging both datasets 
condom_grpa=merge(condom,grpa_short, by="PARTID", all.x=TRUE)
dim(condom_grpa)
names(condom_grpa)
head(condom_grpa)
```

Redcap data
```{r}
#Loading redcap data
grpa_redcap = read.csv("CCPE_RedCap_GRPA_Condom_Data.csv", header = TRUE, na.strings=c(98,-88,-99))
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

colnames(binded_data)[colnames(binded_data)=="YOB"] <- "Age"
binded_data$Age=2019-binded_data$Age
names(binded_data)
head(binded_data)

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender
```
Next steps
Make YOB of their age
Check for errors: can use describe, summary
```{r}
binded_data_error = subset(binded_data, EffectivePreventSTD > 7) 
binded_data_error
#participant ID #1767 invalid response for EffectivePrevent STD --> answer of 33 
#find out where their data is living/where the hard copy of the data is 
binded_data_error1=subset(binded_data, ConvenientToUse > 7) 
binded_data_error1
#participant ID #1671 invalid response for ConvenientToUse --> answer of 22 
summary(binded_data)

```




