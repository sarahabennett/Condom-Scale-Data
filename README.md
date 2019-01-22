# Condom-Scale-Data
---
title: "Condom_GRPA_merged_data"
output: html_document
---

Cleaning data
```{r}
#gpra baseline data, renaming ID 
grpa=read.csv("CCPE_GRPA_Baseline_Condom.csv", header=TRUE, na.strings=c(98,-88,-99))
grpa_short=data.frame(PARTID=grpa$PARTID, R_BLACK_N=grpa$R_BLACK_N, R_WHITE_N=grpa$R_WHITE_N, GENDER=grpa$GENDER, SEX_PR=grpa$SEX_PR, YOB=grpa$YOB, SPEAK_ENG=grpa$SPEAK_ENG, EDLEVEL_N=grpa$EDLEVEL_N, WRGSEX_UNP_A=grpa$WRGSEX_UNP_A, RSKANYSEX_UNP=grpa$RSKANYSEX_UNP, RSKSEX_ALCDRG=grpa$RSKSEX_ALCDRG, CNTRL_REFUSEMOOD=grpa$CNTRL_REFUSEMOOD, CNTRL_WAITCNDM=grpa$CNTRL_WAITCNDM, CNTRL_TREAT=grpa$CNTRL_TREAT, CNTRL_SEXPRAC=grpa$CNTRL_SEXPRAC, CNTRL_ASKCNDM=grpa$CNTRL_ASKCNDM, CNTRL_REFUSECNDM=grpa$CNTRL_REFUSECNDM, KNOW_HIV=grpa$KNOW_HIV, MENTLH30D=grpa$MENTLH30D, SEX_HAD=grpa$SEX_HAD, SEX_ANY30D=grpa$SEX_ANY30D, LASTSEX_UNP=grpa$LASTSEX_UNP, SEX_MNY_3MOS=grpa$SEX_MNY_3MOS, WHENSEX4STFF_UNP=grpa$WHENSEX4STFF_UNP, WHENSEXHIVSTD_UNP=grpa$WHENSEXHIVSTD_UNP, WHENSEXDRUG_UNP=grpa$WHENSEXDRUG_UNP, WHENSEXALCDRG=grpa$WHENSEXALCDRG, ANYABUSE_3M=grpa$ANYABUSE_3M, SEXUNWANT_12M=grpa$SEXUNWANT_12M, MSTATUS_N=grpa$MSTATUS_N, LIVE_N=grpa$LIVE_N )

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
grpa_redcap_short=data.frame(PARTID=grpa_redcap$PARTID, R_BLACK_N=grpa_redcap$R_BLACK_N, R_WHITE_N=grpa_redcap$R_WHITE_N, GENDER=grpa_redcap$GENDER, SEX_PR=grpa_redcap$SEX_PR, YOB=grpa_redcap$YOB, SPEAK_ENG=grpa_redcap$SPEAK_ENG, EDLEVEL_N=grpa_redcap$EDLEVEL_N, WRGSEX_UNP_A=grpa_redcap$WRGSEX_UNP_A, RSKANYSEX_UNP=grpa_redcap$RSKANYSEX_UNP, RSKSEX_ALCDRG=grpa_redcap$RSKSEX_ALCDRG, CNTRL_REFUSEMOOD=grpa_redcap$CNTRL_REFUSEMOOD, CNTRL_WAITCNDM=grpa_redcap$CNTRL_WAITCNDM, CNTRL_TREAT=grpa_redcap$CNTRL_TREAT, CNTRL_SEXPRAC=grpa_redcap$CNTRL_SEXPRAC, CNTRL_ASKCNDM=grpa_redcap$CNTRL_ASKCNDM, CNTRL_REFUSECNDM=grpa_redcap$CNTRL_REFUSECNDM, KNOW_HIV=grpa_redcap$KNOW_HIV, MENTLH30D=grpa_redcap$MENTLH30D, SEX_HAD=grpa_redcap$SEX_HAD, SEX_ANY30D=grpa_redcap$SEX_ANY30D, LASTSEX_UNP=grpa_redcap$LASTSEX_UNP, SEX_MNY_3MOS=grpa_redcap$SEX_MNY_3MOS, WHENSEX4STFF_UNP=grpa_redcap$WHENSEX4STFF_UNP, WHENSEXHIVSTD_UNP=grpa_redcap$WHENSEXHIVSTD_UNP, WHENSEXDRUG_UNP=grpa_redcap$WHENSEXDRUG_UNP, WHENSEXALCDRG=grpa_redcap$WHENSEXALCDRG, ANYABUSE_3M=grpa_redcap$ANYABUSE_3M, SEXUNWANT_12M=grpa_redcap$SEXUNWANT_12M, MSTATUS_N=grpa_redcap$MSTATUS_N, LIVE_N=grpa_redcap$LIVE_N,  PreventPregnancy=grpa_redcap$RAPREVENTPREGNANCY, EffectivePreventSTD=grpa_redcap$RAEFFECTIVEPREVENTSTD, EffectivePreventHIV=grpa_redcap$RAEFFECTIVEPREVENTHIV, Comfortable=grpa_redcap$RACOMFORTABLE, ConvenientToUse=grpa_redcap$RACONVENIENTTOUSE, SexualPleasure=grpa_redcap$RASEXUALPLEASURE, Obtain=grpa_redcap$RAOBTAIN, Freinds=grpa_redcap$RAFRIENDS, SexualPartner=grpa_redcap$RASEXUALPARTNER, ExcitingDull=grpa_redcap$RAEXCITINGDULL, Embarrassing=grpa_redcap$RAEMBARRASSING, DiscussWithPartner=grpa_redcap$RADISCUSSWITHPARTNER, EasyHardToUse=grpa_redcap$RAEASYHARDTOUSE, NeatMessy=grpa_redcap$RANEATMESSY)

#binding grpa_redcap shortened dataset with the merged condom_grpa dataset 
binded_data=rbind(grpa_redcap_short, condom_grpa)

colnames(binded_data)[colnames(binded_data)=="YOB"] <- "Age"
binded_data$Age=2019-binded_data$Age

#gender: there are 2 four's, meaning 2 people responded saying they are unsure of their gender
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
#Apply function; percentages for all MISSING variables 
apply(binded_data, 2, function(col)sum(is.na(col))/length(col))
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

#Total score for condom variables 
Total_score=apply(Condom_Scale_TotalScore, 1, sum)
head(Total_score)
mean(Total_score, na.rm=TRUE)

#creating three total variables and attaching them to the dataset 
Perceived_effectiveness= data.frame(PreventPregnancy=binded_data$PreventPregnancy, EffectivePreventSTD=binded_data$EffectivePreventSTD, EffectivePreventHIV=binded_data$EffectivePreventHIV)

Affective=data.frame(Comfortable=binded_data$Comfortable, SexualPleasure=binded_data$SexualPleasure, Freinds=binded_data$Freinds, SexualPartner=binded_data$SexualPartner, ExcitingDull=binded_data$ExcitingDull, NeatMessy=binded_data$NeatMessy)

Manageability=data.frame(ConvenientToUse=binded_data$ConvenientToUse, EasyHardToUse=binded_data$EasyHardToUse, Embarrassing=binded_data$Embarrassing, DiscussWithPartner=binded_data$DiscussWithPartner, EasyHardToUse=binded_data$EasyHardToUse)

PE_Total_Score=sum(Perceived_effectiveness, na.rm=TRUE)
PE_Total_Score


Affective_Total_Score=sum(Affective_Total_Score, na.rm=TRUE)
Affective_Total_Score

Manageability_Total_Score=sum(Manageability, na.rm=TRUE)
Manageability_Total_Score


#Putting Total_score variable back into the original dataset 
binded_data1 =  cbind(binded_data, Total_score,PE_Total_Score,Affective_Total_Score,Manageability_Total_Score)
names(binded_data1)

#age ranges 
range(binded_data1$Age, na.rm = TRUE)
```

```{r}
#creating three total variables and attaching them to the dataset 
Perceived_effectiveness= data.frame(PreventPregnancy=binded_data$PreventPregnancy, EffectivePreventSTD=binded_data$EffectivePreventSTD, EffectivePreventHIV=binded_data$EffectivePreventHIV)

Affective=data.frame(Comfortable=binded_data$Comfortable, SexualPleasure=binded_data$SexualPleasure, Freinds=binded_data$Freinds, SexualPartner=binded_data$SexualPartner, ExcitingDull=binded_data$ExcitingDull, NeatMessy=binded_data$NeatMessy)

Manageability=data.frame(ConvenientToUse=binded_data$ConvenientToUse, EasyHardToUse=binded_data$EasyHardToUse, Embarrassing=binded_data$Embarrassing, DiscussWithPartner=binded_data$DiscussWithPartner, EasyHardToUse=binded_data$EasyHardToUse)

PE_Total_Score=sum(Perceived_effectiveness, na.rm=TRUE)
PE_Total_Score


Affective_Total_Score=sum(Affective_Total_Score, na.rm=TRUE)
Affective_Total_Score

Manageability_Total_Score=sum(Manageability, na.rm=TRUE)
Manageability_Total_Score

```


Reiability 
```{r}
omegaItems = na.omit((binded_data[,7:19]))
dim(omegaItems)
head(binded_data)
omegaResults = omega(omegaItems)
summary(omegaResults)
```
Splitting Data
```{r}
inTrain = createDataPartition(omegaItems$PreventPregnancy, p = .50, list = FALSE)
efaData = omegaItems[inTrain,]
cfaData = omegaItems[-inTrain,]

```
EFA
```{r}
efa1 = fa(r = efaData, nfactors = 1, fm = "gls", cor = "poly")

efa2 = fa(r = efaData, nfactors = 2, fm = "gls", cor = "poly")
fa.diagram(efa2)

efa3 = fa(r = efaData, nfactors = 3, fm = "gls", cor = "poly")
fa.diagram(efa3)
summary(efa3)

anova(efa1, efa2)
anova(efa1, efa3)
anova(efa2, efa3)
```
VSS
```{r}
vss(efaData, n = 3, rotate = "oblimin", fm = "mle", cor = "poly")
```
paran
```{r}
paran(efaData, centile = 95, iterations = 1000, graph = TRUE, cfa = TRUE)

```




CFA
```{r}
library(lavaan)
model1 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull + NeatMessy
          M =~ ConvenientToUse + EasyHardToUse + Embarrassing + DiscussWithPartner + Obtain'

fit1 = cfa(model1, estimator  = "MLR", missing = "ML", data = binded_data)
summary(fit1, fit.measures = TRUE, standardized = TRUE)

head(binded_data)


model2 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV + Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull + NeatMessy + ConvenientToUse + EasyHardToUse + Embarrassing + DiscussWithPartner + Obtain'

fit2 = cfa(model2, estimator  = "MLR", missing = "ML", data = binded_data)
summary(fit2, fit.measures = TRUE, standardized = TRUE)

### CFA 

model3 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull  +ConvenientToUse
          M =~ EasyHardToUse + Embarrassing + DiscussWithPartner + Obtain'

fit3 = cfa(model3, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit3, fit.measures = TRUE, standardized = TRUE)

head(binded_data)



model4 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull 
          M =~ EasyHardToUse + Embarrassing + DiscussWithPartner + Obtain + ConvenientToUse'

fit4 = cfa(model4, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit4, fit.measures = TRUE, standardized = TRUE)



model5 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull  +ConvenientToUse
          M =~ EasyHardToUse + Embarrassing + DiscussWithPartner'

fit5 = cfa(model5, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit5, fit.measures = TRUE, standardized = TRUE)


model6 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull
          M =~ EasyHardToUse + Embarrassing + DiscussWithPartner + ConvenientToUse'

fit6 = cfa(model6, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit6, fit.measures = TRUE, standardized = TRUE)


model7 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
          A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull
          M =~ EasyHardToUse + Embarrassing + DiscussWithPartner + ConvenientToUse
          PE ~~ M                                                                        
          '

fit7 = cfa(model7, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit7, fit.measures = TRUE, standardized = TRUE)


```
Final Model
```{r}
model5 = 'PE =~ PreventPregnancy + EffectivePreventSTD + EffectivePreventHIV
A =~ Comfortable+ SexualPleasure + Freinds + SexualPartner + ExcitingDull +ConvenientToUse
M =~ EasyHardToUse + Embarrassing + DiscussWithPartner'



fit5 = cfa(model5, estimator  = "MLR", missing = "ML", data = cfaData)
summary(fit5, fit.measures = TRUE, standardized = TRUE)
```
