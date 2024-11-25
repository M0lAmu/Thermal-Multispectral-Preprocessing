setwd("C:/Users/thato/OneDrive/Desktop/MSc/DATA/Data Analysis")
install.packages("psych")
install.packages("tidyverse")
install.packages("readxl")
install.packages("reshape2")
install.packages("modeldata")
install.packages("Bonferroni")
library("psych")
library("tidyverse")
library("readxl")
library("reshape2")
library("modeldata")
library("Bonferroni")

read_excel("Temporary_VIs.xlsx")
Data <- read_excel("Temporary_VIs.xlsx")
View(Data)
summary(Data)

Data$REP <- as.factor(Data$REP)
Data$Treatment <- as.factor(Data$Treatment)
Data$Class <- as.factor(Data$Class)
class(Data$Treatment)
summary(Data)



#Get group summary statistics

describeBy(Data, Data$Class)
describeBy(Data, Data$Treatment)


#filter by Phenological stage

filter(Data, Class == "NSZ")
data_NSZ <- filter(Data, Class == "NSZ")
View(data_NSZ)

filter(Data, Class == "NTF")
data_NTF <- filter(Data, Class == "NTF")
View(data_NTF)



#Get grouped summary statistics by phenological stage, and Treatment at each phenological stage

describeBy(data_NSZ, data_NSZ$Treatment)

describeBy(data_NTF, data_NTF$Treatment)


##filter by treatment for Nut Sizing

filter(data_NSZ, Treatment == "1")
data_NSZ_1 <- filter(data_NSZ, Treatment == "1")
View(data_NSZ_1)

filter(data_NSZ, Treatment == "2")
data_NSZ_2 <- filter(data_NSZ, Treatment == "2")
View(data_NSZ_2)

filter(data_NSZ, Treatment == "3")
data_NSZ_3 <- filter(data_NSZ, Treatment == "3")
View(data_NSZ_3)

filter(data_NSZ, Treatment == "4")
data_NSZ_4 <- filter(data_NSZ, Treatment == "4")
View(data_NSZ_4)

filter(data_NSZ, Treatment == "5")
data_NSZ_5 <- filter(data_NSZ, Treatment == "5")
View(data_NSZ_5)


#Visualisation for Control and Treatment 3 during Nut sizing: ExcessGreen
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$ExG_mean), length(data_NSZ_3$ExG_mean)))
ExG <- c(data_NSZ_1$ExG_mean, data_NSZ_3$ExG_mean)

df <- data.frame(TreatmentGroups_NSZ, ExG)
boxplot(ExG ~ TreatmentGroups_NSZ, data = df)


#Visualisation for Control and Treatment 3 during Nut sizing: maskedExcessGreen
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mExG_mean), length(data_NSZ_3$mExG_mean)))
mExG <- c(data_NSZ_1$mExG_mean, data_NSZ_3$mExG_mean)

df2 <- data.frame(TreatmentGroups_NSZ, mExG)
boxplot(mExG ~ TreatmentGroups_NSZ, data = df2)


#Visualisation for Control and Treatment 3 during Nut sizing: NDVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$NDVI_mean), length(data_NSZ_3$NDVI_mean)))
NDVI <- c(data_NSZ_1$NDVI_mean, data_NSZ_3$NDVI_mean)

df3 <- data.frame(TreatmentGroups_NSZ, NDVI)
boxplot(NDVI ~ TreatmentGroups_NSZ, data = df3)


#Visualisation for Control and Treatment 3 during Nut sizing: mNDVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mNDVI_mean), length(data_NSZ_3$mNDVI_mean)))
mNDVI <- c(data_NSZ_1$mNDVI_mean, data_NSZ_3$mNDVI_mean)

df4 <- data.frame(TreatmentGroups_NSZ, mNDVI)
boxplot(mNDVI ~ TreatmentGroups_NSZ, data = df4)


#Visualisation for Control and Treatment 3 during Nut sizing: GNDVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$GNDVI_mean), length(data_NSZ_3$GNDVI_mean)))
GNDVI <- c(data_NSZ_1$GNDVI_mean, data_NSZ_3$GNDVI_mean)

df5 <- data.frame(TreatmentGroups_NSZ, GNDVI)
boxplot(GNDVI ~ TreatmentGroups_NSZ, data = df5)


#Visualisation for Control and Treatment 3 during Nut sizing: mGNDVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mGNDVI_mean), length(data_NSZ_3$mGNDVI_mean)))
mGNDVI <- c(data_NSZ_1$mGNDVI_mean, data_NSZ_3$mGNDVI_mean)

df6 <- data.frame(TreatmentGroups_NSZ, mGNDVI)
boxplot(mGNDVI ~ TreatmentGroups_NSZ, data = df6)


#Visualisation for Control and Treatment 3 during Nut sizing: GRVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$GRVI_mean), length(data_NSZ_3$GRVI_mean)))
GRVI <- c(data_NSZ_1$GRVI_mean, data_NSZ_3$GRVI_mean)

df7 <- data.frame(TreatmentGroups_NSZ, GRVI)
boxplot(GRVI ~ TreatmentGroups_NSZ, data = df7)


#Visualisation for Control and Treatment 3 during Nut sizing: mGRVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mGRVI_mean), length(data_NSZ_3$mGRVI_mean)))
mGRVI <- c(data_NSZ_1$mGRVI_mean, data_NSZ_3$mGRVI_mean)

df8 <- data.frame(TreatmentGroups_NSZ, mGRVI)
boxplot(mGRVI ~ TreatmentGroups_NSZ, data = df8)


#Visualisation for Control and Treatment 3 during Nut sizing: NDRE
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$NDRE_mean), length(data_NSZ_3$NDRE_mean)))
NDRE <- c(data_NSZ_1$NDRE_mean, data_NSZ_3$NDRE_mean)

df9 <- data.frame(TreatmentGroups_NSZ, NDRE)
boxplot(NDRE ~ TreatmentGroups_NSZ, data = df9)


#Visualisation for Control and Treatment 3 during Nut sizing: mNDRE
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mNDRE_mean), length(data_NSZ_3$mNDRE_mean)))
mNDRE <- c(data_NSZ_1$mNDRE_mean, data_NSZ_3$mNDRE_mean)

df10 <- data.frame(TreatmentGroups_NSZ, mNDRE)
boxplot(mNDRE ~ TreatmentGroups_NSZ, data = df10)



#Visualisation for Control and Treatment 3 during Nut sizing: REGI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$REGI_mean), length(data_NSZ_3$REGI_mean)))
REGI <- c(data_NSZ_1$REGI_mean, data_NSZ_3$REGI_mean)

df11 <- data.frame(TreatmentGroups_NSZ, REGI)
boxplot(REGI ~ TreatmentGroups_NSZ, data = df11)


#Visualisation for Control and Treatment 3 during Nut sizing: mREGI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                c(length(data_NSZ_1$mREGI_mean), length(data_NSZ_3$mREGI_mean)))
mREGI <- c(data_NSZ_1$mREGI_mean, data_NSZ_3$mREGI_mean)

df12 <- data.frame(TreatmentGroups_NSZ, mREGI)
boxplot(mREGI ~ TreatmentGroups_NSZ, data = df12)



#Visualisation for Control and Treatment 3 during Nut sizing: mARVI
TreatmentGroups_NSZ <- rep(c("Control", "Treatment 3"),
                  c(length(data_NSZ_1$mARVI_mean), length(data_NSZ_3$mARVI_mean)))
mARVI <- c(data_NSZ_1$mARVI_mean, data_NSZ_3$mARVI_mean)

df13 <- data.frame(TreatmentGroups_NSZ, mARVI)
boxplot(mARVI ~ TreatmentGroups_NSZ, data = df13)







#filter by treatment for Nut filling

filter(data_NTF, Treatment == "1")
data_NTF_1 <- filter(data_NTF, Treatment == "1")
View(data_NTF_1)

filter(data_NTF, Treatment == "2")
data_NTF_2 <- filter(data_NTF, Treatment == "2")
View(data_NTF_2)

filter(data_NTF, Treatment == "3")
data_NTF_3 <- filter(data_NTF, Treatment == "3")
View(data_NTF_3)

filter(data_NTF, Treatment == "4")
data_NTF_4 <- filter(data_NTF, Treatment == "4")
View(data_NTF_4)

filter(data_NTF, Treatment == "5")
data_NTF_5 <- filter(data_NTF, Treatment == "5")
View(data_NTF_5)


###VISUALISATIONS EXPLICITLY FOR 04 MARCH 2024
#Visualisation for Control and Treatment 4 during Nut Filling: ExG
TreatmentGroups <- rep(c("a", "b"),
                  c(length(data_NTF_1$ExG_mean), length(data_NTF_4$ExG_mean)))
VI <- c(data_NTF_1$ExG_mean, data_NTF_4$ExG_mean)

df14 <- data.frame(TreatmentGroups, VI)
boxplot(VI ~ TreatmentGroups, data = df14)


#Visualisation for Control and Treatment 4 during Nut Filling: mExG
TreatmentGroups2 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mExG_mean), length(data_NTF_4$mExG_mean)))
VI2 <- c(data_NTF_1$mExG_mean, data_NTF_4$mExG_mean)

df15 <- data.frame(TreatmentGroups2, VI2)
boxplot(VI2 ~ TreatmentGroups2, data = df15)


#Visualisation for Control and Treatment 4 during Nut Filling: NDVI
TreatmentGroups4 <- rep(c("a", "b"),
                       c(length(data_NTF_1$NDVI_mean), length(data_NTF_4$NDVI_mean)))
VI4 <- c(data_NTF_1$NDVI_mean, data_NTF_4$NDVI_mean)

df16 <- data.frame(TreatmentGroups4, VI4)
boxplot(VI4 ~ TreatmentGroups4, data = df16)


#Visualisation for Control and Treatment 4 during Nut Filling: mNDVI
TreatmentGroups5 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mNDVI_mean), length(data_NTF_4$mNDVI_mean)))
VI5 <- c(data_NTF_1$mNDVI_mean, data_NTF_4$mNDVI_mean)

df17 <- data.frame(TreatmentGroups5, VI5)
boxplot(VI5 ~ TreatmentGroups5, data = df17)


#Visualisation for Control and Treatment 4 during Nut Filling: GNDVI
TreatmentGroups6 <- rep(c("a", "b"),
                       c(length(data_NTF_1$GNDVI_mean), length(data_NTF_4$GNDVI_mean)))
VI6 <- c(data_NTF_1$GNDVI_mean, data_NTF_4$GNDVI_mean)

df18 <- data.frame(TreatmentGroups6, VI6)
boxplot(VI6 ~ TreatmentGroups6, data = df18)


#Visualisation for Control and Treatment 4 during Nut Filling: mGNDVI
TreatmentGroups7 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mGNDVI_mean), length(data_NTF_4$mGNDVI_mean)))
VI7 <- c(data_NTF_1$mGNDVI_mean, data_NTF_4$mGNDVI_mean)

df19 <- data.frame(TreatmentGroups7, VI7)
boxplot(VI7 ~ TreatmentGroups7, data = df19)


#Visualisation for Control and Treatment 4 during Nut Filling: GRVI
TreatmentGroups8 <- rep(c("a", "b"),
                       c(length(data_NTF_1$GRVI_mean), length(data_NTF_4$GRVI_mean)))
VI8 <- c(data_NTF_1$GRVI_mean, data_NTF_4$GRVI_mean)

df20 <- data.frame(TreatmentGroups8, VI8)
boxplot(VI8 ~ TreatmentGroups8, data = df20)


#Visualisation for Control and Treatment 4 during Nut Filling: mGRVI
TreatmentGroups9 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mGRVI_mean), length(data_NTF_4$mGRVI_mean)))
VI9 <- c(data_NTF_1$mGRVI_mean, data_NTF_4$mGRVI_mean)

df21 <- data.frame(TreatmentGroups9, VI9)
boxplot(VI9 ~ TreatmentGroups9, data = df21)


#Visualisation for Control and Treatment 4 during Nut Filling: NDRE
TreatmentGroups10 <- rep(c("a", "b"),
                       c(length(data_NTF_1$NDRE_mean), length(data_NTF_4$NDRE_mean)))
VI10 <- c(data_NTF_1$NDRE_mean, data_NTF_4$NDRE_mean)

df22 <- data.frame(TreatmentGroups10, VI10)
boxplot(VI10 ~ TreatmentGroups10, data = df22)


#Visualisation for Control and Treatment 4 during Nut Filling: mNDRE
TreatmentGroups11 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mNDRE_mean), length(data_NTF_4$mNDRE_mean)))
VI11 <- c(data_NTF_1$mNDRE_mean, data_NTF_4$mNDRE_mean)

df23 <- data.frame(TreatmentGroups11, VI11)
boxplot(VI11 ~ TreatmentGroups11, data = df23)


#Visualisation for Control and Treatment 4 during Nut Filling: REGI
TreatmentGroups12 <- rep(c("a", "b"),
                       c(length(data_NTF_1$REGI_mean), length(data_NTF_4$REGI_mean)))
VI12 <- c(data_NTF_1$REGI_mean, data_NTF_4$REGI_mean)

df24 <- data.frame(TreatmentGroups12, VI12)
boxplot(VI12 ~ TreatmentGroups12, data = df24)


#Visualisation for Control and Treatment 4 during Nut Filling: mREGI
TreatmentGroups13 <- rep(c("a", "b"),
                       c(length(data_NTF_1$mREGI_mean), length(data_NTF_4$mREGI_mean)))
VI13 <- c(data_NTF_1$mREGI_mean, data_NTF_4$mREGI_mean)

df25 <- data.frame(TreatmentGroups13, VI13)
boxplot(VI13 ~ TreatmentGroups13, data = df25)







##VISUALISATION FOR ALL TREATMENTS ON 04 MARCH 2024

TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                       c(length(data_NTF_1$ExG_mean), length(data_NTF_2$ExG_mean), length(data_NTF_3$ExG_mean),length(data_NTF_4$ExG_mean), length(data_NTF_5$ExG_mean)))
ExG <- c(data_NTF_1$ExG_mean, data_NTF_2$ExG_mean, data_NTF_3$ExG_mean, data_NTF_4$ExG_mean, data_NTF_5$ExG_mean)

df26 <- data.frame(TreatmentGroups_04032024, ExG)
boxplot(ExG ~ TreatmentGroups_04032024, data = df26)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mExG_mean), length(data_NTF_2$mExG_mean), length(data_NTF_3$mExG_mean),length(data_NTF_4$mExG_mean), length(data_NTF_5$mExG_mean)))
mExG <- c(data_NTF_1$mExG_mean, data_NTF_2$mExG_mean, data_NTF_3$mExG_mean, data_NTF_4$mExG_mean, data_NTF_5$mExG_mean)

df27 <- data.frame(TreatmentGroups_04032024, mExG)
boxplot(mExG ~ TreatmentGroups_04032024, data = df27)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$NDVI_mean), length(data_NTF_2$NDVI_mean), length(data_NTF_3$NDVI_mean),length(data_NTF_4$NDVI_mean), length(data_NTF_5$NDVI_mean)))
NDVI <- c(data_NTF_1$NDVI_mean, data_NTF_2$NDVI_mean, data_NTF_3$NDVI_mean, data_NTF_4$NDVI_mean, data_NTF_5$NDVI_mean)

df28 <- data.frame(TreatmentGroups_04032024, NDVI)
boxplot(NDVI ~ TreatmentGroups_04032024, data = df28)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mNDVI_mean), length(data_NTF_2$mNDVI_mean), length(data_NTF_3$mNDVI_mean),length(data_NTF_4$mNDVI_mean), length(data_NTF_5$mNDVI_mean)))
mNDVI <- c(data_NTF_1$mNDVI_mean, data_NTF_2$mNDVI_mean, data_NTF_3$mNDVI_mean, data_NTF_4$mNDVI_mean, data_NTF_5$mNDVI_mean)

df29 <- data.frame(TreatmentGroups_04032024, mNDVI)
boxplot(mNDVI ~ TreatmentGroups_04032024, data = df29)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$GNDVI_mean), length(data_NTF_2$GNDVI_mean), length(data_NTF_3$GNDVI_mean),length(data_NTF_4$GNDVI_mean), length(data_NTF_5$GNDVI_mean)))
GNDVI <- c(data_NTF_1$GNDVI_mean, data_NTF_2$GNDVI_mean, data_NTF_3$GNDVI_mean, data_NTF_4$GNDVI_mean, data_NTF_5$GNDVI_mean)

df30 <- data.frame(TreatmentGroups_04032024, GNDVI)
boxplot(GNDVI ~ TreatmentGroups_04032024, data = df30)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mGNDVI_mean), length(data_NTF_2$mGNDVI_mean), length(data_NTF_3$mGNDVI_mean),length(data_NTF_4$mGNDVI_mean), length(data_NTF_5$mGNDVI_mean)))
mGNDVI <- c(data_NTF_1$mGNDVI_mean, data_NTF_2$mGNDVI_mean, data_NTF_3$mGNDVI_mean, data_NTF_4$mGNDVI_mean, data_NTF_5$mGNDVI_mean)

df31 <- data.frame(TreatmentGroups_04032024, mGNDVI)
boxplot(mGNDVI ~ TreatmentGroups_04032024, data = df31)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$GRVI_mean), length(data_NTF_2$GRVI_mean), length(data_NTF_3$GRVI_mean),length(data_NTF_4$GRVI_mean), length(data_NTF_5$GRVI_mean)))
GRVI <- c(data_NTF_1$GRVI_mean, data_NTF_2$GRVI_mean, data_NTF_3$GRVI_mean, data_NTF_4$GRVI_mean, data_NTF_5$GRVI_mean)

df32 <- data.frame(TreatmentGroups_04032024, GRVI)
boxplot(GRVI ~ TreatmentGroups_04032024, data = df32)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mGRVI_mean), length(data_NTF_2$mGRVI_mean), length(data_NTF_3$mGRVI_mean),length(data_NTF_4$mGRVI_mean), length(data_NTF_5$mGRVI_mean)))
mGRVI <- c(data_NTF_1$mGRVI_mean, data_NTF_2$mGRVI_mean, data_NTF_3$mGRVI_mean, data_NTF_4$mGRVI_mean, data_NTF_5$mGRVI_mean)

df33 <- data.frame(TreatmentGroups_04032024, mGRVI)
boxplot(mGRVI ~ TreatmentGroups_04032024, data = df33)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$NDRE_mean), length(data_NTF_2$NDRE_mean), length(data_NTF_3$NDRE_mean),length(data_NTF_4$NDRE_mean), length(data_NTF_5$NDRE_mean)))
NDRE <- c(data_NTF_1$NDRE_mean, data_NTF_2$NDRE_mean, data_NTF_3$NDRE_mean, data_NTF_4$NDRE_mean, data_NTF_5$NDRE_mean)

df34 <- data.frame(TreatmentGroups_04032024, NDRE)
boxplot(NDRE ~ TreatmentGroups_04032024, data = df34)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mNDRE_mean), length(data_NTF_2$mNDRE_mean), length(data_NTF_3$mNDRE_mean),length(data_NTF_4$mNDRE_mean), length(data_NTF_5$mNDRE_mean)))
mNDRE <- c(data_NTF_1$mNDRE_mean, data_NTF_2$mNDRE_mean, data_NTF_3$mNDRE_mean, data_NTF_4$mNDRE_mean, data_NTF_5$mNDRE_mean)

df35 <- data.frame(TreatmentGroups_04032024, mNDRE)
boxplot(mNDRE ~ TreatmentGroups_04032024, data = df35)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$REGI_mean), length(data_NTF_2$REGI_mean), length(data_NTF_3$REGI_mean),length(data_NTF_4$REGI_mean), length(data_NTF_5$REGI_mean)))
REGI <- c(data_NTF_1$REGI_mean, data_NTF_2$REGI_mean, data_NTF_3$REGI_mean, data_NTF_4$REGI_mean, data_NTF_5$REGI_mean)

df36 <- data.frame(TreatmentGroups_04032024, REGI)
boxplot(REGI ~ TreatmentGroups_04032024, data = df36)



TreatmentGroups_04032024 <- rep(c("Treatment 1", "Treatment 2", "Treatment 3", "Treatment 4", "Treatment 5"),
                                c(length(data_NTF_1$mREGI), length(data_NTF_2$mREGI), length(data_NTF_3$mREGI),length(data_NTF_4$mREGI), length(data_NTF_5$mREGI)))
mREGI <- c(data_NTF_1$mREGI_mean, data_NTF_2$mREGI, data_NTF_3$mREGI, data_NTF_4$mREGI, data_NTF_5$mREGI)

df37 <- data.frame(TreatmentGroups_04032024, mREGI)
boxplot(mREGI ~ TreatmentGroups_04032024, data = df37)







#filter Treatment by Rep?

