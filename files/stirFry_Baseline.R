library(extraDistr)
library(mc2d)
library(qraLm)
source("stirFry_MildHeat.R")
source("stirFry_Storage.R")

#########################################
### Generating the contamination matrix
#########################################
### Parameters
nLots <- 1000
sizeLot <- 5000 # Units: bag of frozen Vegs
unitSize <- 250 #g
P <- 0.1
C0MeanLog <- 1.0 #log10 CFU/g
C0SdLog <- 0.3
# distribution of initial concentration log CFU/g in frozen vegs at retail
hist (rnorm(10000, C0MeanLog, C0SdLog))

### Production of contamination matrix
set.seed(291) #applying LotGen(), a function of qraLm package
data <- LotGen(nLots=nLots, sizeLot=sizeLot, unitSize = unitSize,
               P=P, C0MeanLog=C0MeanLog,C0SdLog=C0MeanLog)

#########################################
### Stir frying of frozen vegs
#########################################
### Parameters
tempBlanch <- 70 #ºC
timeBlanch <- 3  #min

### Simulation
set.seed(101)
res1 <- stirFry_MildHeat(data,
                       tempBlanch = tempBlanch,
                       timeBlanch = timeBlanch,
                       logDrefMean=-1.79,
                       logDrefSd=0.247,
                       zT=5.891)

### Results (CFU per 250 g)
hist(res1$N)

### We must convert to CFU/g
Conc1 = res1$N/unitSize
hist(Conc1)
hist(log10(Conc1)) #log CFU/g
mean(Conc1)  #CFU/g
quantile(Conc1, probs=c(0.025, 0.50, 0.975)) #CFU/g

### However, a portion of 150 g is consumed
### So, let us express in CFU per 150 g
Conc1.Portion = Conc1*150
hist(Conc1.Portion)
hist(log10(Conc1.Portion), xlab="Dose (log10 CFU/portion)") #log CFU/150 g
mean(Conc1.Portion)  # CFU/150 g
quantile(Conc1.Portion, probs=c(0.025, 0.50, 0.975)) # CFU/150 g

#################################################
### Storage of left over of stir fried vegetables
################################################

### Parameters
Temp = 6  #ºC
time = 24 # hours

### Simulation
set.seed(102)
res2 <- stirFry_Storage(res1,
                        Temp = Temp,
                        time = time,
                        MPD = 8, #log CFU/g
                        Tmin = -1.18, #ºC
                        meanEGR5 = 0.0117,
                        sdEGR5 = 0.00816, 
                        servingSize = 100, #amount of left overs
                        pDefrost = 1)

### Results (CFU per 100 g of left over)
hist(res2$N)
hist(log10(res2$N), xlab="Dose (log10 CFU/portion)", xlim=c(0,4)) #log CFU/100 g
mean(res2$N)  #CFU/100 g
quantile(res2$N, probs=c(0.025, 0.50, 0.975)) #CFU/100 g

