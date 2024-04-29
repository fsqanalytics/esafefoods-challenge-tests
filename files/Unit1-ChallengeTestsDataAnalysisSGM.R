## Fitting Secondary Growth Models
## Statistical Examination of Challenge Test Data and Modelling
## Read the listeria.csv with the read.csv() function and save it as data frame (df)

df <- read.csv("../data/listeria.csv", sep=";", header=TRUE)

## Check the dataset using the structure (str()) function
str(df)

## Check the data shape {background-image="../images/background.png" background-size="contain"}
plot(GR ~ Temp, data=df,
     xlab="Temperature (ºC)", ylab="Growth rate")

plot(sqrtGR ~ Temp, data=df,
     xlab="Temperature (ºC)", ylab="Growth rate")

## Install the predmicror package from GitHub repository
devtools::install_github("fsqanalytics/predmicror")

## Load the predmicror package
library(predmicror)

## Cardinal model for temperature
CMTI

## Starting values
start.values = list(Tmax=47, Tmin=0.1, MUopt=0.020, Topt=37)

## Fit using gsl_nls() function

library(predmicror)
library(gslnls)
fit <- gsl_nls(sqrtGR ~ CMTI(Temp,Tmax,Tmin,MUopt,Topt),
           data=df,
           start =  start.values
           )
fit

## Check fitting results {background-image="../images/background.png" background-size="contain"}
summary(fit)
## Extract the model coefficients using the coef() function
coefs <- coef(fit)
coefs

## Extract the confidence intervals of the parameters using the confint() function
confint(fit)

## Plot the fitted values
## Define a vector with auxiliary temperature data

new.temp=seq(0.7,44, by=0.2)

## Use the predict() function to compute the prediction interval
fits <- predict(fit,
                newdata = data.frame(Temp=new.temp),
                interval = "prediction", level = 0.95)
str(fits)

## Check the fits object
head(fits[, ])

## Plot prediction interval
plot(sqrtGR ~ Temp, data=df, ylim=c(0,.2))
lines(new.temp, fits[, 1], col="blue")
lines(new.temp, fits[, 2], col="red")
lines(new.temp, fits[, 3], col="red")

