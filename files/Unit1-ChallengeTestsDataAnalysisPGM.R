## Fitting Primary Growth Models
## Statistical Examination of Challenge Test Data and Modelling

## load the data set
df <- read.csv("./data/ecoli.csv", sep=";", header=TRUE)

## Check the first lines of the file
head(df)


## Check the structure of the data set
str(df)

## Check the levels of the variables: `Condition` and `Temp`
table(df$Condition)
table(df$Temp)
table(df$Condition, df$Temp)

## Check the data shape {background-image="../images/background.png" background-size="contain"}
library(ggplot2)
ggplot(data = df, aes(x=Time, y=lnN, group=factor(Temp), colour=factor(Temp))) +
  geom_point()

## Download and install the `R` software: https://cran.r-project.org/
## Install the `predmicror` package: https://fsqanalytics.github.io/predmicror/
library(devtools)
devtools::install_github("fsqanalytics/predmicror")

## Huang full model: https://fsqanalytics.github.io/predmicror/reference/HuangFM.html

## Defining the starting values
start.values = list(Y0=0.0, Ymax=22.0, MUmax=1.7, lag=5.0) 

## Fit using gsl_nls() function
library(predmicror)
library(gslnls)
fit <- gsl_nls(lnN ~ HuangFM(Time, Y0, Ymax, MUmax, lag),
           data=df[df$Condition==1 & df$Repetition==1, ],
           start =  start.values
)
fit

## Check fitting results
summary(fit)

## Extract the confidence intervals
confint(fit)

## Plot the fitted values
new.times=seq(0,24, by=0.1)

## Compute the confidence interval
fits <- predict(fit,
                newdata = data.frame(Time=new.times),
                interval = "prediction", level = 0.95)
str(fits)

head(fits[, ])

## Plot confidence interval
plot(lnN ~ Time, data=df[df$Condition==1 & df$Repetition==1, ], ylim=c(5,20))
lines(new.times, fits[, 1], col="blue")
lines(new.times, fits[, 2], col="red")
lines(new.times, fits[, 3], col="red")
