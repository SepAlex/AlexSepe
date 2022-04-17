

library(nnfor)
library(forecast)
library(vars)
library(tswge)
library(naniar)
library(DataExplorer)
library(GGally)


hk_eco <- read.csv(file.choose(), header = TRUE)

head(hk_eco)


##### Doing EDA on data

hk_eda <- hk_eco

str(hk_eda)

summary(hk_eda)

ggpairs(hk_eda[2:14])



#choosing hk_eda$HSI...volume to run Univarate Model
hk <- (hk_eco$HSI...volume)

#scale data by dividing by 1 million to normalize data
hk.scale <- hk_eco$HSI...volume/1000000

#plotting different timeframes in data to game insight
plotts.sample.wge(hk.scale)# True Acf is slowly dumping  with spectral density showing evidence of a wandering behavior at peak f=0 as seen in the realization
# with a hint of cyclic behavior gives an indication of a dominant (1-B) present hence we proceed to difference the data


# Taking First difference to remove (1-B) from full  scaled data
hk.df1= artrans.wge(hk.scale,phi.tr = c(rep(0,364),1)) # differencing to remove seasonality to make the data more stationary

hk.df2 = artrans.wge(hk.df1,phi.tr = 1)#differencing the data to remove (1-b) factor and the sample auto-correlation appears not to look white noise  and but shows cyclic behavior.

# we will proceed to model the data

# Model identification for the transformed or difference data using AIC and BIC approach
aic5.wge(hk.df1, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.df1, p= 0:13, q= 0:5, type = "bic")

# Since AIC value is lower compared to that of BIC, we use ARMA(3,4) from the AIC results

#ljung.wge(hk.df2)$pval #FTR Ho
#ljung.wge(hk.df2, K = 48)$pval #FTR Ho

# Use estimate p & q to get estimates of phi and thetas
est.hk <- est.arma.wge(hk.df1,p=2, q=2)

est.hk$avar

mean(hk.scale) #1492.191

#Final Model
# (1-B^12)(1 - 1.7383B + 0.7394B^2)(X_t + 1492.191) = (1 - 1.3420B + 0.3625B^2)a_t  $avar = 538864.1



###### Short and Long term forecast

#ARMA(2,2)
ped.scale.sht1 = fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
plot(ped.scale.sht1$f, type = "l")
plot(seq(1,4233,1), hk.scale, type = "l",xlim = c(0,4598), ylab = "HSI_volume", main = "HSI_volume 365 days Short Term Forecast Using ARMA(2,2)")
lines(seq(4234,4598,1), ped.scale.sht1$f, type = "l", col = "red") # stationary models have natural attraction towards the mean

#ARIMA(2,1,2)
ped.scale.sht2  = fore.aruma.wge(hk.scale, s=365,est.hk$phi,theta = est.hk$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
plot(ped.scale.sht2$f, type = "l")
plot(seq(1,4233,1), hk.scale, type = "l",xlim = c(0,4598), ylab = "HSI_volume", main = "HSI_volume 365 days Short Term Forecast Using ARIMA(2,1,2)")
lines(seq(4234,4598,1), ped.scale.sht2$f, type = "l", col = "red")


#ARMA(2,2)
ped.scale.lng1 = fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
plot(ped.scale.lng1$f, type = "l")
plot(seq(1,4233,1), hk.scale, type = "l",xlim = c(0,6058), ylab = "HSI_volume", main = "HSI_volume 5 Years Long Term Forecast Using ARMA(2,2)")
lines(seq(4234,6058,1), ped.scale.lng1$f, type = "l", col = "red") # stationary models have natural attraction towards the mean

#ARIMA(2,1,2)
ped.scale.lng2 = fore.aruma.wge(hk.scale, s=365,est.hk$phi,theta = est.hk$theta, n.ahead = 1825, lastn = TRUE,limits = FALSE)
plot(ped.scale.lng2$f, type = "l")
plot(seq(1,4233,1), hk.scale, type = "l",xlim = c(0,6058), ylab = "HSI_volume", main = "HSI_volume 5 Years Long Term Forecast Using ARIMA(2,1,2)")
lines(seq(4234,6058,1), ped.scale.lng2$f, type = "l", col = "red")



###### Calculating ASE for short and Long terms observations

#Calculating ASE for short term 365 days observations ARMA Model

n = 365
x = hk.scale
x.pred.scale1 = fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
ASE.scale1 = mean((x[(length(x)-n+1):(length(x))]-x.pred.scale1$f)^2)
print(paste0('Short Term ARMA(2,2) ASE: ',round(ASE.scale1,2)))

# Short Term ARMA ASE: 271120.67

# Calculating ASE for short Terms 365 days observations of ARIMA model
n = 365
x = hk.scale
x.pred.scale2 = fore.aruma.wge(hk.scale,s=365,est.hk$phi,theta = est.hk$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
ASE.scale2 = mean((x[(length(x)-n+1):(length(x))]-x.pred.scale2$f)^2)
print(paste0('Short Term ARIMA(2,1,2) ASE: ',round(ASE.scale2,2)))

# Short Term ARIMA ASE Seasonal: 18729.93

# Calculating ASE for Long Term 1825 days (5 years) observations of ARMA model
n1 = 1825
x = hk.scale
x.pred.scale3 = fore.arma.wge(hk.scale, phi = est.hk$phi,theta = est.hk$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
ASE.scale3 = mean((x[(length(x)-n1+1):(length(x))]-x.pred.scale3$f)^2)
print(paste0('Long Term ARMA(2,2) ASE: ',round(ASE.scale3,2)))

#Long Term ARMA ASE: 2413770.67

# Calculating ASE for Long Term 1825 days (5 years) observations of ARIMA model
n1 = 1825
x = hk.scale
x.pred.scale4 = fore.aruma.wge(hk.scale,s=365,est.hk$phi,theta = est.hk$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
ASE.scale4 = mean((x[(length(x)-n1+1):(length(x))]-x.pred.scale4$f)^2)
print(paste0('Long Term ARIMA(2,1,2) ASE: ',round(ASE.scale4,2)))



#Long Term ARIMA ASE: 538957.77


#Rolling Window ASE

Rolling_Window_ASE = function(series, trainingSize, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
{
  ASEHolder = numeric()

  for( i in 1:(length(series)-(trainingSize + horizon) + 1))
  {

    forecasts = fore.aruma.wge(series[i:(i+(trainingSize-1))],phi = phis, theta = thetas, s = s, d = d,n.ahead = horizon)

    ASE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts$f)^2)

    ASEHolder[i] = ASE

  }

  ASEHolder
  hist(ASEHolder)
  WindowedASE = mean(ASEHolder)

  print("The Summary Statistics for the Rolling Window ASE Are:")
  print(summary(ASEHolder))
  print(paste("The Rolling Window ASE is: ",WindowedASE))
  return(WindowedASE)
}

#ARMA Rolling ASE
Rolling_Window_ASE(hk.scale, phis = est.hk$phi,thetas = est.hk$theta, horizon = 365, trainingSize = 500)

#The Short Term Rolling Window ASE is: 645583.3


#ARIMA rOLLING ASE
Rolling_Window_ASE(hk.scale,d=1, s = 365, phis = est.hk$phi,thetas = est.hk$theta, horizon = 1825, trainingSize = 1000)

#The Long Term Rolling Window ASE is: 1558748



#Compare Generated Realizations

hk.arma.gen = gen.arma.wge(1000,phi = est.hk$phi,theta = est.hk$theta, vara = est.hk$avar)

plotts.sample.wge(hk.arma.gen)


hk.arima.gen = gen.arima.wge(1000,d=1,s=365, phi = est.hk$phi,theta = est.hk$theta, vara = est.hk$avar)

plotts.sample.wge(hk.arima.gen)




###################### using Private.Domestic..Price.Index ##########################################################

hk.ind <- hk_eda$Private.Domestic..Price.Index.


#plotting different time frames in data to game insight
plotts.sample.wge(hk.ind)# realization is shows wandering behavior with no sign of cyclic behavior
# it also has its True sample autocorrelation showing identical rho_k = 1
#spectral density shows a dominant (1-B) behavior at peak f=0 hence an indication for an ARIMA model we then proceed to take the first different

# Taking First difference to remove (1-B^12) from full data
hk.ind.df1= artrans.wge(hk.ind,phi.tr = c(rep(0,364),1))# taking the 12th difference of the data (1-B^12) to stationarize the data but shows seasonality at lag 12
#hence we remove the seasonality in the data with this differencing like in the airline model to make more stationary

acf(hk.ind.df1) # Acfs are slowing dumping with some dominant (1-B) Factor

hk.ind.df2= artrans.wge(hk.ind.df1,phi.tr = 1) #taking the second difference of the first difference to remove the (1-B) factor shows significance at lag 1 with hint of cyclic behavior

acf(hk.ind.df2)


plotts.sample.wge(hk.ind.df2) # sample Autocorrelations shows cyclic behavior at lag 1 and a slight sinusoidal behavior damping in the sample auto-correlation
#spectral density of the realization has a high peak at f=0 but has ripples between f=0 and f=0.5 worth noting

parzen.wge(hk.ind.df2) # wander with a strong indication of cyclic behavior

# Find a model for the transformed or difference data using AIC and BIC with the assumption data is stationary
aic5.wge(hk.ind.df2, p= 0:13, q= 0:5, type = "aic")

aic5.wge(hk.ind.df2, p= 0:13, q= 0:5, type = "bic")


# Use estimate p & q to get estimates of phi and thetas
est.hk.ind <- est.arma.wge(hk.ind.df2,p=11, q=4)

est.hk.ind$avar

mean(hk.ind)



#Final Model
# (1+1.1624B+0.9835B^2 (1 + 0.0322B - 0.3918B^2 - 1.1102B^3 + 0.0967B^4 + 0.0845B^5 + 0.0623B^6 + 0.0676B^7 + 0.0547B^8 + 0.0711B^9 + 0.0540B^10 + 0.0511B^11)(X_t + 197.065) = (1 - 0.7450B - 0.3282B^2 - 0.7976B^3 + 0.9230B^4)a_t  $avar = 0.1306278


###### Short and Long term forecast

#ARMA(11,4)
ped.ind.sht1 = fore.arma.wge(hk.ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
plot(ped.ind.sht1$f, type = "l")
plot(seq(1,4233,1), hk.ind, type = "l",xlim = c(0,4598), ylab = "Private.Domestic.Price.Index", main = "Private.Domestic.Price.Index 365 days Short Term Forecast using ARMA(11,4)")
lines(seq(4234,4598,1), ped.ind.sht1$f, type = "l", col = "red")

#ARIMA(11,2,4)
ped.ind.sht2  = fore.aruma.wge(hk.ind,d=1,s=365,est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
plot(ped.ind.sht2$f, type = "l")
plot(seq(1,4233,1), hk.ind, type = "l",xlim = c(0,4598), ylab = "Private.Domestic.Price.Index", main = "Private.Domestic.Price.Index 365 days Short Term Forecast using ARIMA(11,2,4)")
lines(seq(4234,4598,1), ped.ind.sht2$f, type = "l", col = "red")


#ARMA(11,4)
ped.ind.lng1 = fore.arma.wge(hk.ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
plot(ped.ind.lng1$f, type = "l")
plot(seq(1,4233,1), hk.ind, type = "l",xlim = c(0,6058), ylab = "Private.Domestic.Price.Index", main = "Private.Domestic.Price.Index 5 Years Long Term Forecast Using ARMA(11,4)")
lines(seq(4234,6058,1), ped.ind.lng1$f, type = "l", col = "red")

#ARIMA(11,2,4)
ped.ind.lng2 = fore.aruma.wge(hk.ind,d= 1, s=365,est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
plot(ped.ind.lng2$f, type = "l")
plot(seq(1,4233,1), hk.ind, type = "l",xlim = c(0,6058), ylab = "Private.Domestic.Price.Index", main = "Private.Domestic.Price.Index 5 Years Long Term Forecast Using ARIMA(11,2,4)")
lines(seq(4234,6058,1), ped.ind.lng2$f, type = "l", col = "red")



##### Calculating ASE for short and Long terms observations

#Calculating ASE for short term 365 days observations

# Calculating ASE for Short Term 365 days observations of ARMA model
n = 365
x = hk.ind
x.pred.ind1 = fore.arma.wge(hk.ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
ASE.ind1 = mean((x[(length(x)-n+1):(length(x))]-x.pred.ind1$f)^2)
print(paste0('Short Term ARMA(11,4) ASE: ',round(ASE.ind1,2)))

# Short Term ARMA(11,4) ASE: 29985.79

# Calculating ASE for Short Term 365 days observations of ARIMA model
n = 365
x = hk.ind
x.pred.ind2 = fore.aruma.wge(hk.ind,s=365,est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 365, lastn = FALSE,limits = FALSE)
ASE.ind2 = mean((x[(length(x)-n+1):(length(x))]-x.pred.ind2$f)^2)
print(paste0('Short Term ARIMA(11,2,4) ASE: ',round(ASE.ind2,2)))

# Short Term ARIMA(11,2,4) ASE: 104.97


# Calculating ASE for Long Term 1825 days observations of ARMA model
n1 = 1825
x = hk.ind
x.pred.ind3 = fore.arma.wge(hk.ind, phi = est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
ASE.ind3 = mean((x[(length(x)-n1+1):(length(x))]-x.pred.ind3$f)^2)
print(paste0('Long Term ARMA(11,4) ASE: ',round(ASE.ind3,2)))

#Long Term ARMA(11,4) ASE: 14794.49


# Calculating ASE for Long Term 1825 days observations of ARIMA model
n1 = 1825
x = hk.ind
x.pred.ind4 = fore.aruma.wge(hk.ind,s=365,est.hk.ind$phi,theta = est.hk.ind$theta, n.ahead = 1825, lastn = FALSE,limits = FALSE)
ASE.ind4 = mean((x[(length(x)-n1+1):(length(x))]-x.pred.ind4$f)^2)
print(paste0('Long Term ARIMA(11,2,4) ASE: ',round(ASE.ind4,2)))



#Long Term ARIMA(11,2,4) ASE: 8384.81


#Rolling Window ASE

Rolling_Window_ASE = function(series, trainingSize, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
{
  ASEHolder = numeric()

  for( i in 1:(length(series)-(trainingSize + horizon) + 1))
  {

    forecasts = fore.aruma.wge(series[i:(i+(trainingSize-1))],phi = phis, theta = thetas, s = s, d = d,n.ahead = horizon)

    ASE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts$f)^2)

    ASEHolder[i] = ASE

  }

  ASEHolder
  hist(ASEHolder)
  WindowedASE = mean(ASEHolder)

  print("The Summary Statistics for the Rolling Window ASE Are:")
  print(summary(ASEHolder))
  print(paste("The Rolling Window ASE is: ",WindowedASE))
  return(WindowedASE)
}

#ARMA Rolling ASE
Rolling_Window_ASE(hk.ind, phis = est.hk.ind$phi,thetas = est.hk.ind$theta, horizon = 365, trainingSize = 500)

#The Short Term Rolling Window ASE is: 1677.426


#ARIMA rOLLING ASE
Rolling_Window_ASE(hk.ind, d=2, s = 365, phis = est.hk.ind$phi,thetas = est.hk.ind$theta, horizon = 1825, trainingSize = 1000)

#The Short Term Rolling Window ASE is: 38524175


### Testing parameterized model estimates

#Compare Generated Realizations and Spectral Densities

hk.ind.arma.gen = gen.arma.wge(5000,phi = est.hk.ind$phi,theta = est.hk.ind$theta, vara = est.hk.ind$avar)

plotts.sample.wge(hk.ind.arma.gen)


hk.ind.arima.gen = gen.arima.wge(5000,d=2, s=365, phi = est.hk.ind$phi,theta = est.hk.ind$theta, vara = est.hk.ind$avar)

plotts.sample.wge(hk.ind.arima.gen)



#Compare Spectral Densities using parameters from ARIMA(11,2,4) model generated
sims = 2
SpecDen = parzen.wge(hk.ind, plot = "FALSE")
plot(SpecDen$freq,SpecDen$pzgram, type = "l", lwd = 2)

for( i in 1: sims)
{
  SpecDen2 = parzen.wge(gen.arima.wge(10000,d=2, s=365, phi = est.hk.ind$phi,theta = est.hk.ind$theta, vara = est.hk.ind$avar, plot = "FALSE"), plot = "FALSE")
  lines(SpecDen2$freq,SpecDen2$pzgram, lwd = 2, col = "red")
}

## Estimated parameters from generated model are pretty consistent with with the actual parameters from the data


####### end of code

#########################################################################################################################

####### var Model

#########################################################################################################################
library(RColorBrewer)
library(vars)



hk_vol <- hk_eco$HSI...volume
hk_ind <- hk_eco$Private.Domestic..Price.Index.
hk_sal_amt <- hk_eco$First.hand.sales.amount




##### Taking First Difference of data variables
hk_vol_df1 <- artrans.wge(hk_vol,phi.tr = c(rep(0,364),1))
hk_ind_df1 <- artrans.wge(hk_ind,phi.tr = c(rep(0,364),1))
hk_sal_amt_df1 <- artrans.wge(hk_sal_amt,phi.tr = c(rep(0,364),1))


######## Taking Second difference of data Variables
hk_vol_df2 <- artrans.wge(hk_vol_df1,phi.tr = 1)
hk_ind_df2 <- artrans.wge(hk_ind_df1,phi.tr = 1)
hk_sal_amt_df2 <- artrans.wge(hk_sal_amt_df1,phi.tr = 1)

################################## VAR Model for Short Term

#VARSelect on Differenced Data chooses BIC = 12
VARselect(cbind(hk_vol_df2, hk_ind_df2, hk_sal_amt_df2),lag.max = 20, type = "both")

###BIC= (SC(n)) = 7.205915e+01 for p = 12

#VAR with p = 11
hk.VAR.st = VAR(cbind(hk_vol_df2, hk_ind_df2, hk_sal_amt_df2),type = "both",p = 11)


hk.pred.st = predict(hk.VAR.st,n.ahead=365,lastn= TRUE, ci = 0.95)


#We have predicted differences .... calculate actual Hong Kong HSI Volune  
startingpoint.st = hk_vol[3869:4233]
hk.vol.fore.st = hk.pred.st$fcst$hk_vol[,1:3] + startingpoint.st


#Plot
plot(seq(1,4233,1), hk_vol, type = "l",xlim = c(0,4233),  xlab = "Number of Housing", ylab = "Hong Kong HSI Volume", main = "VAR Model Short Term Forecast for Hong Kong Housing Index Volume")
lines(seq(3869,4233,1), hk.vol.fore.st[,1], type = "l", col = "red")

fanchart(hk.pred.st, colors = brewer.pal(n = 8, name = "Blues"))


#Find Short Term VAR ASE 


VAR.ASE.st = mean((hk_vol[3869:4233] - hk.vol.fore.st[,1])^2)

VAR.ASE.st 


################################# VAR Model for Long Term 

#VARSelect on Differenced Data chooses BIC = 12
VARselect(cbind(hk_vol_df2, hk_ind_df2, hk_sal_amt_df2),lag.max = 15, type = "both")

#BIC= 7.213653e+01 for p = 11

#VAR with p = 15
hk.VAR.lg = VAR(cbind(hk_vol_df2, hk_ind_df2, hk_sal_amt_df2),type = "both",p = 11)

hk.pred.lg = predict(hk.VAR.lg,n.ahead=1825,lastn= TRUE, ci = 0.95)


#We have predicted differences .... calculate actual Hong Kong HSI Volune  
startingpoint.lg = hk_vol[2409:4233]
hk_hsi_for.lg = hk.pred.lg$fcst$hk_vol[,1:3] + startingpoint.lg


#Plot

plot(seq(1,4233,1), hk_vol, type = "l",xlim = c(0,4233), xlab = "Number of Housing", ylab = "Hong Kong HSI Volume", main = "VAR Model Long Term Forecast for Hong Kong Housing Index Volume")
lines(seq(2409,4233,1), hk_hsi_for.lg[,1], type = "l", col = "blue")

fanchart(hk.pred.lg, colors = brewer.pal(n = 8, name = "Blues")) # Change color pallet to make distinguishable.



#Find Long Term VAR ASE 


VAR.ASE.lg = mean((hk_vol[2409:4233] - hk_hsi_for.lg[,1])^2)

VAR.ASE.lg 




################################################### End of Code for VAR Model ########################################################

######################################################################################################################################

########################################################## MLP Model #################################################################

####### Data Prep for Selected Features from Hong Kong Housing Market Data
hk_eco1 <- hk_eco[2:14]

hk.st <- hk_eco[1:3867,]

hk.lg <- hk_eco[1:2408,]

hk.all <- hk_eco[1:4233,]



hk.all.df <- data.frame(hk.ind=ts(hk.all$Private.Domestic..Price.Index.), hk.sale.amt=ts(hk.all$First.hand.sales.amount))

hk.mlp.df.st <- data.frame(hk.ind=ts(hk.st$Private.Domestic..Price.Index.), hk.sale.amt=ts(hk.st$First.hand.sales.amount))

hk.mlp.df.lg <- data.frame(hk.ind=ts(hk.lg$Private.Domestic..Price.Index.), hk.sale.amt=ts(hk.lg$First.hand.sales.amount))

?mlp()

ts.hk <- ts(hk.st$HSI...volume)
###### Fitting the MLP Model on Hong Kong Housing Index Volume for Short Term 
fit.hk.mlp.st <- mlp(ts.hk,difforder= c(1,365), reps = 10,comb = 'mean', xreg = hk.mlp.df.st )  

fit.hk.mlp.st 

#### Plot the fitted model on Hong Kong Housing Index Volume using MLP Approach Short Term
plot(fit.hk.mlp.st)


tss.hk <- ts(hk.lg$HSI...volume)
###### Fitting the MLP Model on Hong Kong Housing Index Volume for Long Term 
fit.hk.mlp.lg <- mlp(tss.hk,difforder= c(1,365),reps = 20,comb = 'mean', xreg = hk.mlp.df.lg )  

fit.hk.mlp.lg  



####### Plot the fitted model on Hong Kong Housing Index Volume using MLP Approach Long Term 
plot(fit.hk.mlp.lg)


######## Short Term Forecast using MLP Approach
hk.fore.st = forecast(fit.hk.mlp.st, h = 365,xreg = hk.all.df)

plot(hk.fore.st)

### Long Term Forecast using MLP Approach
hk.fore.lg = forecast(fit.hk.mlp.lg, h = 1825, xreg = hk.all.df)

plot(hk.fore.lg)



########### Calculating short term ASE for Hong Kong House Index Volume
ASE.mlp.st = mean((hk_eco$HSI...volume[3869:4233] - hk.fore.st$mean)^2)

ASE.mlp.st


##### Calculating long term ASE for Hong Kong House Index Volume
ASE.mlp.lg = mean((hk_eco$HSI...volume[2409:4233] - hk.fore.lg$mean)^2)

ASE.mlp.lg


#### Plot for short term forecast 
plot(seq(1,4233,1), hk_eco$HSI...volume, type = "l",xlim = c(0,4233), xlab = "Number of Housing", ylab = "Hong Kong HSI Volume", main = "MLP Short Term Forecast for Hong Kong Housing Index Volume ")
lines(seq(3869,4233,1), hk.fore.st$mean, type = "l", col = "red")


#### Plot for long term forecast 
plot(seq(1,4233,1), hk_eco$HSI...volume, type = "l",xlim = c(0,4233), xlab = "Number of Housing", ylab = "Hong Kong HSI Volume", main = "MLP Long Term Forecast for Hong Kong Housing Index Volume")
lines(seq(2409,4233,1), hk.fore.lg$mean, type = "l", col = "purple")


########################## Short Term Ensemble #####################################################

ensemble.st  = (hk.vol.fore.st[,1] + hk.fore.st$mean)/2


##### Short Ensemble ASE
ASE.ens.st = mean((hk_eco$HSI...volume[3869:4233] - ensemble.st)^2)

ASE.ens.st

############## Plot Short Term Ensemble
plot(hk_eco$HSI...volume, type = "l")
lines(seq(3869,4233,1),ensemble.st,col = "green")


##############################Plot for Short Term Ensemble Forecast 
plot(seq(1,4233,1), hk_eco$HSI...volume, type = "l",xlim = c(0,4233), ylab = "Hong Kong HSI Volume", main = "Short Term Ensemble Forecast")
lines(seq(3869,4233,1), ensemble.st, type = "l", col = "red")


############################### Long Term Ensemble #####################################################

ensemble.lg  = (hk_hsi_for.lg[,1] + hk.fore.lg$mean)/2


######### Long Ensemble ASE
ASE.ens.lg = mean((hk_eco$HSI...volume[3869:4233] - ensemble.lg)^2)

ASE.ens.lg

#### Plot Short Term Ensemble
plot(hk_eco$HSI...volume, type = "l")
lines(seq(2409,4233,1),ensemble.lg,col = "green")


########### Plot for Short Term Ensemble Forecast 
plot(seq(1,4233,1), hk_eco$HSI...volume, type = "l",xlim = c(0,4233), ylab = "Hong Kong HSI Volume", main = "Long Term Ensemble Forecast")
lines(seq(2409,4233,1), ensemble.lg, type = "l", col = "purple")

