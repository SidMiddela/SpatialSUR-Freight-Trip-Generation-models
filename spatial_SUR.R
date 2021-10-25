r2pool = vector()
ipool = vector()
AICpool = vector()
LMpool=vector()
rhopool=vector()
for (i in 2:50) ## loop for finding the best fit by i
{
library(spsur)

###data preparation

source("D:/R Codes from Lab system/packages.R")
source("D:/R Codes from Lab System/FA_data.r")
att$leased_all <- att$leased_3w + att$leased_bi + att$leased_pi + att$leased_tr
  
att1 <- data.frame(1:1074)

att1$lon <- att$lon
att1$lat <- att$lat
att1$to_b3 <- att$to_b3
att1$to_pi <- att$to_pi
att1$to_tr <- att$to_tr
att1$C4772 <-att$C4772
att1$C47713 <-att$C47713
att1$owned_3w <- att$owned_3w
att1$owned_bi <- att$owned_bi
att1$owned_tr <- att$owned_tr
att1$leased_tr <- att$leased_tr
att1$i_ea <- att$i_ea
att1$est_area <- att$est_area
att1$home_del_pro <- att$home_del_pro
att1$cons <- att$cons
att1$RetAss <- att$RetAss
att1$oper_age <- att$oper_age
att1$oper_time <- att$oper_time
att1$i_ea <- att$i_ea
att1$i_emp <- as.factor(att$i_emp)
att1$leased_pi <- att$leased_pi
att1$owned_pi <- att$owned_pi
att1$aload_docks <- att$aload_docks
att1$no_load_docks <- att$no_load_docks
att1$off_street_park   <- att$off_street_park
att1$C5610 <-att$C5610
att1$C4752 <-att$C4752
att1$ind <- att$ind
att1$owned_all <- att1$owned_3w + att1$owned_bi + att1$owned_pi + att1$owned_tr
att1$leased_all <- att$leased_all

att1 <- att1[-1]
att1 <- att1[complete.cases(att1),] ### Pair-wise deletion to ensure complete data - required for SUR models

lo <- att1$lon
la <- att1$lat

Y1 <- (att1$to_b3)
Y2 <- (att1$to_pi)
Y3 <- (att1$to_tr)


library(spdep)
library(spsur)

#### Code for obtaining the nearest neighbour matrix
xy=cbind(lo, la)
nn4 <- knearneigh(xy, k=i, longlat = TRUE)
nc.nn4 <- knn2nb(nn4)
list.nc.nn4 <- nb2listw(nc.nn4, style="W", zero.policy=T)
# mat.nc.nn4 <- listw2mat(list.nc.nn4) ## can be used to convert the list object to matrix

## SpSUR model formula - Y1 represents FTA from bikes and ThWs, Y2 -pickups and Y3 - trucks
Tformula <- Y1 | Y2 | Y3 ~ owned_bi+owned_3w+off_street_park| 
  C47713+leased_pi+i_ea+
  owned_pi+Fst+off_street_park| owned_tr +
  leased_tr + cons + RetAss + est_area +C5610
LMs <- lmtestspsur(Tformula, data = att1, listw = list.nc.nn4)
LMs

spcsur.slm.to <- spsurml(formula=Tformula, data = att1, type = "slm", listw = list.nc.nn4)
summary(spcsur.slm.to)
AIC(spcsur.slm.to)

ipool <- append(ipool, i, 
                 after=length(ipool))
AICpool <- append(AICpool, round(AIC(spcsur.slm.to),0), 
                 after=length(AICpool))
r2pool <- append(r2pool, round(as.numeric(spcsur.slm.to$R2[1]),3), 
       after=length(r2pool))
print(round(as.numeric(spcsur.slm.to$R2[1]),3))
LMpool <- append(LMpool,paste(round(as.numeric(LMs[[1]][3]),2),
      round(as.numeric(LMs[[2]][3]),2),
      round(as.numeric(LMs[[3]][3]),2),
      round(as.numeric(LMs[[4]][3]),2),
      round(as.numeric(LMs[[5]][3]),2)),after=length(LMpool))
rhopool <- append(rhopool, paste(round(as.numeric(spcsur.slm.to[[8]][1]),2),
  round(as.numeric(spcsur.slm.to[[8]][2]),2),
  round(as.numeric(spcsur.slm.to[[8]][3]),2)),after=length(rhopool))
}

spcsur.sim <- spsurml(formula=Tformula, data = att1, type = "sim", listw = list.nc.nn4)
summary(spcsur.sim)
AIC(spcsur.sim)
spsur::impactspsur(spcsur.sim)

## Spatial Multiple Linear Regression Models for FTA bikes and three-wheelers
OLS_b3 <- lagsarlm(to_b3~owned_3w+owned_bi+off_street_park,
                   data = att1, listw = list.nc.nn4, tol.solve = 1e-30)
summary(OLS_b3)

## Non-spatial Multiple Linear Regression Models for FTA bikes and three-wheelers
OLS_Ab3 <- lm(to_b3~owned_3w+owned_bi+C4772+off_street_park,
                   data = att1)
summary(OLS_Ab3)

impacts.sarlm(OLS_b3, listw=list.nc.nn4)

## Spatial Multiple Linear Regression Models for FTA pickups
OLS_pi <- lagsarlm(to_pi ~ C47713+leased_pi+i_ea+owned_pi+
                     Fst+C4772+off_street_park,
                   data = att1, listw = list.nc.nn4, tol.solve = 1e-30)
summary(OLS_pi)
impacts.sarlm(OLS_pi, listw=list.nc.nn4)

## Non-spatial Multiple Linear Regression Models for FTA pickups
OLS_Api <- lm(to_pi ~ C47713+leased_pi+i_ea+owned_pi+
                     Fst+C4772+off_street_park,
                   data = att1)
summary(OLS_Api)

## Spatial Multiple Linear Regression Models for FTA trucks
OLS_tr <- lagsarlm(to_tr~owned_tr +
                     leased_tr + cons + RetAss + est_area +C5610,data=att1,
                   listw = list.nc.nn4, tol.solve = 1e-30)
summary(OLS_tr)
impacts.sarlm(OLS_tr, listw=list.nc.nn4)

## Non-spatial Multiple Linear Regression Models for FTA trucks
OLS_Atr <- lm(to_tr~owned_tr +
                leased_tr + cons + RetAss + est_area +C5610,data=att1)

summary(OLS_Atr)

rmse(OLS_Ab3$fitted.values, OLS_Ab3$model$to_b3)
rmse(OLS_Api$fitted.values, OLS_Api$model$to_pi)
rmse(OLS_Atr$fitted.values, OLS_Atr$model$to_tr)

mean(att1$to_b3)
mean(att1$to_pi)
mean(att1$to_tr)

library(Metrics)

## Calculation of goodness of fit - RMSE

bikesRMSEsur <- data.frame(spcsur.slm.to$fitted.values[1:1021])
names(bikesRMSEsur) <- c("fit")
bikesRMSEsur$Y <- Y1
bikesRMSEsur <- bikesRMSEsur[abs(bikesRMSEsur$fit - bikesRMSEsur$Y)<10,]
rmse(bikesRMSEsur$fit, bikesRMSEsur$Y)

piRMSEsur <- data.frame(spcsur.slm.to$fitted.values[1022:2042])
names(piRMSEsur) <- c("fit")
piRMSEsur$Y <- Y2
piRMSEsur <- piRMSEsur[abs(piRMSEsur$fit - piRMSEsur$Y)<10,]
rmse(piRMSEsur$fit, piRMSEsur$Y)

trRMSEsur <- data.frame(spcsur.slm.to$fitted.values[2043:3063])
names(trRMSEsur) <- c("fit")
trRMSEsur$Y <- Y3
trRMSEsur <- trRMSEsur[abs(trRMSEsur$fit - trRMSEsur$Y)<10,]
rmse(trRMSEsur$fit, trRMSEsur$Y)

rmse(spcsur.slm.to$fitted.values[1:1021], Y1)
rmse(spcsur.slm.to$fitted.values[1022:2042], Y2)
rmse(spcsur.slm.to$fitted.values[2043:3063], Y3)