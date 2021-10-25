library(spsur)

source("D:/R Codes from Lab system/packages.R")
source("D:/R Codes from Lab System/FP_data.r")

fatt1 <- data.frame(1:747)
fatt1$from_b3 <- fatt$from_b3
fatt1$from_tr <- fatt$from_tr
fatt1$from_pi <- fatt$from_pi

fatt1$cour <- fatt$cour
fatt1$only_est <- fatt$only_est
fatt1$oper_time <- fatt$oper_time
fatt1$oper_age <- fatt$oper_age
fatt1$ind <- fatt$ind
fatt1$no_load_docks <- fatt$no_load_docks
fatt1$home_del_pro <- fatt$home_del_pro
fatt1$off_street_park <- fatt$off_street_park
fatt1$C5610 <- fatt$C5610
fatt1$i_emp <- fatt$i_emp
fatt1$i_ea  <- fatt$i_ea
fatt1$C4772 <- fatt$C4772
fatt1$C4752 <- fatt$C4752
fatt1$C47713 <- fatt$C47713
fatt1$C4721 <- fatt$C4721
fatt1$owned_3w  <- fatt$owned_3w
fatt1$owned_bi  <- fatt$owned_bi
fatt1$owned_pi  <- fatt$owned_pi
fatt1$owned_tr  <- fatt$owned_tr
fatt1$leased_tr <- fatt$leased_tr
fatt1$cons      <- fatt$cons
fatt1$full_emp  <- fatt$full_emp
fatt1$part_emp  <- fatt$part_emp
fatt1$len_m_100 <- fatt$len_m_100

fatt1$lo <- fatt$lon
fatt1$la <- fatt$lat

fatt1 <- fatt1[-1]
fatt1 <- fatt1[complete.cases(fatt1),]

lo <- fatt1$lo
la <- fatt1$la

source("D:/R Codes from Lab System/kNN_SW.r")
xy=cbind(lo, la)
nn4 <- knearneigh(xy, k=8, longlat = TRUE)
nc.nn4 <- knn2nb(nn4)
list.nc.nn4 <- nb2listw(nc.nn4, style="W", zero.policy=T)
mat.nc.nn4 <- listw2mat(list.nc.nn4)

Y1 <- fatt1$from_b3
Y2 <- fatt1$from_pi
Y3 <- fatt1$from_tr

Tformula <- Y1 | Y2 | Y3 ~ ind+no_load_docks+C5610+i_emp+
  C4772+owned_3w+owned_bi|
  i_ea+cons+i_emp+C47713+C4772+C5610+C4721+
  owned_pi+C4772+home_del_pro|
  i_emp+no_load_docks + C4752
library(spsur)
LMs <- lmtestspsur(Tformula, data = fatt1, listw = list.nc.nn4)
LMs

spcsur.sim <- spsurml(formula=Tformula, data = fatt1, type = "sim", listw = list.nc.nn4)
summary(spcsur.sim)
AIC(spcsur.sim)

spcsur.slm <- spsurml(formula=Tformula, data = fatt1, type = "slm",
                        listw = list.nc.nn4)
summary(spcsur.slm)
AIC(spcsur.slm)

spsur::impactspsur(spcsur.slm, listw=list.nc.nn4)

bikesRMSEsur <- data.frame(spcsur.slm$fitted.values[1:712])
names(bikesRMSEsur) <- c("fit")
bikesRMSEsur$Y <- Y1
bikesRMSEsur <- bikesRMSEsur[abs(bikesRMSEsur$fit - bikesRMSEsur$Y)<10,]
rmse(bikesRMSEsur$fit, bikesRMSEsur$Y)

piRMSEsur <- data.frame(spcsur.slm$fitted.values[713:1424])
names(piRMSEsur) <- c("fit")
piRMSEsur$Y <- Y2
piRMSEsur <- piRMSEsur[abs(piRMSEsur$fit - piRMSEsur$Y)<10,]
rmse(piRMSEsur$fit, piRMSEsur$Y)

trRMSEsur <- data.frame(spcsur.slm$fitted.values[1425:2136])
names(trRMSEsur) <- c("fit")
trRMSEsur$Y <- Y3
trRMSEsur <- trRMSEsur[abs(trRMSEsur$fit - trRMSEsur$Y)<10,]
rmse(trRMSEsur$fit, trRMSEsur$Y)

rmse(spcsur.slm$fitted.values[1:712], Y1)
rmse(spcsur.slm$fitted.values[713:1424], Y2)
rmse(spcsur.slm$fitted.values[1425:2136], Y3)


bikesRMSE <- data.frame(OLS_fb3$fitted.values)
names(bikesRMSE) <- c("fit")
bikesRMSE$Y <- Y1
bikesRMSE <- bikesRMSE[abs(bikesRMSE$fit - bikesRMSE$Y)<10,]
rmse(bikesRMSE$fit, bikesRMSE$Y)

piRMSE <- data.frame(OLS_fpi$fitted.values)
names(piRMSE) <- c("fit")
piRMSE$Y <- Y2
piRMSE <- piRMSE[abs(piRMSE$fit - piRMSE$Y)<10,]
rmse(piRMSE$fit, piRMSE$Y)

trRMSE <- data.frame(OLS_ftr$fitted.values)
names(trRMSE) <- c("fit")
trRMSE$Y <- Y3
trRMSE <- trRMSE[abs(trRMSE$fit - trRMSE$Y)<10,]
rmse(trRMSE$fit, trRMSE$Y)


rmse(spcsur.sim$fitted.values[1:712], Y1)
rmse(spcsur.sim$fitted.values[713:1424], Y2)
rmse(spcsur.sim$fitted.values[1425:2136], Y3)


