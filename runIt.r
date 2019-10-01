source("brook90Functions.r")
source("brook90Variables.r")

dataPath <- "./Input_data"
source("brook90ReadData.r")

source("brook90Run.r")

## ----chunkplot, results='hide'-------------------------------------------
maxSF   <- max(c(timeseries_flow[1:NDAYS], timeseries_mesfld,timeseries_evp), na.rm = T)
maxPR   <- max(timeseries_prec, na.rm = T)
par(mar = c(4, 4, 3, 4) + 0.2)
plot(1:NDAYS, timeseries_evp,
     type = 'l', col = "darkseagreen2",
     ylim = c(0, 14),
     xaxs = "i", yaxs = "i",
     xlab = "", ylab = "",
     lwd=2)
grid (10,10, lty = 6, col = "grey89")
lines(1:NDAYS, timeseries_mesfld, col = "steelblue",lwd=2)
lines(1:NDAYS, timeseries_flow[1:NDAYS], col = "indianred3",lwd=2)

par(new = TRUE)
plot(x = 1:NDAYS, y = rep(0, length(timeseries_prec)),
     type = "n", ylim = c(5 * maxPR, 0),
     xaxs = "i", yaxs = "i",
     axes = FALSE, xlab = "", ylab = "")
segments(x0 = 1:NDAYS, y0 = rep(0, length(timeseries_prec)),
         x1 = 1:NDAYS, y1 = timeseries_prec,
         lend = 2, lwd =1, col="deepskyblue3")

yrAxis  <- seq(0, ceiling(maxPR), length.out = 5)
axis(4, at = yrAxis, labels = paste0(yrAxis))
#       mtext(y = yrAxis, par(usr)[1], labels = yrAxis)
mtext("Precipitation [mm/d]", side = 4, line = 2, adj = 1)
mtext("Day of the Year 1999 [d]", side = 1, line = 2, adj = 1)
mtext("Streamflow or Evapotranspiration [mm/d]", side = 2, line = 2, adj = 1)

legend("topright",
       inset=c(0.05,0.05),
       xpd=TRUE,
       legend=c("Simulated Streamflow [mm/d]","Observed Streamflow [mm/d]","Observed Precipitation [mm/d]","Simulated Evapotranspiration [mm/d]"),
       col=c("indianred3","steelblue","deepskyblue3","darkseagreen2"),
       lty=c("solid","solid","solid","solid"),
       cex=0.9,
       lwd=c(2,2,1,2),
       y.intersp = 0.8,
       bty="n")
