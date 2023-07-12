#####Nomogram#####
#Nomogram
regplot(
  mod.sur,
  points = T,
  plots = c("violin", "spikes"),
  failtime = c(365 * 3, 730, 365),
  prfail = F,
  center = T,
  title = "",
  subticks = T,
  cexscales = 0.7,
  cexcats = 0.7
  
)
p <- recordPlot()
pdf(paste0("TCGA_NOMOgram.pdf"),
    width = 8,
    height = 5)
print(p)
dev.off()
rm(p)
#Calibration 
pdf(paste0("TCGA_Calibration.pdf"),
    width = 5,
    height = 5)
f <-
  cph(
    as.formula(fml),
    x = T,
    y = T,
    surv = T,
    data = datamod,
    time.inc = 365
  )
cal <-
  calibrate(
    f,
    cmethod = "KM",
    method = "boot",
    u = 365,
    m = (nrow(datamod) / 4),
    B = 1000
  )
plot(
  cal,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "Nomogram-predicted futime (%)",
  ylab = "Observed futime (%)",
  lwd = 1.5,
  col = "#7CC767",
  sub = F
)
f <-
  cph(
    as.formula(fml),
    x = T,
    y = T,
    surv = T,
    data = datamod,
    time.inc = 730
  )
cal <-
  calibrate(
    f,
    cmethod = "KM",
    method = "boot",
    u = 730,
    m = (nrow(datamod) / 4),
    B = 1000
  )
plot(
  cal,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lwd = 1.5,
  col = "#0066FF",
  sub = F,
  add = T
)
f <-
  cph(
    as.formula(fml),
    x = T,
    y = T,
    surv = T,
    data = datamod,
    time.inc = 365 * 3
  )
cal <-
  calibrate(
    f,
    cmethod = "KM",
    method = "boot",
    u = 365 * 3,
    m = (nrow(datamod) / 4),
    B = 1000
  )
plot(
  cal,
  xlim = c(0, 1),
  ylim = c(0, 1),
  xlab = "",
  ylab = "",
  lwd = 1.5,
  col = "#FF0000",
  sub = F,
  add = T
)
legend(
  'bottomright',
  c('1-year', '2-year', '3-year'),
  col = c("#7CC767", "#0066FF", "#FF0000"),
  lwd = 1.5,
  bty = 'n'
)
dev.off()