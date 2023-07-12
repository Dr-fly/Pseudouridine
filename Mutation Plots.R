#####Mutation Plots#####
#Mutation Violin Plot
ggbetweenstats(
  data = MutationBurden,
  x = Type,
  y = mutation_total,
  xlab = "Risk",
  ylab = "Mutation Burden",
  type = "nonparametric",
  centrality.plotting = T,
  plot.type = "violin",
  results.subtitle = FALSE,
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.6),
    alpha = 0.6,
    size = 1.5,
    stroke = 0
  )
)
#Survival Plots for risk and genes
Survplot <- function(survplotdata, classNM,name) {
  require(survminer)
  require(survival)
  if (classNM == "risk") {
    
  }
  else if (classNM %in% colnames(survplotdata)) {
    survplotdata[, classNM] <-
      ifelse(survplotdata[, classNM] >= median(survplotdata[, classNM]),
             'High',
             'Low')
  }
  pdf(
    file = paste0("Survival-",name,"-", classNM, ".pdf"),
    width = 5,
    height = 5,
    onefile = F
  )
  fml <- as.formula(paste0("Surv(overall_survival, deceased) ~", classNM))
  fit <- surv_fit(fml, data = survplotdata)
  print(
    ggsurvplot(
      fit,
      data = survplotdata,
      risk.table = TRUE,
      pval = TRUE,
      break.time.by = 500,
      ggtheme = theme_minimal(),
      risk.table.y.text.col = TRUE,
      risk.table.y.text = FALSE,
      palette = c(mypal[1], mypal[4])
    )
  )
  dev.off()
}
#Immune Score Plotting
library(ggstatplot)
ggbetweenstats(
  data = mf1,
  x = Risk,
  y = StromalScore,
  subtitle = "Stromal Score",
  xlab = "",
  ylab = "Stromal Score",
  type = "nonparametric",
  centrality.plotting = T,
  bf.message = F,
  pairwise.comparisons = F,
  plot.type = "violin",
  results.subtitle = FALSE,
  point.args = list(
    position = ggplot2::position_jitterdodge(dodge.width = 0.6),
    alpha = 0.8,
    size = 1.5,
    stroke = 0
  )
) + theme(plot.subtitle = element_text(size = 15,
                                       face = "bold",
                                       color = "Black"))
ggscatterstats(
  data = mf1,
  ## data frame from which variables are taken
  x = risk,
  ## predictor/independent variable
  y = StromalScore,
  ## dependent variable
  xlab = "Risk Score",
  ## label for the x-axis
  ylab = "Stromal Score",
  ## label for the y-axis
  point.label.args = list(
    alpha = 0.7,
    size = 4,
    color = "grey50"
  ),
  xfill = mycolors[1],
  yfill = mycolors[2],
  ## fill for marginals on the x-axis
  point.args = list(
    size = 2.2,
    alpha = 0.8,
    stroke = 0,
    na.rm = TRUE,
    color = mycolors[4]
  ),
  smooth.line.args = list(
    size = 1.5,
    color = "blue",
    method = "lm",
    formula = y ~ x
  ),
  bf.message = F
) + scale_y_log10()