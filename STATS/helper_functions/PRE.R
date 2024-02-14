# Intellectual property of Pr. Dominique Muller; UGA LIP/PC2S, Grenoble, France

PRE <- function(model) {
  library(car)
  A <- Anova(model, type = 3)
  A
  R <- A[1]
  R
  df <- A[2]
  df
  dim <- dim(R)
  dim
  d <- dim[1]
  d
  SSE <- R[d, 1]
  SSE
  d1 <- d - 1
  SST <- sum(R[2:d, ])
  SST
  ETASQ <- R / SST
  ETASQ
  ETASQP <- R / (SSE + R)
  ETASQP
  ANOVA <- cbind(A, ETASQP, ETASQ)
  ANOVA <- round(ANOVA, digits = 3)
  colnames(ANOVA) <- c("SumSq", "df", "F", "p", "h2p", "h2")
  ANOVA
  ANOVA[1, c(5, 6)] <- c(" ", " ")
  ANOVA[d, c(3, 4, 5)] <- c(" ", " ", " ")
  ANOVA
  print(ANOVA)
}
