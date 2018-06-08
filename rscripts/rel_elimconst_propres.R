# Quick exploration into the relationship between elimination rate constant and
# proportional residual
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  res <- c(0, 0.108, 0.214, 0.667, 0.866, 0.285, 0.135, 0.374, 0.678)
  K <- c(0, 0.57, 1.23, 5.7, 10.43, 1.71, 0.724, 2.41, 5.87)
  MMfun <- function(x, xmax, km) {
    xmax*x/(km + x)
  }
  expfun <- function(x, m) {
    1 - exp(-m*x)
  }
  Kseq <- seq(0, 12, by = 0.1)
  # reshat <- MMfun(Kseq, 1, 3.7)
  reshat <- expfun(Kseq, 1/5.2)
  df <- data.frame(K, res)
  dfhat <- data.frame(Kseq, reshat)
  
  plot(df)
  
  p <- NULL
  p <- ggplot()
  p <- p + geom_point(aes(x = K, y = res), data = df)
  # p <- p + geom_abline(slope = 0.19, linetype = "dashed")
  p <- p + geom_line(aes(x = Kseq, y = reshat), data = dfhat, linetype = "dashed")
  p
