## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----mod_calc,highlight=TRUE--------------------------------------------------
suppressMessages(library(tsmarch))
suppressMessages(library(xts))
data(globalindices)
Sys.setenv(TZ = "UTC")
train_set <- 1:1600
series <- 1:10
y <- as.xts(globalindices[, series])
train <- y[train_set,]
test <- escc_test(y, lags = 2, constant = TRUE)

## ----print_test,warning=FALSE,message=FALSE-----------------------------------
as_flextable(test, include.decision = TRUE, use.symbols = TRUE, footnote.reference = TRUE)

