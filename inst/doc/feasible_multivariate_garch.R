## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "##",
  R.options = list(width = 60)
)

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(..., dynamics = 'constant')

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(..., dynamics = 'dcc')

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(...) |> estimate(...)

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(...) |> estimate(...) |> vcov(...)

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(distribution = 'mvt')

## ----eval=FALSE-------------------------------------------
#  cgarch_modelspec(...)

## ----eval=FALSE-------------------------------------------
#  cgarch_modelspec(..., dynamics = 'adcc')

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(...) |> estimate(...) |> predict(...)
#  dcc_modelspec(...) |> estimate(...) |> simulate(...)

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(...) |> estimate(...) |> tsaggregate(...)

## ----eval=FALSE-------------------------------------------
#  dcc_modelspec(...) |> estimate(...) |> newsimpact(...) |> plot(...)

## ----eval=FALSE-------------------------------------------
#  gogarch_modelspec(...)

## ----eval=FALSE-------------------------------------------
#  gogarch_modelspec(...) |> estimate(...)

## ----eval=FALSE-------------------------------------------
#  tscov(...)
#  tscor(...)
#  tscoskew(...)
#  tscoskurt(...)

## ----eval=FALSE-------------------------------------------
#  gogarch_modelspec(...) |> estimate(...) |> predict(...) |> tsaggregate(...)

## ----eval=FALSE-------------------------------------------
#  gogarch_modelspec(...) |> estimate(...) |> predict(...) |> tsconvolve(...) |> dfft(...)

## ----eval=FALSE-------------------------------------------
#  gogarch_modelspec(...) |> estimate(...) |> newsimpact(...) |> plot(...)

