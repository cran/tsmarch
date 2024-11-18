## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----mod_calc,highlight=TRUE--------------------------------------------------
suppressMessages(library(tsmarch))
suppressMessages(library(tstests))
suppressMessages(library(xts))
suppressMessages(library(shape))
suppressMessages(library(tsgarch))
suppressMessages(library(tsdistributions))
data(globalindices)
Sys.setenv(TZ = "UTC")
train_set <- 1:1600
test_set <- 1601:1698
series <- 1:5
y <- as.xts(globalindices[, series])
train <- y[train_set,]
mu <- colMeans(train)
train <- sweep(train, 2, mu, "-")
test <- y[test_set,]
test <- sweep(test, 2, mu, "-")
oldpar <- par(mfrow = c(1,1))

## ----gogarch_model------------------------------------------------------------
gogarch_mod <- gogarch_modelspec(train, distribution = "nig", model = "gjrgarch", components = 4) |> estimate()
summary(gogarch_mod)

## ----newsimpact,fig.width=6,fig.height=4--------------------------------------
gogarch_mod |> newsimpact(type = "correlation", pair = c(2,5), factor = c(1,3)) |> plot()

## ----var_calc,highlight=TRUE--------------------------------------------------
h <- 98
w <- rep(1/5, 5)
gogarch_filter_mod <- gogarch_mod
var_value <- rep(0, 98)
actual <- as.numeric(coredata(test) %*% w)

# first prediction without filtering update
var_value[1] <- predict(gogarch_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
for (i in 2:h) {
  gogarch_filter_mod <- tsfilter(gogarch_filter_mod, y = test[i - 1,])
  var_value[i]  <- predict(gogarch_filter_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
}

## -----------------------------------------------------------------------------
as_flextable(var_cp_test(actual, var_value, alpha = 0.1), include.decision = TRUE)

## ----comoments,highlight=TRUE-------------------------------------------------
V <- tscov(gogarch_mod)
S <- tscoskew(gogarch_mod, standardized = TRUE, folded = TRUE)
K <- tscokurt(gogarch_mod, standardized = TRUE, folded = TRUE)

## ----fig.width=6,fig.height=3-------------------------------------------------
p <- predict(gogarch_mod, h = 25, nsim = 1000, seed = 100)
K_p <- tscokurt(p, standardized = TRUE, folded = TRUE, distribution = TRUE, index = 1:25)
K_p <- t(K_p[1,1,1,1,,])
colnames(K_p) <- as.character(p$forc_dates)
class(K_p) <- "tsmodel.distribution"
L <- list(original_series = xts(K[1,1,1,1,], as.Date(gogarch_mod$spec$target$index)), distribution = K_p)
class(L) <- "tsmodel.predict"
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8)
plot(L, gradient_color = "orange", interval_color = "cadetblue", median_color = "black", median_type = 2, median_width = 1, 
     n_original = 100, main = "Kurtosis [AEX]", xlab = "", cex.main = 0.8)
par(oldpar)

## ----port_comoments,fig.width=6,fig.height=3----------------------------------
port_moments_estimate <- tsaggregate(gogarch_mod, weights = w)
port_moments_predict <- tsaggregate(p, weights = w, distribution = TRUE)
L <- list(original_series = port_moments_estimate$skewness, distribution = port_moments_predict$skewness)
class(L) <- "tsmodel.predict"
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8)
plot(L, gradient_color = "orange", interval_color = "cadetblue", median_color = "black", median_type = 2, median_width = 1, 
     n_original = 100, main = "Portfolio Skewness", xlab = "", cex.main = 0.8)
par(oldpar)

## ----convolution,fig.width=6,fig.height=5-------------------------------------
p <- predict(gogarch_mod, h = 98, nsim = 1000)
port_f_moments <- do.call(cbind, tsaggregate(p, weights = w, distribution = FALSE))
pconv <- tsconvolve(p, weights = w, fft_support = NULL, fft_step = 0.0001, fft_by = 0.00001, distribution = FALSE)
p_c_moments <- matrix(0, ncol = 4, nrow = 98)
for (i in 1:98) {
  df <- dfft(pconv, index = i)
  mu <- pconv$mu[i]
  f_2 <- function(x) (x - mu)^2 * df(x)
  f_3 <- function(x) (x - mu)^3 * df(x)
  f_4 <- function(x) (x - mu)^4 * df(x)
  sprt <- attr(pconv$y[[i]],"support")
  p_c_moments[i,2] <- sqrt(integrate(f_2, sprt[1], sprt[2], abs.tol = 1e-8, subdivisions = 500)$value)
  p_c_moments[i,3] <- integrate(f_3, sprt[1], sprt[2], abs.tol = 1e-8, subdivisions = 500)$value/p_c_moments[i,2]^3
  p_c_moments[i,4] <- integrate(f_4, sprt[1], sprt[2], abs.tol = 1e-8, subdivisions = 500)$value/p_c_moments[i,2]^4
}
par(mar = c(2,2,2,2), mfrow = c(3,1), pty = "m")
matplot(cbind(as.numeric(port_f_moments[,2]), p_c_moments[,2]), type = "l", lty = c(1,3), lwd = c(2, 2), col = c("grey","tomato1"), main = "Sigma", xaxt = "n")
grid()
matplot(cbind(as.numeric(port_f_moments[,3]), p_c_moments[,3]), type = "l", lty = c(1,3), lwd = c(2, 2), col = c("grey","tomato1"), main = "Skewness", xaxt = "n")
grid()
matplot(cbind(as.numeric(port_f_moments[,4]), p_c_moments[,4]), type = "l", lty = c(1,3), lwd = c(2, 2), col = c("grey","tomato1"), main = "Kurtosis")
grid()
par(oldpar)

## ----port_convolution,fig.width=6,fig.height=4--------------------------------
p <- predict(gogarch_mod, h = 98, nsim = 5000)
pconv <- tsconvolve(p, weights = w, fft_support = NULL, fft_step = 0.0001, fft_by = 0.00001, distribution = FALSE)
q_seq <- seq(0.025, 0.975, by = 0.025)
q_surface = matrix(NA, ncol = length(q_seq), nrow = 98)
for (i in 1:98) {
  q_surface[i,] <- qfft(pconv, index = i)(q_seq)
}
par(mar = c(1.8,1.8,1.1,1), pty = "m")
col_palette <- drapecol(q_surface, col = femmecol(100), NAcol = "white")
persp(x = 1:98, y = q_seq, z = q_surface,  col = col_palette, theta = 45, 
      phi = 15, expand = 0.5, ltheta = 25, shade = 0.25, 
      ticktype = "simple", xlab = "Time", ylab = "Quantile", 
      zlab = "VaR", cex.axis = 0.8, main = "Value at Risk Prediction Surface")
par(oldpar)

## ----pit----------------------------------------------------------------------
pit_value <- pit(pconv, actual)
as_flextable(shortfall_de_test(pit_value, alpha = 0.1), include.decision = TRUE)

## ----garch_dcc_model----------------------------------------------------------
garch_model <- lapply(1:5, function(i) {
  garch_modelspec(train[,i], model = "gjrgarch") |> estimate(keep_tmb = TRUE)
})
garch_model <- to_multi_estimate(garch_model)
names(garch_model) <- colnames(train)

## ----dcc_model, message=FALSE-------------------------------------------------
dcc_mod <- dcc_modelspec(garch_model, dynamics = "adcc", distribution = "mvt") |> estimate()
dcc_mod |> summary()

## ----dcc_newsimpact, fig.width=6,fig.height=4---------------------------------
newsimpact(dcc_mod, pair = c(1,2)) |> plot()

## ----dcc_var_calc,highlight=TRUE----------------------------------------------
h <- 98
w <- rep(1/5, 5)
dcc_filter_mod <- dcc_mod
var_value <- rep(0, 98)
actual <- as.numeric(coredata(test) %*% w)
# first prediction without filtering update
var_value[1] <- predict(dcc_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
for (i in 2:h) {
  dcc_filter_mod <- tsfilter(dcc_filter_mod, y = test[i - 1,])
  var_value[i]  <- predict(dcc_filter_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
}
as_flextable(var_cp_test(actual, var_value, alpha = 0.1), include.decision = TRUE)

## ----dcc_weighted, fig.width=6,fig.height=4-----------------------------------
p <- predict(dcc_mod, h = 98, nsim = 5000)
simulated_aggregate <- tsaggregate(p, weights = w, distribution = TRUE)
# we don't have any conditional mean dynamics but uncertainty around zero from the simulation
weighted_mu <- t(apply(p$mu, 1, rowMeans)) %*% w
H <- tscov(p, distribution = FALSE)
weighted_sigma <- sqrt(sapply(1:98, function(i) w %*% H[,,i] %*% w))
shape <- unname(coef(dcc_mod)["shape"])
simulated_var <- unname(apply(simulated_aggregate$mu, 2, quantile, 0.05))
analytic_var <- qstd(0.05, mu = weighted_mu, sigma = weighted_sigma, shape = shape)
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8, cex.main = 0.8)
plot(as.Date(p$forc_dates), simulated_var, type = "l", ylab = "", xlab = "", main = "Value at Risk [5%]", ylim = c(-0.039, -0.033))
lines(as.Date(p$forc_dates), analytic_var, col = 2, lty = 2)
legend("topright", c("Simulated","Analytic"), col = 1:2, lty = 1:2, bty = "n")
par(oldpar)

## ----garch_copula_model-------------------------------------------------------
distributions <- c(rep("jsu",4), rep("sstd",1))
garch_model <- lapply(1:5, function(i) {
  garch_modelspec(train[,i], model = "gjrgarch", distribution = distributions[i]) |> estimate(keep_tmb = TRUE)
})
garch_model <- to_multi_estimate(garch_model)
names(garch_model) <- colnames(train)

## ----cgarch_model, message=FALSE----------------------------------------------
cgarch_mod <- cgarch_modelspec(garch_model, dynamics = "adcc", 
                               transformation = "parametric", 
                               copula = "mvt") |> estimate()
cgarch_mod |> summary()

## ----cgarch_newsimpact, fig.width=6,fig.height=4------------------------------
newsimpact(cgarch_mod, pair = c(1,2)) |> plot()

## ----cgarch_var_calc,highlight=TRUE-------------------------------------------
h <- 98
w <- rep(1/5, 5)
cgarch_filter_mod <- cgarch_mod
var_value <- rep(0, 98)
actual <- as.numeric(coredata(test) %*% w)
# first prediction without filtering update
var_value[1] <- predict(cgarch_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
for (i in 2:h) {
  cgarch_filter_mod <- tsfilter(cgarch_filter_mod, y = test[i - 1,])
  var_value[i]  <- predict(cgarch_filter_mod, h = 1, nsim = 5000, seed = 100) |> value_at_risk(weights = w, alpha = 0.1)
}
as_flextable(var_cp_test(actual, var_value, alpha = 0.1), include.decision = TRUE)

## ----cgarch_weighted, fig.width=6,fig.height=4--------------------------------
p <- predict(cgarch_mod, h = 98, nsim = 5000)
simulated_aggregate <- tsaggregate(p, weights = w, distribution = TRUE)
simulated_var <- unname(apply(simulated_aggregate$mu, 2, quantile, 0.05))
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8, cex.main = 0.8)
plot(as.Date(p$forc_dates), simulated_var, type = "l", ylab = "", xlab = "", main = "Value at Risk [5%]")
par(oldpar)

## -----------------------------------------------------------------------------
arima_model <- lapply(1:5, function(i){
  arima(train[,i], order = c(6,0,0), method = "ML")
})
.residuals <- do.call(cbind, lapply(arima_model, function(x) as.numeric(residuals(x))))
colnames(.residuals) <- colnames(train)
.residuals <- xts(.residuals, index(train))
.fitted <- train - .residuals
.predicted <- do.call(cbind, lapply(1:5, function(i){
  as.numeric(predict(arima_model[[i]], n.ahead = 25)$pred)
}))
colnames(.predicted) <- colnames(train)

## -----------------------------------------------------------------------------
dcc_mod_mean <- dcc_modelspec(garch_model, dynamics = "adcc", distribution = "mvt", cond_mean = .fitted) |> estimate()
all.equal(fitted(dcc_mod_mean), .fitted)

## -----------------------------------------------------------------------------
p <- predict(dcc_mod_mean, h = 25, cond_mean = .predicted, nsim = 5000, seed = 100)
simulated_mean <- as.matrix(t(apply(p$mu, 1, rowMeans)))
colnames(simulated_mean) <- colnames(train)
all.equal(simulated_mean, .predicted)

## ----fig.width=6,fig.height=4-------------------------------------------------
res <- p$mu
arima_pred <- lapply(1:5, function(i){
  # we eliminate the mean prediction from the simulated predictive distribution
  # to obtain the zero mean innovations
  res_i <- scale(t(res[,i,]), scale = FALSE, center = TRUE)
  sim_p <- do.call(rbind, lapply(1:5000, function(j) {
    arima.sim(model = list(ar = coef(arima_model[[i]])[1:6]), n.start = 20, n = 25, innov = res_i[j,], start.innov = as.numeric(tail(.residuals[,i],20))) |> as.numeric() + coef(arima_model[[i]])[7]
  }))
  return(sim_p)
})
arima_pred <- array(unlist(arima_pred), dim = c(5000, 25, 5))
arima_pred <- aperm(arima_pred, c(2, 3, 1))

simulated_mean <- as.matrix(t(apply(arima_pred, 1, rowMeans)))
colnames(simulated_mean) <- colnames(train)
par(mfrow = c(3,2), mar = c(2,2,2,2))
for (i in 1:5) {
  matplot(cbind(simulated_mean[,i], .predicted[,i]), type = "l", lty = c(1,3), lwd = c(2, 2), col = c("grey","tomato1"), ylab = "", xaxt = "n")
  grid()
}
par(oldpar)

## ----fig.width=6,fig.height=4-------------------------------------------------
i <- 1
sim_1a <- t(p$mu[,i,])
sim_1b <- t(arima_pred[,i,])
colnames(sim_1a) <- colnames(sim_1b) <- as.character(p$forc_dates)
class(sim_1a) <- class(sim_1b) <- "tsmodel.distribution"
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8)
plot(sim_1a, gradient_color = "whitesmoke", interval_color = "orange", median_color = "orange")
plot(sim_1b, add = TRUE, gradient_color = "whitesmoke", interval_color = "steelblue", median_color = "steelblue", median_type = 2)
par(oldpar)

## ----fig.width=6,fig.height=4-------------------------------------------------
j <- 2
sim_2a <- t(p$mu[,j,])
sim_2b <- t(arima_pred[,j,])
colnames(sim_2a) <- colnames(sim_2b) <- as.character(p$forc_dates)
class(sim_2a) <- class(sim_2b) <- "tsmodel.distribution"
C_a <- sapply(1:25, function(i) cor(sim_1a[,i], sim_2a[,i]))
C_b <- sapply(1:25, function(i) cor(sim_1b[,i], sim_2b[,i]))
par(mar = c(2,2,1.1,1), pty = "m", cex.axis = 0.8, cex.main = 0.8)
matplot(cbind(C_a, C_b), type = "l", lty = c(1,3), lwd = c(2, 2), col = c("grey","tomato1"), ylab = "", main = "Pairwise Correlation")
grid()
par(oldpar)

