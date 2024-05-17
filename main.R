if (!require("rugarch")) install.packages("rugarch", dependencies = TRUE)
library(rugarch)
if (!require("forecast")) install.packages("forecast", dependencies = TRUE)
library(forecast)
if (!require("tseries")) install.packages("tseries", dependencies = TRUE)
library(tseries)


data = read.csv("^GDAXI.csv")
data = na.omit(data)

data$log_rets = c(NA, diff(log(data$GDAXI.Close)))
data = na.omit(data)
data$Index = as.Date(data$Index)

adf.test(data$log_rets)

# Examine the time series of returns
acf(data$log_rets, main = "ACF of Log Returns")
pacf(data$log_rets, main = "PACF of Log Returns")
# no concrete pattern, we can see significant 4 and 7 lags for both of them

# Define the range of p and q values to consider
p_values <- c(0, 1, 2, 4, 7)
q_values <- c(0, 1, 2, 4, 7)

# Initialize matrices to store AIC and BIC values
aic_matrix <- matrix(NA, nrow = length(p_values), ncol = length(q_values))
bic_matrix <- matrix(NA, nrow = length(p_values), ncol = length(q_values))

# Iterate over all combinations of p and q values
for (i in seq_along(p_values)) {
  for (j in seq_along(q_values)) {
    p <- p_values[i]
    q <- q_values[j]
    
    # Fit the ARMA(p, q) model
    arma_model <- Arima(data$log_rets, order = c(p, 0, q))
    
    # Store AIC and BIC values in the matrices
    aic_matrix[i, j] <- AIC(arma_model)
    bic_matrix[i, j] <- BIC(arma_model)
  }
}

# Convert matrices to data frames for better printing
aic_df <- as.data.frame(aic_matrix)
colnames(aic_df) <- q_values
rownames(aic_df) <- p_values

bic_df <- as.data.frame(bic_matrix)
colnames(bic_df) <- q_values
rownames(bic_df) <- p_values

# Print the AIC and BIC matrices
cat("AIC values:\n")
print(aic_df)

cat("\nBIC values:\n")
print(bic_df)

# Find the minimum AIC and BIC values and corresponding p and q
min_aic <- min(aic_matrix)
min_bic <- min(bic_matrix)

aic_indices <- which(aic_matrix == min_aic, arr.ind = TRUE)
bic_indices <- which(bic_matrix == min_bic, arr.ind = TRUE)

cat("\nMinimum AIC:", min_aic, "at (p, q) =", p_values[aic_indices[1]], q_values[aic_indices[2]], "\n")
cat("Minimum BIC:", min_bic, "at (p, q) =", p_values[bic_indices[1]], q_values[bic_indices[2]], "\n")

# different results -> combined
comb_matrix = aic_matrix + bic_matrix
comb_df = data.frame(comb_matrix)
colnames(comb_df) <- q_values
rownames(comb_df) <- p_values
min_score = min(comb_matrix)
min_index = which(comb_matrix == min_score, arr.ind = TRUE)
cat("\nMinimum combined score:", min_score, "at (p, q) =", p_values[min_index[1]], q_values[min_index[2]], "\n")

arma_order = c(0,0)


# Function to fit models with different garchOrder
fit_models <- function(data, garchOrder, arma_order = c(0,0)) {
  spec_garch <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = garchOrder),
    mean.model = list(armaOrder = arma_order, include.mean = TRUE),
    distribution.model = "norm"
  )
  
  spec_tarch <- ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = garchOrder),
    mean.model = list(armaOrder = arma_order, include.mean = TRUE),
    distribution.model = "norm"
  )
  
  fit_garch <- ugarchfit(spec = spec_garch, data = data$log_rets)
  fit_tarch <- ugarchfit(spec = spec_tarch, data = data$log_rets)
  
  return(list(garch = fit_garch, tarch = fit_tarch))
}


# Try different garchOrder and select the best model
garchOrders <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))
model_fits <- lapply(garchOrders, function(order) fit_models(data, order))

# Extract AIC, BIC, and log-likelihood values for GARCH models
garch_aic <- sapply(model_fits, function(fit) infocriteria(fit$garch)[1])
garch_bic <- sapply(model_fits, function(fit) infocriteria(fit$garch)[2])
garch_loglik <- sapply(model_fits, function(fit) fit$garch@fit$LLH)

# Extract AIC, BIC, and log-likelihood values for TARCH models
tarch_aic <- sapply(model_fits, function(fit) infocriteria(fit$tarch)[1])
tarch_bic <- sapply(model_fits, function(fit) infocriteria(fit$tarch)[2])
tarch_loglik <- sapply(model_fits, function(fit) fit$tarch@fit$LLH)

# Create the model comparison table
model_comparison <- data.frame(
  garchOrder = sapply(garchOrders, paste, collapse = ","),
  GARCH_AIC = garch_aic,
  GARCH_BIC = garch_bic,
  GARCH_LogLik = garch_loglik,
  TARCH_AIC = tarch_aic,
  TARCH_BIC = tarch_bic,
  TARCH_LogLik = tarch_loglik
)
print(model_comparison)

# Select the best model based on the chosen criterion (e.g., AIC)
best_garch_order <- garchOrders[[which.min(model_comparison$GARCH_AIC)]]
best_tarch_order <- garchOrders[[which.min(model_comparison$TARCH_AIC)]]

# Fit the best models
best_garch_fit <- ugarchfit(spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = best_garch_order),
                                              mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                              distribution.model = "norm"),
                            data = data$log_rets)

best_tarch_fit <- ugarchfit(spec = ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = best_tarch_order),
                                              mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                              distribution.model = "norm"),
                            data = data$log_rets)

if (!require("zoo")) install.packages("zoo", dependencies = TRUE)
library(zoo)
window_size <- 30
rolling_volatility <- rollapply(data$log_rets, width = window_size, FUN = sd, by.column = TRUE, align = "right", fill = NA)
data$sigma = rolling_volatility
data$sigma_garch = as.vector(sigma(best_garch_fit))
data$sigma_tarch = as.vector(sigma(best_tarch_fit))


if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
# Close Price
ggplot(data, aes(Index, GDAXI.Close)) +
  geom_line() +
  labs(
    x = "Time Index",
    y = "DAX Price")
# Sigma + returns
ggplot(data[31:dim(data)[1],], aes(Index, sigma)) + 
  geom_line() +
  geom_line(aes(y = sigma_garch), color = "red") +
  geom_line(aes(y = sigma_tarch), color = "blue")
# Sigma + realized sigma
ggplot(data[31:nrow(data), ], aes(x = Index)) + 
  geom_line(aes(y = log_rets, color = "Log Returns")) +
  geom_line(aes(y = sigma_garch, color = "GARCH")) +
  geom_line(aes(y = sigma_tarch, color = "TARCH")) +
  labs(
       x = "Time Index",
       y = "Volatility") +
  scale_color_manual(values = c("Log Returns" = "black", "GARCH" = "red", "TARCH" = "blue"),
                     labels = c("Log Returns", "GARCH", "TARCH"))

