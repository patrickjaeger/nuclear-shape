# https://stackoverflow.com/questions/31851936/exponential-curve-fitting-in-r

t <- 1:100      # these are your time points
a <- 10         # assume the size at t = 0 is 10
r <- -0.1        # assume a growth constant
y <- a*exp(r*t) # generate some y observations from our exponential model

# visualise
par(mfrow = c(1, 2))
plot(t, y)      # on the original scale
plot(t, log(y)) # taking the log(y)



set.seed(12) # for reproducible results

# errors constant across time - additive
y_add <-  a*exp(r*t) + rnorm(length(t), sd = 5000) # or: rnorm(length(t), mean = a*exp(r*t), sd = 5000)

# errors grow as y grows - multiplicative (constant on the log-scale)
y_mult <- a*exp(r*t + rnorm(length(t), sd = 1))  # or: rlnorm(length(t), mean = log(a) + r*t, sd = 1)

# visualise
par(mfrow = c(1, 2))
plot(t, y_add, main = "additive error")
lines(t, a*exp(t*r), col = "red") 
plot(t, y_mult, main = "multiplicative error")
lines(t, a*exp(t*r), col = "red")



add_nls <- nls(y_add ~ a*exp(r*t), start = list(a = 10, r = -0.2))
coef(add_nls)



plot(t, resid(add_nls))
abline(h = 0, lty = 2)



mult_lm <- lm(log(y_mult) ~ t)
coef(mult_lm)



mult_nls <- nls(y_mult ~ a*exp(r*t), start = list(a = 20, r = -0.2))
coef(mult_nls)



# get the model's coefficients
lm_coef <- coef(mult_lm)
nls_coef <- coef(mult_nls)

# make the plot
dev.off()
plot(t, y_mult)
lines(t, a*exp(r*t), col = "brown", lwd = 5)
lines(t, exp(lm_coef[1])*exp(lm_coef[2]*t), col = "dodgerblue", lwd = 2)
lines(t, nls_coef[1]*exp(nls_coef[2]*t), col = "orange2", lwd = 2)
legend("topleft", col = c("brown", "dodgerblue", "orange2"), 
       legend = c("known model", "lm fit", "nls fit"), lwd = 3)



plot(t, resid(mult_nls))
abline(h = 0, lty = 2)

plot(t, resid(mult_lm))
abline(h = 0, lty = 2)
