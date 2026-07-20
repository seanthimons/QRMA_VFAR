iter <- c(16,100,500,1000)
elap <- c(2.2,20.46,457.28, 1703)
plot(iter, elap)

model <- lm(log(elap)~log(iter))

summary(model)

#ln(elap) = (-3.9792) + 1.6238(iter)
  


(exp(-3.9792)*(10000)^1.6238)/60

# 16.24624 hrs

iter <- c(16,
          104,
          504,
          1000,
          5000)
elap <- c(2.47,
          3.35,
          19.42,
          64.14,
          1436.3)
plot(iter, elap)

model <- lm(log(elap)~log(iter))

summary(model)

(exp(-3.1076)*(10000)^1.1027)/60

#~19.1 min, 98.08% decrease
