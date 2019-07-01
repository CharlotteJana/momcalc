# raw moments of a one dimensional gamma distribution
symbolicMoments(distribution = "gamma", missingOrders = as.matrix(1:3, ncol = 1),
                mean = "μ", var = "σ")

# raw moments of a one dimensional lognormal distribution 
symbolicMoments(distribution = "lognormal", missingOrders = as.matrix(1:2, ncol = 1),
                mean = 2, var = 1, simplify = FALSE)

# evaluate the result
symbolicMoments(distribution = "lognormal", missingOrders = as.matrix(1:2, ncol = 1),
                mean = 2, var = 1, simplify = TRUE)

#### central moments of a four dimensional normal distribution ####

missingOrders <- matrix(c(4, 0, 0, 0,
                          3, 1, 0, 0,
                          2, 2, 0, 0,
                          2, 1, 1, 0,
                          1, 1, 1, 1), ncol = 4, byrow = TRUE)

cov <-  matrix(c("σ11", "σ12", "σ13", "σ14", 
                 "σ12", "σ22", "σ23", "σ24", 
                 "σ13", "σ23", "σ33", "σ34", 
                 "σ14", "σ24", "σ34", "σ44"), ncol = 4, byrow = TRUE)

symbolicMoments("normal", missingOrders, mean = "μ", cov = cov)