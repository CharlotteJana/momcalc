#v1 Warum funktioniert der Test für lognormal, n = 2 nicht?

context("symbolic moments")

test_that("errors are thrown when wrong input", {
  
  cov <-  matrix(c(1, -1, 3, 2, 4, 5, 3, 5, 6), nrow = 3) # not symmetric
  expect_error(symbolicMoments(distribution = 'normal', 
                               cov = cov, 
                               missingOrders = 1:3))
  expect_identical(symbolicMoments(distribution = 'zero', 
                                   cov = cov, 
                                   missingOrders = 1:3), 
                   0) # cov is not needed for 'zero' and therfore not validated
  
  cov[2] <- 2 # now it is symmetric
  expect_error(symbolicMoments(distribution = 'normal', 
                               cov = cov, 
                               missingOrders = 1:2)) # wrong dimensions 
  expect_error(symbolicMoments(distribution = 'gamma', 
                               cov = cov, 
                               missingOrders = 1:3, 
                               mean = 1:2)) # wrong dimensions of mean
  expect_error(symbolicMoments(distribution = 'bla', 
                               cov = cov, 
                               missingOrders = 1:3)) # wrong distribution
})

test_that("distribution = 'zero' works", {
  moment <- symbolicMoments(distribution = 'zero', 
                            missingOrders = matrix(1:6, nrow = 2))
  expect_equal(moment, c(0, 0))
})

test_that("distribution = 'NA' works", {
  moment <- symbolicMoments(distribution = 'NA', 
                            missingOrders = matrix(1:3, nrow = 1))
  expect_equal(moment, NA)
})

test_that("distribution = 'normal' works", {
  
  cov <-  matrix(c("a", "b", "c", "d", 
                   "b", "e", "f", "h", 
                   "c", "f", "g", "i", 
                   "d", "h", "i", "j"), ncol = 4, byrow = TRUE)

  # order is odd -> moment should be zero
  moment <- symbolicMoments(distribution = 'normal', 
                            missingOrders = c(1,2,6,4), 
                            cov = cov)
  expect_equal(moment, 0)
  
  # n = 4 and order is even (the example can be found on Wikipedia)
  missingOrders <- matrix(c(4, 0, 0, 0,
                            3, 1, 0, 0,
                            2, 2, 0, 0,
                            2, 1, 1, 0,
                            1, 1, 1, 1), ncol = 4, byrow = TRUE)
  
  # with simplify = FALSE
  mom1 <- symbolicMoments(distribution = 'normal', 
                          missingOrders = missingOrders, 
                          cov = cov, simplify = FALSE)
  mom2 <- list(
    quote(3 * "a"^2),
    quote(3 * ("a"^1 * "b"^1)),
    quote(2 * "b"^2 + 1 * ("a"^1 * "e"^1)),
    quote(2 * ("b"^1 * "c"^1) + 1 * ("a"^1 * "f"^1)),
    quote(1 * ("d"^1 * "f"^1) + 1 * ("c"^1 * "h"^1) + 1 * ("b"^1 * "i"^1))
  )
  expect_equal(mom1, mom2)
  
  # with simplify = TRUE
  mom1 <- symbolicMoments(distribution = 'normal', 
                          missingOrders = missingOrders, 
                          cov = cov, simplify = TRUE)
  mom2 <- c(quote(3 * a^2),
            quote(3 * (a * b)),
            quote(2 * b^2 + a * e),
            quote(2 * (b * c) + a * f),
            quote(b * i + c * h + d * f))
  expect_equal(mom1, mom2)
  
  # n = 1 and different orders
  mom1 <- symbolicMoments(distribution = 'normal', var = 4, simplify = FALSE,
                          missingOrders = as.matrix(1:8, ncol = 1))
  mom2 <- symbolicMoments(distribution = 'normal', var = 4, simplify = TRUE,
                          missingOrders = as.matrix(1:8, ncol = 1))
  mom1 <- sapply(mom1, eval)
  mom3 <- actuar::mnorm(1:8, mean = 0, sd = 2)
  expect_equal(mom1, mom3)
  expect_equal(mom2, mom3)
})

test_that("distribution = 'lognormal' works for n = 1", {

  order <- 5
  meanlog <- 1
  sdlog <- 0.5
  mom1 <- actuar::mlnorm(order, meanlog = meanlog, sdlog = sdlog)
  
  mean <- exp(meanlog+0.5*sdlog^2)
  var <- exp(2*meanlog+sdlog^2)*(exp(sdlog^2)-1)
  mom2 <- symbolicMoments(distribution = 'lognormal', 
                          missingOrders = order, 
                          var = var, mean = mean, simplify = FALSE)
  mom2 <- sapply(mom2, eval)
  mom3 <- symbolicMoments(distribution = 'lognormal', 
                          missingOrders = order, 
                          var = var, mean = mean, simplify = TRUE)
  
  expect_equal(mom1, mom2)
  expect_equal(mom1, mom3, tolerance = 1e-04)
})

test_that("distribution = 'lognormal' works for n = 2", {
  
  cov <-  matrix(c(6, 2, 2, 4), nrow = 2)
  mean <- c(2, 4)
  order <- c(2, 2)
  
  mom1 <- 4*log(cov[1,2] + mean[1]*mean[2]) +
          log(cov[1,1] + mean[1]^2) - 4*log(mean[1]) +
          log(cov[2,2] + mean[2]^2) - 4*log(mean[2])
  mom1 <- exp(mom1)
  
  mom2 <- symbolicMoments(distribution = 'lognormal', missingOrders = order, 
                          cov = cov, mean = mean, simplify = FALSE)
  expect_identical(mom1, eval(mom2[[1]]))
  
})

test_that("distribution = 'gamma' works for n = 1", {
  
  order <- c(3, 8)
  alpha <- 1.5
  beta <- 2
  
  mom1 <- actuar::mgamma(order, shape = alpha, scale = beta)
  mom2 <- symbolicMoments(distribution = "gamma", 
                          missingOrders = as.matrix(order, ncol = 1),
                          var = alpha*beta^2, 
                          mean = alpha*beta, 
                          simplify = TRUE)
  
  expect_equal(mom1, mom2)
})

test_that("distribution = 'gamma' works for n = 2", {
  # this example is based on [Lak+15], Appendix B
  
  cov <-  matrix(c(6, 2, 2, 4), nrow = 2)
  mean <- c(3, 4)
  order <- c(2, 1)
  
  # in this case we have beta = c(2, 1) and
  # A₁₁ = 1/2, A₁₂ = A₂₁ = 1, A₂₂ = 3
  
  # the following summands should appear in the result:
  # (A₁₁)₂*A₁₂ = (1/2 + 1)*(1/2)*1 = 0.75
  # 2*A₁₁*(A₁₂)₂ = 2*(1/2)*(1 + 1)*(1) = 2
  # (A₁₂)₃ = (1 + 2)*(1 + 1)*(1) = 6
  # (A₁₁)₂*A₂₂ = (1/2 + 1)*(1/2)*3 = 2.25
  # 2*A₁₁*A₁₂*A₂₂ = 2*(1/2)*1*3 = 3
  # (A₁₂)₂*A₂₂ = 3*(1 + 1)*(1) = 6
  # the sum is 20
  # multplied by beta[1]²*beta[2] = 4 leads to 80
  
  moment <- symbolicMoments(distribution = 'gamma', missingOrders = order, 
                            cov = cov, mean = mean, simplify = TRUE)
  expect_identical(moment[[1]], 80)
  
  # test if 'gamma' works for more than one given order
  moments <- symbolicMoments(distribution = 'gamma', 
                             missingOrders = rbind(c(2,1), c(1,3)), 
                             cov = cov, mean = mean, simplify = TRUE)
  expect_identical(moments, c(80, 540))
})

