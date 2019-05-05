context("dBEGG")

test_that("dBEGG equals normal distribution for special parameters", {
  x <- seq(-5, 10, by = 0.1)
  expect_equal(dnorm(x), 
               dBEGG(x, alpha = 2, beta = 1, delta0 = 0, 
                     delta1 = 0, eta = 2, eps = 0))
  expect_equal(mnorm(order = 1:4), 
               mBEGG(order = 1:4, 
                     alpha = 2, beta = 1, delta0 = 0, 
                     delta1 = 0, eta = 2, eps = 0))
  
})




