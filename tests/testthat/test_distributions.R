context("distributions")

test_that("mtrunc works as expected", {
  expect_identical(mtrunc(spec = "exp", order = 5), 
                   actuar::mexp(order = 5))
  expect_identical(mtrunc(spec = "lnorm", order = 3, meanlog = 4), 
                   actuar::mlnorm(order = 3, meanlog = 4))
  expect_equal(mtrunc(spec = "norm", order = 1:2, a = -2, b = 2),
               c(0, 0.774), tolerance = 1e-3)
  expect_error(mtrunc(spec = "unif", order = 1:3, a = 3, b = 1))
  expect_error(mtrunc(spec = "unif", order = 1:3, a = 3, b = 1, 
                      min = 0, max = 1))
})

test_that("dtrunc works as expected", {
  # see http://lagrange.math.siu.edu/Olive/ch4.pdf (with λ = 1/rate)
  x <- seq(0, 2, 0.1)
  expect_identical(dtrunc(x, spec = "exp", b = 5),
                  (exp(-x)/(1-exp(-5))))
  expect_identical(dtrunc(x, spec = "exp", rate = 2, b = 5),
                   2*exp(-x*2)/(1-exp(-10)))
  expect_warning(d <- dtrunc(x, spec = "norm", a = 20, b = 30))
  expect_identical(d, rep(0, length(x)))
  expect_error(dtrunc(x, spec = "unif", a = 30, b = 0))
})

test_that("dmix works as expected", {
  x <- seq(0, 2, 0.1)
  expect_identical(dmix(distrib = list(list(spec = "norm", mean = 3)))(x),
                   dnorm(x, mean = 3))
  expect_identical(dmix(distrib = list(list(spec = "lnorm")), lower = 2)(x),
                   dtrunc(x, spec = "lnorm", a = 2))
  expect_equal(dmix(distrib = list(list(spec = "unif", min = 0, max = 1-1e-9),
                                   list(spec = "unif", min = 1, max = 2)))(x),
                   dunif(x, min = 0, max = 2))
})

test_that("mmix works as expected", {
  expect_identical(mmix(distrib = list(list(spec = "norm", mean = 3)), 
                   order = 2),
                   actuar::mnorm(mean = 3, order = 2))
  expect_identical(mmix(distrib = list(list(spec = "lnorm")), 
                   order = 2:5, lower = 2),
                   mtrunc(order = 2:5, spec = "lnorm", a = 2))
  expect_equal(mmix(distrib = list(list(spec = "unif", min = 0, max = 1-1e-9),
                                   list(spec = "unif", min = 1, max = 2)), 
                    order = 1:4),
               actuar::munif(order = 1:4, min = 0, max = 2))
})
