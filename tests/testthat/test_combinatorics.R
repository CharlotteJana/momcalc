context("combinatorics")

test_that("dfactorial works as expected", {
  expect_equal(dfactorial(5), 5*3*1)
  expect_equal(dfactorial(8), 8*6*4*2)
  expect_equal(dfactorial(0), 1)
  expect_error(dfactorial(-3))
  
})

test_that("multinomial works as expected", {
  expect_equal(multinomial(3, 2), 3*2/2)
  expect_equal(multinomial(5, c(2,4)), (5*4*3*2*1)/(2*1*4*3*2*1))
})
