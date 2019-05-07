context("momentList")
#t1 tests für validate_momentList

test_that("returned object is of class momentList", {
  
  mList <- new_momentList(rawMomentOrders = cbind(c(3, 1, 5)),
                          rawMoments = list("A", "m1", "B"),
                          centralMomentOrders = cbind(c(0, 1, 2)),
                          centralMoments = list(1, 0, "C"))
  expect_identical(class(mList), "momentList")
  
  mList <- momentList(rawMomentOrders = cbind(c(3, 1, 5)),
                      rawMoments = list("A", "m1", "B"),
                      centralMomentOrders = cbind(c(0, 1, 2)),
                      centralMoments = list(1, 0, "C"))
  expect_identical(class(mList), "momentList")
  
  mList <- validate_momentList(mList)
  expect_identical(class(mList), "momentList")
  
  mList <- transformMoment(order = 3, type = "central", momentList = mList)
  expect_identical(class(mList), "momentList")
})

test_that("mean works", {
  
  mList <- momentList(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                      rawMoments = list("A", "m1", "m2", "m3", "B"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  
  mean <- mean(mList)
  expect_equal(mean, list("m1", "m2", "m3"))
  
  mList <- momentList(rawMomentOrders = rbind(c(1, 0, 0)),
                      rawMoments = list("m1"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  expect_equal(mean(mList), list("m1", NA, NA))
})

test_that("cov works for n > 1", {
  
  mList <- momentList(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                      rawMoments = list("A", "m1", "m2", "m3", "B"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  
  cov <- cov(mList)
  expect_equal(cov, list(list("c", "e", "k"),
                         list("e", "g", "m"),
                         list("k", "m", "s")))
  
  # with use of transformMoment
  mList <- momentList(rawMomentOrders = expand.grid(list(0:2, 0:2)),
                      rawMoments = as.list(letters[1:9]), 
                      warnings = FALSE)
  
  cov <- cov(mList)
  expect_equal(cov, list(list(quote(b^2 * (a - 2) + c), quote (b*d * (a - 2) + e)),
                         list(quote(b*d * (a - 2) + e), quote (d^2 * (a - 2) + g))))
})

test_that("cov works for n = 1", {
  
  mList <- structure(list(rawMomentOrders = cbind(c(3, 1, 5)),
                          rawMoments = list("A", "m1", "B"),
                          centralMomentOrders = cbind(c(0, 1, 2)),
                          centralMoments = list(1, 0, "C")),
                     class = "momentList")
  
  cov <- cov(mList)
  expect_equal(cov, list(list("C")))
  
  # with transform moment
  mList <- momentList(rawMomentOrders = cbind(c(0, 1, 2)),
                      rawMoments = list("A", "m1", "B"),
                      centralMomentOrders = cbind(c(0, 1)),
                      centralMoments = list(1, 0),
                      warnings = FALSE)

  cov <- cov(mList)
  expect_equal(cov, list(list(quote(B + m1^2 * (A - 2)))))
})
