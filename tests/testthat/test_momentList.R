#t3 test that validate creates trivial moments: names of data.frames are not identical

context("momentList")

test_that("returned object is of class momentList", {
  
  mList <- new_momentList(rawMomentOrders = cbind(c(3, 1, 5)),
                          rawMoments = list("A", "m1", "B"),
                          centralMomentOrders = cbind(c(0, 1, 2)),
                          centralMoments = list(1, 0, "C"))
  expect_identical(class(mList), "momentList")
  
  mList <- momentList(rawMomentOrders = data.frame(cbind(c(3, 1, 5))),
                      rawMoments = list("A", "m1", "B"),
                      centralMomentOrders = cbind(c(0, 1, 2)),
                      centralMoments = list(1, 0, "C"))
  expect_identical(class(mList), "momentList")
  
  mList <- validate_momentList(mList)
  expect_identical(class(mList), "momentList")
  
  mList <- transformMoment(order = 3, type = "central", momentList = mList)
  expect_identical(class(mList), "momentList")
})

test_that("extractMean works", {
  
  mList <- momentList(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                      rawMoments = list("A", "m1", "m2", "m3", "B"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  
  mean <- extractMean(mList)
  expect_equal(mean, list("m1", "m2", "m3"))
  
  mList <- momentList(rawMomentOrders = rbind(c(1, 0, 0)),
                      rawMoments = list("m1"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  expect_equal(extractMean(mList), list("m1", NA, NA))
})

test_that("extractCov works for n > 1", {
  
  mList <- momentList(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                      rawMoments = list("A", "m1", "m2", "m3", "B"),
                      centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                      centralMoments = append(as.list(letters), "lastElement"),
                      warnings = FALSE)
  
  cov <- extractCov(mList)
  expect_equal(cov, list(list("c", "e", "k"),
                         list("e", "g", "m"),
                         list("k", "m", "s")))
  
  # with use of transformMoment
  mList <- momentList(rawMomentOrders = expand.grid(list(0:2, 0:2)),
                      rawMoments = as.list(letters[1:9]), 
                      warnings = FALSE)
  
  cov <- extractCov(mList)
  expect_equal(cov, list(list(quote(b^2 * (a - 2) + c), quote (b*d * (a - 2) + e)),
                         list(quote(b*d * (a - 2) + e), quote (d^2 * (a - 2) + g))))
})

test_that("extractCov works for n = 1", {
  
  mList <- structure(list(rawMomentOrders = cbind(c(3, 1, 5)),
                          rawMoments = list("A", "m1", "B"),
                          centralMomentOrders = cbind(c(0, 1, 2)),
                          centralMoments = list(1, 0, "C")),
                     class = "momentList")
  
  cov <- extractCov(mList)
  expect_equal(cov, list(list("C")))
  
  # with transform moment
  mList <- momentList(rawMomentOrders = cbind(c(0, 1, 2)),
                      rawMoments = list("A", "m1", "B"),
                      centralMomentOrders = cbind(c(0, 1)),
                      centralMoments = list(1, 0),
                      warnings = FALSE)

  cov <- extractCov(mList)
  expect_equal(cov, list(list(quote(B + m1^2 * (A - 2)))))
})

test_that("validateMomentlist throws errors", {
  
  expect_error(momentList())
  expect_error(validate_momentList(new_momentList()))
  expect_error(validate_momentList(new_momentList(rawMomentOrders = diag(3))))
  expect_error(validate_momentList(new_momentList(centralMomentOrders = diag(2))))
  expect_error(validate_momentList(new_momentList(rawMomentOrders = diag(3),
                                                  rawMoments = list("a"))))
  expect_error(validate_momentList(new_momentList(centralMoments = list("A"))))
  expect_error(validate_momentList(new_momentList(rawMomentOrders = diag(2),
                                                  centralMoments = list("A", "B"))))
  expect_error(validate_momentList(new_momentList(rawMomentOrders = diag(3),
                                                  rawMoments = list("a", "b", "c"),
                                                  centralMomentOrders = diag(2),
                                                  centralMoments = list("A", "B"))))
  expect_error(validate_momentList(new_momentList(rawMomentOrders = rbind(diag(2), diag(2)),
                                                  rawMoments = list("a", "b", "a", "b"))))
})

test_that("validateMomentlist creates trivial moments", {
  
  
  mList <- validate_momentList(new_momentList(rawMomentOrders = diag(2),
                                              rawMoments = list("a", "b")))
  mList2 <- new_momentList(rawMomentOrders = rbind(c(0, 0), diag(2)),
                           rawMoments = list(1, "a", "b"),
                           centralMomentOrders = data.frame(rbind(c(0, 0), diag(2))),
                           centralMoments <- list(1, 0, 0))
  expect_equal(mList, mList2, check.names = FALSE)
  
  expect_warning(validate_momentList(new_momentList(centralMomentOrders = diag(2),
                                                    centralMoments = list("A", "B"))))
                               
  
})