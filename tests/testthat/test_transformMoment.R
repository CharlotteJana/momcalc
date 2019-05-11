context("transformMoment")

test_that("transformMoment works for n = 3, type = 'raw' and all central moments given", {
  
  mList <- structure(list(rawMomentOrders = diag(3),
                          rawMoments = list("m1", "m2", "m3"),
                          centralMomentOrders = expand.grid(list(0:1, 0:1, 0:2)),
                          centralMoments = as.list(letters[1:12])),
                     class = "momentList")
  
  res1 <- transformMoment(order = c(1,1,2),
                          type = 'raw',
                          closure = "",
                          momentList = mList,
                          simplify = FALSE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^1 * "m3"^2) * ("a") + 
                     1 * ("m1"^0 * "m2"^1 * "m3"^2) * ("b") + 
                     1 * ("m1"^1 * "m2"^0 * "m3"^2) * ("c") + 
                     1 * ("m1"^0 * "m2"^0 * "m3"^2) * ("d") + 
                     2 * ("m1"^1 * "m2"^1 * "m3"^1) * ("e") + 
                     2 * ("m1"^0 * "m2"^1 * "m3"^1) * ("f") + 
                     2 * ("m1"^1 * "m2"^0 * "m3"^1) * ("g") + 
                     2 * ("m1"^0 * "m2"^0 * "m3"^1) * ("h") + 
                     1 * ("m1"^1 * "m2"^1 * "m3"^0) * ("i") + 
                     1 * ("m1"^0 * "m2"^1 * "m3"^0) * ("j") + 
                     1 * ("m1"^1 * "m2"^0 * "m3"^0) * ("k") + 
                     1 * ("m1"^0 * "m2"^0 * "m3"^0) * ("l"))
  
  res2 <- structure(list(rawMomentOrders = rbind(c(0, 0, 0),
                                                 diag(3), 
                                                 c(1, 1, 2)),
                         rawMoments = list(1, "m1", "m2", "m3", momentRaw),
                         centralMomentOrders = mList$centralMomentOrders,
                         centralMoments = mList$centralMoments),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  
  
})


test_that("transformMoment works for n = 3, type = 'central' and all raw moments given", {
  
  mList <- structure(list(rawMomentOrders = expand.grid(list(0:1, 0:1, 0:2)),
                          rawMoments = as.list(letters[1:12]),
                          centralMomentOrders = matrix(, nrow = 0, ncol = 3),
                          centralMoments = list()),
                     class = "momentList")
  
  res1 <- transformMoment(order = c(1,1,2),
                          type = 'central',
                          closure = "",
                          momentList = mList,
                          simplify = FALSE)
  
  momentCentr <- quote(1 * ("b"^1 * "c"^1 * "e"^2) * ("a") + 
                       -1 * ("b"^0 * "c"^1 * "e"^2) * ("b") + 
                       -1 * ("b"^1 * "c"^0 * "e"^2) * ("c") + 
                       1 * ("b"^0 * "c"^0 * "e"^2) * ("d") + 
                       -2 * ("b"^1 * "c"^1 * "e"^1) * ("e") + 
                       2 * ("b"^0 * "c"^1 * "e"^1) * ("f") + 
                       2 * ("b"^1 * "c"^0 * "e"^1) * ("g") + 
                       -2 * ("b"^0 * "c"^0 * "e"^1) * ("h") + 
                       1 * ("b"^1 * "c"^1 * "e"^0) * ("i") + 
                       -1 * ("b"^0 * "c"^1 * "e"^0) * ("j") + 
                       -1* ("b"^1 * "c"^0 * "e"^0) * ("k") + 
                       1 * ("b"^0 * "c"^0 * "e"^0) * ("l"))
  
  res2 <- structure(list(rawMomentOrders = mList$rawMomentOrders,
                         rawMoments = mList$rawMoments,
                         centralMomentOrders = rbind(c(0, 0, 0),
                                                     diag(3), 
                                                     c(1, 1, 2)),
                         centralMoments = list(1, 0, 0, 0, momentCentr)),
                    class = "momentList")
  expect_equal(res1, res2, check.attributes = FALSE)
})

test_that("transformMoment works with only the last moment(s) missing", {
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cOrders <- cOrders[-6, ]
  cMoments <- list(1, 0, 0, "a", "b")
  
  rOrders <- expand.grid(list(0:2, 0:2))
  rOrders <- rOrders[rowSums(rOrders) <= 2, ]
  rMoments <- list(1, "m1", "A", "m2", "B", "C")
  
  mList <- structure(list(rawMomentOrders = rOrders,
                          rawMoments = rMoments,
                          centralMomentOrders = cOrders,
                          centralMoments = cMoments),
                     class = "momentList")
  
  #------- (1,2) missing in rOrders and cOrders, type = 'raw' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  res1Simplified <- transformMoment(order = c(1,2),
                                    type = 'raw',
                                    closure = "zero",
                                    momentList = mList,
                                    simplify = TRUE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^2) * (1) + 
                     1 * ("m1"^0 * "m2"^2) * (0) + 
                     2 * ("m1"^1 * "m2"^1) * (0) + 
                     2 * ("m1"^0 * "m2"^1) * ("a") + 
                     1 * ("m1"^1 * "m2"^0) * ("b") + 
                     1 * ("m1"^0 * "m2"^0) * (0))
  
  momentRawSimplified <- quote(b * m1 + m2 * (2 * a + m1 * m2))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  res2Simplified <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 2)),
                         rawMoments = append(rMoments, momentRawSimplified),
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  
  expect_equal(res1, res2, check.attributes = FALSE)
  expect_equal(res1Simplified, res2Simplified, check.attributes = FALSE)
  
  #------- (1,2) missing in rOrders and cOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  res2 <- structure(list(rawMomentOrders = rOrders,
                         rawMoments = rMoments,
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) only missing in rOrders, type = 'raw' ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cMoments <- list(1, 0, 0, "a", "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^2) * (1) + 
                       1 * ("m1"^0 * "m2"^2) * (0) + 
                       2 * ("m1"^1 * "m2"^1) * (0) + 
                       2 * ("m1"^0 * "m2"^1) * ("a") + 
                       1 * ("m1"^1 * "m2"^0) * ("b") + 
                       1 * ("m1"^0 * "m2"^0) * ("c"))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = cOrders,
                         centralMoments = cMoments),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) only missing in rOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  expect_equal(res1, mList, check.attributes = FALSE)
  
  #------- (1,2) only missing in cOrders, type = 'raw' ----------
  
  cOrders <- cOrders[-6, ]
  cMoments <- cMoments[-6]
  rOrders <- rbind(rOrders, c(1, 2))
  rMoments <- append(rMoments, "D")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  mList$rawMomentOrders <- rOrders
  mList$rawMoments <- rMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  expect_equal(res1, mList, check.attributes = FALSE)

  #------- (1,2) only missing in cOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentCentr <- quote(-1 * ("m1"^1 * "m2"^2) * (1) + 
                       1 * ("m1"^0 * "m2"^2) * ("m1") + 
                       2 * ("m1"^1 * "m2"^1) * ("m2") + 
                       -2 * ("m1"^0 * "m2"^1) * ("B") + 
                       -1 * ("m1"^1 * "m2"^0) * ("C") + 
                       1 * ("m1"^0 * "m2"^0) * ("D"))
  
  res2 <- structure(list(rawMomentOrders = rOrders,
                         rawMoments = rMoments,
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, momentCentr)),
                    class = "momentList")
                    
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) not missing at all ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cMoments <- list(1, 0, 0, "a", "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  res2 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList)
  
  expect_equal(res1, mList, check.attributes = FALSE)
  expect_equal(res2, mList, check.attributes = FALSE)
})

test_that("transformMoment works with a missing moment in the middle and at the end", {
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cOrders <- cOrders[-6, ]
  cMoments <- list(1, 0, 0, "a", "b")
  
  rOrders <- matrix(c(0, 0,
                      1, 0,
                      2, 0,
                      0, 1,
                      0, 2,
                      1, 2), ncol = 2, byrow = TRUE)
  rMoments <- list(1, "m1", "A", "m2", "B", "C")
  
  mList <- structure(list(rawMomentOrders = rOrders,
                          rawMoments = rMoments,
                          centralMomentOrders = cOrders,
                          centralMoments = cMoments),
                     class = "momentList")
  
  #------- (1,2) missing in cOrders and (1, 1) missing in rOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^1) * (1) + 
                     1 * ("m1"^0 * "m2"^1) * (0) + 
                     1 * ("m1"^1 * "m2"^0) * (0) + 
                     1 * ("m1"^0 * "m2"^0) * ("a"))
  
  momentCentr <- bquote(-1 * ("m1"^1 * "m2"^2) * (1) + 
                       1 * ("m1"^0 * "m2"^2) * ("m1") + 
                       2 * ("m1"^1 * "m2"^1) * ("m2") + 
                       -2 * ("m1"^0 * "m2"^1) * (.(momentRaw)) + 
                       -1 * ("m1"^1 * "m2"^0) * ("B") + 
                       1 * ("m1"^0 * "m2"^0) * ("C"))
  
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 1)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, momentCentr)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) and (1,1) missing in cOrders and (1,1) missing in rOrders, type = 'central' ----------
  
  cOrders <- matrix(c(0, 0,
                      1, 0,
                      0, 1,
                      0, 2), ncol = 2, byrow = TRUE)
  cMoments <- list(1, 0, 0, "b")
  
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^1) * (1) + 
                       1 * ("m1"^0 * "m2"^1) * (0) + 
                       1 * ("m1"^1 * "m2"^0) * (0) + 
                       1 * ("m1"^0 * "m2"^0) * (0))
  
  momentCentr <- bquote(-1 * ("m1"^1 * "m2"^2) * (1) + 
                          1 * ("m1"^0 * "m2"^2) * ("m1") + 
                          2 * ("m1"^1 * "m2"^1) * ("m2") + 
                          -2 * ("m1"^0 * "m2"^1) * (.(momentRaw)) + 
                          -1 * ("m1"^1 * "m2"^0) * ("B") + 
                          1 * ("m1"^0 * "m2"^0) * ("C"))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 1)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,1), c(1,2)),
                         centralMoments = append(cMoments, list(0, momentCentr))),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,1) missing in cOrders and (1, 2) missing in rOrders, type = 'raw' ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cOrders <- cOrders[-4, ]
  cMoments <- list(1, 0, 0, "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  rOrders <- expand.grid(list(0:1, 0:2))
  rOrders <- rOrders[-6, ]
  rMoments <- list(1, "m1", "m2", "A", "B")
  mList$rawMomentOrders <- rOrders
  mList$rawMoments <- rMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentCentr <- quote(1 * ("m1"^1 * "m2"^1) * (1) + 
                       -1 * ("m1"^0 * "m2"^1) * ("m1") + 
                       -1 * ("m1"^1 * "m2"^0) * ("m2") + 
                       1 * ("m1"^0 * "m2"^0) * ("A"))
  
  momentRaw <- bquote(1 * ("m1"^1 * "m2"^2) * (1) + 
                      1 * ("m1"^0 * "m2"^2) * (0) + 
                      2 * ("m1"^1 * "m2"^1) * (0) + 
                      2 * ("m1"^0 * "m2"^1) * (.(momentCentr)) + 
                      1 * ("m1"^1 * "m2"^0) * ("b") + 
                      1 * ("m1"^0 * "m2"^0) * ("c"))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1,2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,1)),
                         centralMoments = append(cMoments, momentCentr)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,1) missing in cOrders and (1,1) and (1,2) missing in rOrders, type = 'raw' ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cOrders <- cOrders[-4, ]
  cMoments <- list(1, 0, 0, "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  rOrders <- expand.grid(list(0:1, 0:2))
  rOrders <- rOrders[c(-4, -6), ]
  rMoments <- list(1, "m1", "m2", "A")
  mList$rawMomentOrders <- rOrders
  mList$rawMoments <- rMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList,
                          simplify = FALSE)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^2) * (1) + 
                     1 * ("m1"^0 * "m2"^2) * (0) + 
                     2 * ("m1"^1 * "m2"^1) * (0) + 
                     2 * ("m1"^0 * "m2"^1) * (0) + 
                     1 * ("m1"^1 * "m2"^0) * ("b") + 
                     1 * ("m1"^0 * "m2"^0) * ("c"))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1,2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,1)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
})

test_that("transformMoment needs raw moments of order 1 as input", {
  mList <- structure(list(rawMomentOrders = 2*diag(3),
                          rawMoments = list("m1", "m2", "m3"),
                          centralMomentOrders = expand.grid(list(0:1, 0:1, 0:2)),
                          centralMoments = as.list(letters[1:12])),
                     class = "momentList")
  
  expect_error(suppressWarnings(
    transformMoment(order = c(1,1,2),
                    type = 'raw',
                    closure = "",
                    momentList = mList)
  ))
  
  expect_error(suppressWarnings(
    transformMoment(order = c(5,3,2),
                    type = 'central',
                    closure = "",
                    momentList = mList)
  ))
})