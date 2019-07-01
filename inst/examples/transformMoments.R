# Calculate the raw moment of order (1,1,2) for a three-dimensional random variable X:

mList <- momentList(rawMomentOrders = diag(3),
                    rawMoments = list("m1", "m2", "m3"),
                    centralMomentOrders = expand.grid(list(0:1,0:1,0:2)),
                    centralMoments = as.list(c(1, 0, 0, "a", 0, letters[2:8])))

transformMoment(order = c(1,1,2), type = 'raw', 
                momentList = mList, simplify = TRUE)