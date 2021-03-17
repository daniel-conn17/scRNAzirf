test_that("zirf_fit works on matrices", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- sim_dat$x
  z <- sim_dat$z
  y <- sim_dat$y
  zirf_mat <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
  expect_match(class(zirf_mat), "zirf_fit")
})

test_that("zirf_fit works on data.frame", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- as.data.frame(sim_dat$x)
  z <- as.data.frame(sim_dat$z)
  y <- sim_dat$y
  zirf_df <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
  expect_match(class(zirf_df), "zirf_fit")
})

test_that("zirf_fit works when column names start with numbers (matrices)", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- sim_dat$x
  colnames(x) <- c("34Agg")
  z <- sim_dat$z
  y <- sim_dat$y
  zirf_mat <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
  expect_match(class(zirf_mat), "zirf_fit")
})

test_that("zirf_fit works when column names start with numbers (data.frames)", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- sim_dat$x
  x <- as.data.frame(x)
  names(x) <- c("34Agg")
  z <- as.data.frame(sim_dat$z)
  y <- sim_dat$y
  zirf_df <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
  expect_match(class(zirf_df), "zirf_fit")
})

test_that("zirf_fit fails when class(x) and class(z) don't match", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- sim_dat$x
  names(x) <- c("34Agg")
  z <- as.data.frame(sim_dat$z)
  y <- sim_dat$y
  zirf_df <- tryCatch(suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1)),
                      error = function(e) "correct error")
  expect_match(zirf_df, "correct error")
})

test_that("zirf_fit fails when class(x) and class(z) don't match", {
  sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
  x <- sim_dat$x
  names(x) <- c("34Agg")
  z <- as.data.frame(sim_dat$z)
  y <- sim_dat$y
  zirf_df <- tryCatch(suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1)),
                      error = function(e) "correct error")
  expect_match(zirf_df, "correct error")
})


sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
x <- sim_dat$x
colnames(x) <- c("34Agg")
x <- sim_dat$x
z <- sim_dat$z
y <- sim_dat$y
y <- as.data.frame(y)
sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
newx <- sim_dat$x
newz <- sim_dat$z
zirf_df <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1, newx = newx,
                                     newz = newz))

sim_dat <- gen_zip(100, beta = c(.1, -.1, 0, .1, 0, -.1, 0), xi = c(-.5, .1))
x <- t(sim_dat$x)
z <- sim_dat$z
count_mat <- matrix(rpois(6*100, 1), 6, 100)
rownames(count_mat) <- paste0("V", 1:6)
ff <- zirf_genie3(exprMatrix = x, countMatrix = count_mat, z = z, nCores = 2,
                  regulators = c("V1", "V2"), targets = c("V1", "V2", "V3", "V4"))


ff <- zirf_genie3(exprMatrix = x, countMatrix = count_mat, z = z, nCores = 2)



#to do:
#write function that computes all zirf models (make sure it uses parallelization
#write functions that incorporate SCENIC functions
#write SCENIC vignette
