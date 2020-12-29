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
  x <- sim_dat$x
  z <- sim_dat$z
  y <- sim_dat$y
  zirf_df <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
  expect_match(class(zirf_df), "zirf_fit")
})


sim_dat <- gen_zip(100, beta = c(.1, -.1), xi = c(-.5, .1))
x <- sim_dat$x
z <- sim_dat$z
y <- sim_dat$y
colnames(x) <- c("34Agg")
xdf <- as.data.frame(x)
zirf_df <- suppressWarnings(zirf_fit(x, z, y, rounds = 10, mtry = 1))
