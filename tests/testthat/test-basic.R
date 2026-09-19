set.seed(10)
ftime_   <- rexp(500)
fstatus_ <- sample(c(2, 4, 7, 8), 500, replace = TRUE)
zmat     <- cbind(z1 = rbinom(500, 1, 0.3),
                  z2 = rnorm(500))

## (alpha, rho, tau, beta1, beta2) for each event
user_init <- c(1,    -1, 0.5,   1,  1,    # event 7 : failcode
               0.01, -1, 0.25, -1,  0,    # event 4 : competing
               1,    -1, 0.25, -1, -1)    # event 8 : competing

newz <- rbind(c(0, 0.15), c(1, -0.30), c(0, 0.70))

fit2  <- pcrr(ftime = ftime_, fstatus = fstatus_, cov = zmat,
              failcode = 7, cencode = 2, init = user_init)
pred2 <- predict(fit2, newz, case = 1)
cure2 <- cure(fit2, newz, case = 1)



test_that("pcrr fits and returns a well-formed object", {
  expect_s3_class(fit2, "pcrr")
  expect_true(fit2$converged)
  expect_equal(length(fit2$coef), length(user_init))
  expect_true(is.finite(fit2$loglik))
  expect_false(anyNA(fit2$coef))
  
  expect_equal(fit2$k, 3L)
  expect_equal(fit2$p, ncol(zmat))
  expect_equal(fit2$mapping, c(7, 4, 8))          # failcode first
  expect_equal(nrow(fit2$case_all), length(fit2$case_model))
})



test_that("pcrr stores the data needed by plot.pcrr", {
  expect_false(is.null(fit2$x))
  expect_false(is.null(fit2$delta))
  expect_false(is.null(fit2$z))
  
  expect_length(fit2$x, 500)
  expect_equal(dim(fit2$delta), c(500L, 3L))
  expect_equal(dim(fit2$z), c(500L, ncol(zmat)))
})



test_that("predict returns the expected nested structure", {
  expect_s3_class(pred2, "predict.pcrr")
  expect_length(pred2$pred, length(pred2$case))
  expect_length(pred2$pred[[1]], length(pred2$event))
  
  cif <- pred2$pred[[1]][[1]]
  expect_true(is.matrix(cif))
  expect_equal(nrow(cif), 200)
  expect_equal(ncol(cif), nrow(newz) + 1)
})


test_that("predicted CIFs are valid probabilities", {
  cif <- pred2$pred[[1]][[1]]
  
  expect_true(all(cif[, -1] >= 0))
  expect_true(all(cif[, -1] <= 1))
  expect_equal(unname(cif[1, -1]), rep(0, nrow(newz)))   # F(0) = 0
  
  for (i in 2:ncol(cif)) {
    expect_true(all(diff(cif[, i]) >= -1e-8))            # nondecreasing
  }
  
  expect_equal(cif[, 1], sort(cif[, 1]))                 # time column sorted
})



test_that("cure returns valid fractions and matching structure", {
  expect_s3_class(cure2, "cure.pcrr")
  
  cf <- cure2$cure[[1]][[1]]
  st <- cure2$status[[1]][[1]]
  tb <- cure2$t_boundary[[1]][[1]]
  
  expect_type(cf, "double")
  expect_length(cf, nrow(newz))
  expect_length(st, nrow(newz))
  expect_length(tb, nrow(newz))
  
  expect_true(all(cf >= 0 | is.na(cf)))
  expect_true(all(cf <= 1 | is.na(cf)))
  expect_true(all(st %in% c("cure", "no cure (asymptotic)",
                            "no cure (finite support)")))
  
  ## NA cure fraction occurs exactly on the finite-support profiles
  expect_equal(is.na(cf), st == "no cure (finite support)")
})



test_that("summary produces finite standard errors", {
  s <- summary(fit2, case = 1)
  expect_s3_class(s, "summary.pcrr")
  
  cf <- s$coef[[1]]
  expect_true(is.matrix(cf))
  expect_equal(nrow(cf), fit2$k * fit2$p)
  expect_false(anyNA(cf[, "se(coef)"]))
  expect_true(all(cf[, "se(coef)"] > 0))
  
  ci <- s$conf_int[[1]]
  expect_true(all(ci[, 3] <= ci[, 4]))          # lower <= upper
})



test_that("print methods run without error", {
  expect_output(print(fit2))
  expect_output(print(summary(fit2, case = 1)))
  expect_output(print(pred2), "obs 1")
  expect_output(print(cure2), "obs")
})



test_that("invalid case and event are rejected", {
  n_case <- nrow(fit2$case_all)
  expect_error(predict(fit2, newz, case = n_case + 1))
  expect_error(predict(fit2, newz, case = 0))
  expect_error(predict(fit2, newz, case = 1.5))
  expect_error(predict(fit2, newz, event = 999))
  expect_error(predict(fit2, newz[, 1, drop = FALSE]))   # wrong ncol
})