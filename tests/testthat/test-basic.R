test_that("basic workflow works", {
    set.seed(10)
  # Random competing risks data
  time  <- rexp(500)
  event <- sample(c(2, 4, 7, 8), 500, replace = TRUE)
  z <- cbind(z1 = rbinom(500, 1, 0.3),
             z2 = rnorm(500))
  
  # Initial values:
  # (alpha, rho, tau, beta1, beta2) for each event
  user_init <- c(1, -1, 0.5,  1,  1,      # event 7 : event of interest (failcode)
                 0.01, -1, 0.25, -1, 0,   # event 4 : competing event
                 1, -1, 0.25, -1, -1)     # event 8 : competing event
  # pcrr test
  fit2 <- pcrr(ftime = time, fstatus = event, cov = z, 
               failcode = 7, cencode = 2, init = user_init)
  expect_s3_class(fit2, "pcrr")
  expect_true(fit2$converged)
  expect_equal(length(fit2$coef), length(user_init))
  expect_true(is.finite(fit2$loglik))
  expect_false(anyNA(fit2$coef))
  
  # predict test
  cov <- rbind(c(0, 0.15), c(1, -0.30), c(0, 0.70))
  
  pred2 <- predict(fit2, cov)
  expect_true(is.matrix(pred2$pred))
  expect_equal(nrow(pred2$pred), 200)
  expect_equal(ncol(pred2$pred), nrow(cov) + 1)
  
  expect_true(all(pred2$pred[, -1] >= 0))
  expect_true(all(pred2$pred[, -1] <= 1))
  
  for(i in 2:ncol(pred2$pred)){
    expect_true(all(diff(pred2$pred[, i]) >= -1e-8))
  }
  
  # cure
  
  cure2 <- cure(fit2, cov)
  
  expect_type(cure2$cure, "double")
  expect_length(cure2$cure, nrow(cov))
  expect_true(all(cure2$cure >= 0 | is.na(cure2$cure)))
  expect_true(all(cure2$cure <= 1 | is.na(cure2$cure)))
})
