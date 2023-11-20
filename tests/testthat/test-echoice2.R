test_that("Sampler runs", {
  data(icecream)
  icecream_est <- vd_est_vdm(icecream,R=200, cores=2)
  mycheck=sd(icecream_est$thetaDraw[1,1,])>0
  
  expect_true(mycheck)
})
