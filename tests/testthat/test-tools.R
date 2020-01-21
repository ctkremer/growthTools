test_that("square root smoother function", {
  expect_equal(sqfunc(0,1,0), 0)
  expect_equal(sqfunc(0,1,1), 1)
  expect_equal(sqfunc(1,1,1), 1.118034)
})
