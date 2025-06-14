test_that("config environment exists", {
    expect_true(exists(".pkg_env", envir = asNamespace("stanflowccus")))
})


test_that("config environment exists", {
    env <- get(".pkg_env", envir = asNamespace("stanflowccus"))
    expect_true(is.list(env$config))
})
