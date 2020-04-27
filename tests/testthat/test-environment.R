test_that("internal environment", {
  expect_is(multiprobit_env_get(), "environment")
})

test_that("internal environment setup", {
  .onLoad(NULL, "multiprobit")
  logg_post <- mp_log_get()
  expect_equal(length(logg_post), 2)
})

test_that("log messages", {
  logg_pre <- mp_log_get()
  mp_log_message("Test log messages")
  logg_post <- mp_log_get()
  expect_equal(length(logg_post), length(logg_pre) + 1)
  expect_match(logg_post[length(logg_post)], "Test log messages$")
})

test_that("setting/getting options", {
  expect_is(
    as.mp_options(),
    "mp_options"
  )
  expect_is(
    as.mp_options(mp_options()),
    "mp_options"
  )
  expect_is(
    as.mp_options(list()),
    "mp_options"
  )
  expect_is(
    as.mp_options(as.list(mp_options_default())),
    "mp_options"
  )
  expect_is(
    mp_options_set(),
    "mp_options"
  )
  expect_is(
    mp_options_get(),
    "mp_options"
  )

  expect_equal(
    suppressWarnings(mp_options_check(mp_options_default())),
    TRUE
  )
  expect_equal(
    suppressWarnings(mp_options_check(mp_options())),
    FALSE
  )
})
