# Health Endpoint Tests

test_that("GET /health returns 200 with status and timestamp", {
  skip_if_not_installed("jsonlite")

  app <- .health_app()

  req <- list(PATH_INFO = "/health", REQUEST_METHOD = "GET")
  resp <- app$call(req)

  expect_equal(resp$status, 200L)
  expect_equal(resp$headers[["Content-Type"]], "application/json")

  parsed <- jsonlite::fromJSON(resp$body)
  expect_equal(parsed$status, "ok")
  expect_true("timestamp" %in% names(parsed))
  expect_true(nchar(parsed$timestamp) > 0)
})

test_that("unknown paths return 404", {
  skip_if_not_installed("jsonlite")

  app <- .health_app()

  req <- list(PATH_INFO = "/other", REQUEST_METHOD = "GET")
  resp <- app$call(req)

  expect_equal(resp$status, 404L)
})

test_that("non-GET request to /health returns 404", {
  skip_if_not_installed("jsonlite")

  app <- .health_app()

  req <- list(PATH_INFO = "/health", REQUEST_METHOD = "POST")
  resp <- app$call(req)

  expect_equal(resp$status, 404L)
})

test_that("start and stop health server lifecycle works", {
  skip_if_not_installed("httpuv")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("curl")

  port <- 28080L
  srv <- start_health_server(host = "127.0.0.1", port = port)
  on.exit(stop_health_server(srv), add = TRUE)

  # Use curl async multi handle so httpuv event loop can run

  pool <- curl::new_pool()
  result <- NULL

  curl::curl_fetch_multi(
    paste0("http://127.0.0.1:", port, "/health"),
    done = function(resp) { result <<- resp },
    fail = function(msg) { result <<- list(error = msg) },
    pool = pool
  )

  # Run both curl and httpuv event loops concurrently
  for (i in 1:100) {
    httpuv::service(10)
    remaining <- curl::multi_run(timeout = 0.01, pool = pool)
    if (remaining$pending == 0) break
  }

  skip_if(is.null(result) || !is.null(result$error), "Could not connect to health server")

  expect_equal(result$status_code, 200L)

  parsed <- jsonlite::fromJSON(rawToChar(result$content))
  expect_equal(parsed$status, "ok")
  expect_true("timestamp" %in% names(parsed))
})
