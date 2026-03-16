#' Build Health App Handler
#'
#' Creates the httpuv app list used by the health server.
#' Exported for testing purposes.
#'
#' @return A list with a \code{call} element suitable for \code{httpuv::startServer}.
#' @keywords internal
.health_app <- function() {
  list(
    call = function(req) {
      if (req$PATH_INFO == "/health" && req$REQUEST_METHOD == "GET") {
        body <- jsonlite::toJSON(
          list(
            status = jsonlite::unbox("ok"),
            timestamp = jsonlite::unbox(format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"))
          ),
          auto_unbox = FALSE
        )
        list(
          status = 200L,
          headers = list("Content-Type" = "application/json"),
          body = body
        )
      } else {
        list(
          status = 404L,
          headers = list("Content-Type" = "application/json"),
          body = jsonlite::toJSON(list(error = jsonlite::unbox("not found")), auto_unbox = FALSE)
        )
      }
    }
  )
}

#' Start Health Check Server
#'
#' Starts a lightweight HTTP server that exposes a \code{/health} endpoint
#' returning JSON status information. Useful for container health checks.
#'
#' @param host Host to bind to (default: "0.0.0.0")
#' @param port Port to listen on (default: 8080)
#'
#' @return An \code{httpuv} server handle (invisible). Use
#'   \code{\link{stop_health_server}} to shut it down.
#'
#' @examples
#' \dontrun{
#' srv <- start_health_server(port = 8080)
#' # GET http://localhost:8080/health -> {"status":"ok","timestamp":"..."}
#' stop_health_server(srv)
#' }
#'
#' @export
start_health_server <- function(host = "0.0.0.0", port = 8080) {
  if (!requireNamespace("httpuv", quietly = TRUE)) {
    stop("Package 'httpuv' is required. Install with: install.packages('httpuv')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Install with: install.packages('jsonlite')")
  }

  app <- .health_app()
  server <- httpuv::startServer(host, port, app)
  message(sprintf("Health server running on http://%s:%d/health", host, port))
  invisible(server)
}

#' Stop Health Check Server
#'
#' Stops a health check server previously started with
#' \code{\link{start_health_server}}.
#'
#' @param server Server handle returned by \code{\link{start_health_server}}
#'
#' @return \code{NULL} (invisible)
#'
#' @examples
#' \dontrun{
#' srv <- start_health_server(port = 8080)
#' stop_health_server(srv)
#' }
#'
#' @export
stop_health_server <- function(server) {
  server$stop()
  message("Health server stopped")
  invisible(NULL)
}
