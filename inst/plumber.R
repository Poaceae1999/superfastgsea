# Plumber API definition for superfastgsea health endpoint.
# Start with: superfastgsea::start_health_server()

#* @apiTitle superfastgsea health API

#* Health check endpoint for container orchestration
#* @get /health
#* @serializer json
function() {
  list(
    status = "ok",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}
