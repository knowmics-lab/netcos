library(R6)
library(log4r)
source("modules/config.R")
source("modules/obj_is_na.R")

TSRLogger <- R6Class(
  "TSRLogger",
  public = list(
    initialize = function(log_file = NA) {
      private$console_logger <- logger(threshold = "DEBUG")
      if (obj_is_na(log_file)) {
        log_file <- tsr_logger_cfg$log_file_name
      }
      private$file_logger <- logger(threshold = "DEBUG", appenders = file_appender(log_file))
      return(NA)
    },
    log = function(message) {
      info(private$file_logger, message)
      info(private$console_logger, message)
    }
  ),
  private = list(
    file_logger = NA,
    console_logger = NA
  )
)