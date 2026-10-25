# Owner-approved five-condition mapping; original specimen/time metadata untouched.
signature_five_conditions <- function(metadata) {
  display_time <- as.character(metadata$time_point)
  exception <- as.character(metadata$sample) == "REVO31-N8"
  if (any(exception & (display_time != "T8H" | as.character(metadata$pressure) != "negative"))) {
    stop("Approved REVO31-N8 display exception requires original T8H/negative metadata.")
  }
  display_time[exception] <- "T10H"
  condition_map <- c(T0H_None = "Control", T4H_positive = "P4", T4H_negative = "N4",
    T10H_positive = "P10", T10H_negative = "N10")
  result <- unname(condition_map[paste(display_time, as.character(metadata$pressure), sep = "_")])
  if (anyNA(result)) stop("Unexpected original time/pressure tuple in approved five-condition view.")
  result
}
