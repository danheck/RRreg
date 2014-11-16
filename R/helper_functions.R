## 
#' @
is2group <- function(model){
  ifelse(model %in% c("CDM","CDMsym","UQTunknown","SLD","mix.unknown"), TRUE, FALSE)
}