
# get all available model names
modelnames <- function() {
  c(
    "Warner", "UQTknown", "UQTunknown", "Mangat", "Kuk",
    "FR", "Crosswise", "Triangular", "CDM", "CDMsym", "SLD", "mix.norm",
    "mix.exp", "mix.unknown", "custom"
  )
}

# check if an RR model needs two groups
is2group <- function(model) {
  model %in% c("CDM", "CDMsym", "UQTunknown", "SLD", "mix.unknown")
}

# check if RR model uses continuous responses
isContinuous <- function(model) {
  model %in% c("mix.exp", "mix.unknown", "mix.norm")
}

# fix second parameter for 2 group models
fix.par2 <- function(model) {
  par2 <- switch(model,
    "SLD" = c(t = 1),
    "UQTunknown" = c(piUQ = .50),
    "CDM" = c(gamma = 0),
    "CDMsym" = c(gamma = 0)
  )
  return(par2)
}
