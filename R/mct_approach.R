
# FIXME: write documentation
mct_inflexionpoints <- function(Q, periods = c(100, 300), n = 30L, return_derivative = FALSE) {
  dQ   <- c(diff(Q), NA)
  dQs  <- roll_mean(dQ, n=n[1L], fill=NA)
  p <- periods[1L]:periods[2L]
  inflex_points <- c(which.max(dQs[p]), which.min(dQs[p])) + p[1L] - 1L
  if (return_derivative) return(list(data = data.frame(dQ = dQ, dQs = dQs), inflexion_points = inflex_points))
  return(inflex_points)
}

# FIXME: write documentation
# TODO: check dims of inputs!
mct_snowstorage <- function(Pcum, Qcum, inflexion_points, 
                            use_tangentes = FALSE, 
                            period = c(1, length(Pcum)),
                            fixed_intercept = FALSE,
                            return_tangents = FALSE) {
  if (use_tangentes) {
    x <- period[1L]:period[2L]
    if (fixed_intercept) {
      Pslp <- c(0, lm(Pcum[x]~x+0)$coefficients) # FIXME: this is slow, find a way using .lm.fit() function
    } else {
      Pslp <- .hsaFastLm(x = x, y = Pcum[x])$coefficients
    }
    b <- Qcum[inflexion_points] - (Pslp[2L] * inflexion_points)
    if (return_tangents) {
      return(list(S = c(Sa = Pslp[1L] - b[1L],
                        Sm = b[2L] - b[1L]),
                  P_slp = Pslp,
                  Qa_slp = c(b[1L], Pslp[2L]),
                  Qm_slp = c(b[2L], Pslp[2L])))
    } else {
      return(c(Sa = Pslp[1L] - b[1L],
               Sm = b[2L] - b[1L]))
    }
  } else {
    return(Pcum[inflexion_points[1L]] - Qcum[inflexion_points[1L]])
  }
}

