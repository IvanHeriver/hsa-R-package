###############################################################################
#
# Hydrological Signatures Analysis
#
# -----------------------------------------------------------------------------
#
# A set of functions to compute various hydrological signatures.
#
###############################################################################
# -----------------------------------------------------------------------------
# Important notes: 
# (1) The 'RcppSignatures.cpp' files is required and should be located in the  
# same folder as this current R file.
# (2) R development tools are also required and should be installed beforehand.
# Note that sourcing this file in R studio will prompt you into installing them.
#
# -----------------------------------------------------------------------------
# Author: Ivan Horner (ivan.horner@irstea.fr), Irstea
# -----------------------------------------------------------------------------
# Last modified: 2020-02-27
# -----------------------------------------------------------------------------
###############################################################################

###############################################################################
# Required packages
###############################################################################

# cat("Loading required packages ..."); flush.console()

# suppressPackageStartupMessages({

# failed = !c(
# "devetools" = require(devtools), # not shure I need that one (it is loaded anyway with Rcpp)
# "dplyr" =     require(dplyr),
# "RcppRoll" =  require(RcppRoll),
# "segmented" = require(segmented),
# "Rcpp" =      require(Rcpp))

# })

# if (any(failed)) {
#   n = sum(failed)
#   if (n > 1) {
#     msg = paste0("The following packages are required: '", paste(names(failed)[failed], collapse = "', '"), "'")
#   } else {
#     msg = paste0("Package '", names(failed)[failed], "' is required.")
#   }
#   msg = paste0(msg, "\n\n")
#   stop(msg)
# }
# cat(" DONE\n"); flush.console()
# cat("Compiling 'RccpSignatures.cpp' ..."); flush.console()
# function_folder <- getSrcDirectory(function(x) x) # trick to get directory of current script file
# sourceCpp(file.path(function_folder, "RcppSignatures.cpp"))
# cat(" DONE\n");flush.console()


###############################################################################
# Performance metrics
###############################################################################

###############################################################################
# Runoff coefficient
###############################################################################
# DONE
###############################################################################
# Flow duration curve
###############################################################################
# DONE
###############################################################################
# Baseflow index
###############################################################################

###############################################################################
# Seasonal dynamic
###############################################################################

###############################################################################
# Recessions
###############################################################################

###############################################################################
# T-Q slopes and max Q timing
###############################################################################

###############################################################################
# MCT snow storage estimation
###############################################################################

###############################################################################
# P-Q snow storage estimation
###############################################################################
