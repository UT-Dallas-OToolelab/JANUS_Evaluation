######################################################################################
# (load.br.matrix)
# Description: this function loads a matrix and mask, with optional parameters
# for associated gallery and probe .csv files, and outputs a list of pertinent
# data and metadata for convenience in further analysis.
#
# Input: matrix filename, mask filename, gallery csv filename (optional),
# probe csv filename (optional)
#
# Output: 6 item list:
# 1) Matrix data (numeric vector of similarity/distance scores)
# 2) Mask data (numeric vector: -1 == match, 127 == non-match, 0 == ignored)
# 3) Matrix dimensions (numeric vector: rows x columns)
# 4) Distance (logical value: TRUE for distance scores, FALSE for simmilarity scores)
# 5) Protocol (2 item list: 1 - gallery csv data, 2 - probe csv data)
# 6) Filename (4 item list: 1 - matrix filename, 2 - mask filename, 3 - gallery
#    csv filename, 4 - probe csv filename)
# 
######################################################################################

load.br.matrix <- function(matrix, mask, gal = FALSE, probe = FALSE) {
  
  #open matrix file
  to.read <- file(matrix, "rb")
  #read header
  header <- readLines(to.read, n = 4)
  
  ###Extract data from header###
  #get score type: "S" for "Similarity", "D" for "Distance"
  if (substr(header[1],1,1) == "D") distance <- TRUE
  else distance <- FALSE
  #check for proper header format: if not, display warning and continue
  if (substr(header[1],2,2) != "2") message("Warning: improperly formatted matrix header")
  #get matrix dimensions
  dim <- as.numeric(strsplit(header[4], split = " ")[[1]][2:3])
  
  #read matrix data
  matrix.data <- readBin(to.read, numeric(), n = (dim[1]*dim[2]), size = 4, endian = "little")
  #close connection to matrix file
  close(to.read)
  
  #open mask file
  to.read <- file(mask, "rb")
  #read header
  header <- readLines(to.read, n = 4)
  #check for consistent score type (distance/simmilarity) between matrix and mask: if not, display warning and continue
  if (substr(header[1],1,1) == "D" && distance == FALSE) message("Warning: matrix/mask score types inconsistent")
  ##check for proper header format: if not, display warning and continue
  if (substr(header[1],2,2) != "2") message("Warning: improperly formatted mask header")
  #check for consistent dimensions between matrix and mask: if not, display warning and continue
  if (!identical(as.numeric(strsplit(header[4], split = " ")[[1]][2:3]), dim)) message("Warning: matrix/mask dimensions inconsistent")
  
  #read mask data
  mask.data <- as.numeric(readBin(to.read, logical(), n = (dim[1]*dim[2]), size = 1, endian = "little"))
  close(to.read)
  
  #if entered: load protocol CSVs
  if (!is.logical(gal)) {
    gal.data <- read.csv(gal)
    gal <- normalizePath(gal)
  }
  if (!is.logical(probe)) {
    probe.data <- read.csv(probe)
    probe <- normalizePath(probe)
  }
  #organize gal and probe in "protocol" list
  if (!exists("gal.data")) gal.data <- FALSE
  if (!exists("probe.data")) probe.data <- FALSE
  protocol <- list(gal = gal.data, probe = probe.data)
  
  output <- list(matrix = matrix.data,
                 mask = mask.data,
                 dimensions = dim,
                 distance = distance,
                 protocol = protocol,
                 filename = list(matrix = normalizePath(matrix),
                                  mask = normalizePath(mask),
                                  gal = gal,
                                  probe = probe))
  
  return(output)
}