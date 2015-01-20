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
  
  #document ignored values
  ignored <- list(n = sum(mask.data == 0),
                  values = matrix.data[mask.data == 0])
  
  #remove ignored values from output
  matrix.data <- matrix.data[mask.data != 0]
  mask.data <- mask.data[mask.data != 0]
  
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
                 n.failed = sum(matrix.data < -3.402*10^38),
                 ignored = ignored,
                 filename = list(matrix = normalizePath(matrix),
                                  mask = normalizePath(mask),
                                  gal = gal,
                                  probe = probe))
  
  return(output)
}

######################################################################################
# (get.roc.points)
# Description: this function takes formatted data and plots the points of an ROC curve
#
# Input: load.br.matrix output
#
# Output: roc points
# 
######################################################################################

get.roc.points <- function(data, n.points = 500, distance = FALSE, response = data$resp) {
  
  #Allow load.br.matrix matrices to be input directly into this function
  if (!is.null(getElement(data, "matrix")) & !is.null(getElement(data, "mask")) & !is.null(getElement(data, "distance"))) {
    distance = data$distance
    match <- character()
    match[data$mask == -1] <- "m"
    match[data$mask == 127] <- "nm"
    
    data <- data.frame(match = match, resp = data$matrix)
  }
  
  #Define a variable just below the (functionally) Inf value to limit the range to actual values
  pre.inf <- 3.402*10^38
  
  #Omit -Inf and +Inf from the range calculation
  range <- range(response[response > (-pre.inf) & response < (pre.inf)])
  
  #Calculate the cutoff scores for each ROC point
  thresh <- seq(min(range),max(range),by=diff(range)/(n.points))
  
  #Initialize variables to count False Alarms and Hits at each threshold
  fa.tally <- numeric()
  hit.tally <- numeric()
  
  #Loop through each threshold value, counting number of values either below or above the threshold
  #depending on the value of "distance"
  for (i in 1:(n.points+1)) {
    if (distance == TRUE) {
      fa.tally[i] <- sum(data$match =='nm'&response<=thresh[i])
      hit.tally[i] <- sum(data$match =='m'&response<=thresh[i])
    } else {
      fa.tally[i] <- sum(data$match =='nm'&response>=thresh[i])
      hit.tally[i] <- sum(data$match =='m'&response>=thresh[i])
    }
  }
  
  #Get number of match/non-match pairs
  n.match <- sum(data$match =='m')
  n.nonmatch <- sum(data$match =='nm')
  
  #Divide tallies by number of match/non-match to get percent scores
  if(distance == TRUE) {
    far <- c(0,(fa.tally/n.nonmatch),1) #add (0, 0) and (1, 1) points for curve plotting
    hr <- c(0,(hit.tally/n.match),1)
  } else {
    far <- c(1,(fa.tally/n.nonmatch),0)
    hr <- c(1,(hit.tally/n.match),0)
  }
  
  #Structure output for plotting
  ROC <- data.frame(False.Alarms = far, Hits = hr)
  
  return(ROC)
}

######################################################################################
# (plot.roc)
# Description: this function uses the plot function with defaults optimized for ROC
# curve plotting. Set "add = TRUE" to add a new line to an existing graph.
#
# Input: load.br.matrix() or get.roc.points() output
#
# Output: roc graph
# 
######################################################################################

plot.roc <- function(data, add = FALSE, ...) {
  
  if (!is.null(getElement(data, "matrix")) & !is.null(getElement(data, "mask")) & !is.null(getElement(data, "distance"))) {
    data <- get.roc.points(data, distance = data$distance)
  }
  
  args <- list(...)
  
  data[data == 0] <- 1*10^-50
  
  args <- c(args, list(x = data))
  if (is.null(getElement(args, "type"))) args$type <- "l"
  if (is.null(getElement(args, "col"))) args$col <- "blue"
  
  if (add == TRUE) do.call(lines, args)
  else do.call(plot, args)
}


######################################################################################
# Get multiple matrices across all scripts with corresponding masks and metadata
#
#TO USE: setwd to whatever challenge set you're using. eg: setwd('~/JANUS_Drive/CS0/')
#
#input: designate which alg you're examining (COTS) as a string
#output: nested list with all matrices, masks, gal, and probe csvs from all splits
######################################################################################

load.multiple.matrices<-function(COTS) {
######################################################################################
  #AH 
  #Get multiple matrices across all scripts with corresponding masks and metadata
  #
  #TO USE: setwd to whatever challenge set you're using. eg: setwd('~/JANUS_Drive/CS0/')
  #
  #input: desigate whith alg you're examining (COTS) as a string
  #output: nested list with all matrices, masks, gal, and probe csvs from all splits
######################################################################################
  
## Get data ##
  
  #initialize lists
  all.splits <- list(A = list(split = list()),
                     B = list(split = list()))
  
  #get data for given algorithm
  cots.name <- COTS
  for (i in 1:10){ 
    test.name <- c(sprintf('%s_A',i),
                   sprintf('%s_B',i))
    
    #Initialize/clear strings for pathnames of matrices and metadata
    mat <- character() ; mask <- character() ; gal <- character() ; probe <- character()
    #Iterate pathnames
    for (e in 1:2) {
      mat <- c(mat, sprintf('./benchmarks/%s/split%s/verify_%s.mtx', cots.name, i, test.name[e]))
      mask <- c(mask, sprintf('./benchmarks/%s/split%s/verify_%s.mask', cots.name, i, test.name[e]))
      gal <- c(gal, sprintf('./protocol/split%s/test_%s_gal.csv', i, test.name[e]))
      probe <- c(probe, sprintf('./protocol/split%s/test_%s_probe.csv', i, test.name[e]))
    }
    
    
    all.splits$A$split[[i]] <- load.br.matrix(mat[1],mask[1],gal[1],probe[1])
    all.splits$B$split[[i]] <- load.br.matrix(mat[2],mask[2],gal[2],probe[2])
    
  }
  
  return(all.splits)
}


######################################################################################
# (plot.roc)
# Description: This function takes lists of loaded matrices and gets a mean ROC of them
#
# Input: list of load.br.matrix() output
#
# Output: mean roc points
# 
######################################################################################

#This function takes lists of matrices and gets a mean ROC of them
mean.roc <- function(data, n.points = 500) {
  
  #Convert raw scores to ROC points
  data <- lapply(data, get.roc.points, n.points = n.points)
  
  #output <- data.frame(False.Alarms = numeric(), Hits = numeric())
  #Get means of returned data
  
  FAR <- matrix(nrow = length(data), ncol = length(data[[1]]$False.Alarms)) ; HR <- FAR
  
  #data is the list of get.roc.points output data frames
  #data[[e]] is the current data frame
  #data[[e]][[i]] is the current point in the current data frame
  for (i in 1:ncol(FAR)) {
    for (e in 1:nrow(FAR)) {
      FAR[e,i] <- data[[e]]$False.Alarms[i]
      HR[e,i] <- data[[e]]$Hits[i]
    }
  }
  
  output <- data.frame(False.Alarms = colMeans(FAR), Hits = colMeans(HR))
  return(output)
}
