load.mtx <- function(matrix, mask, gal = FALSE, probe = FALSE) {

  ######################################################################################
  # Description: this function loads a matrix and mask, with optional parameters
  # for associated gallery and probe .csv files, and outputs a list of useful
  # data and metadata for convenience in further analysis.
  #
  # Input: matrix filename, mask filename, gallery csv filename (optional),
  # probe csv filename (optional)
  #
  # Output: 6 item list:
  # 1) Matrix similarity data (numeric vector of similarity/distance scores)
  # 2) Mask data (numeric vector: -1 == match, 127 == non-match, 0 == ignored)
  # 3) Matrix dimensions (numeric vector: rows x columns)
  # 4) Distance (logical value: TRUE for distance scores, FALSE for similarity scores)
  # 5) Protocol (2 item list: 1 - gallery csv data, 2 - probe csv data)
  # 6) Filename (4 item list: 1 - matrix filename, 2 - mask filename, 3 - gallery
  #    csv filename, 4 - probe csv filename)
  ######################################################################################
  
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


get.roc.points <- function(data, n.points = 500, distance = FALSE, response = data$resp) {
  
  ######################################################################################
  # Description: this function takes formatted data and plots the points of an ROC curve
  #
  # Input: load.mtx() output; or a data frame with a "match" column of "m" and "nm"
  # for "match" and "non-match", and a "resp" column of similarity scores.
  #
  # Output: roc points
  ######################################################################################
  
  #Allow load.br.matrix matrices to be input directly into this function
  if (!is.null(getElement(data, "matrix")) & !is.null(getElement(data, "mask")) & !is.null(getElement(data, "distance"))) {
    distance = data$distance
    match <- character()
    match[data$mask == -1] <- "m"
    match[data$mask == 127] <- "nm"
    
    data <- data.frame(match = match, resp = data$matrix)
    
    response <- data$resp
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


get.roc.subgroup <- function(data, n.points = 500, distance = FALSE, response = data$resp) {

  ######################################################################################
  # Description: this function calculates a separate ROC curve for however many masks
  # have been added to the input data matrix using the add.mask() function.
  #
  # Input: takes a load.mtx() matrix, including masks obtained from add.mask()
  #
  # Output: list of groups of roc points paired with their titles
  ######################################################################################
  
  #initialize output variable
  output <- list()
  
  #for each mask defined in the "group" list in the data matrix
  for(i in 1:length(data$group)) {
    #if the mask applies to the gallery (rows) or probe (columns)
    if (data$group[[i]]$probe == TRUE) {
      #apply the mask by converting the .mask and .mtx scores to matrix form, and removing the masked columns
      group.matrix <- as.vector(matrix(data$matrix, data$dimensions[2], data$dimensions[1])[,data$group[[i]]$mask])
      group.mask <- as.vector(matrix(data$mask, data$dimensions[2], data$dimensions[1])[,data$group[[i]]$mask])
    } else {
      #apply the mask by converting the .mask and .mtx scores to matrix form, and removing the masked rows
      group <- as.vector(matrix(data$matrix, data$dimensions[2], data$dimensions[1])[data$group[[i]]$mask,])
      group.mask <- as.vector(matrix(data$mask, data$dimensions[2], data$dimensions[1])[data$group[[i]]$mask,])
    }
    
    #convert masked data to get.roc.points() format (data frame with "match" and "resp" columns)
    distance = data$distance
    match <- character()
    match[group.mask == -1] <- "m"
    match[group.mask == 127] <- "nm"
    group.data <- data.frame(match = match, resp = group.matrix)
    
    #output a list with the provided subgroup title and roc points (this repeats for each mask in data$group)
    output[[i]] <- list(name = data$group[[i]]$name, points = get.roc.points(group.data, n.points = n.points,
                                                                             distance = distance))
  } #end of loop
  
  return(output)
}


plot.roc <- function(data, add = FALSE, ...) {

  ######################################################################################
  # Description: this function uses the plot function with defaults optimized for ROC
  # curve plotting. Set "add = TRUE" to add a new line to an existing graph.
  #
  # Input: load.br.matrix() or get.roc.points() output
  #
  # Output: roc graph
  ######################################################################################
  
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


plot.rocs <- function(roc.list, add = FALSE, legend = TRUE, col = rainbow(length(roc.list)), ...) {
  
  ######################################################################################
  # Description: this function iterates the plot.roc() function for every ROC subgroup
  # contained in the input variable you feed it, drawing a new line for each group.
  #
  # Input: load.mtx() + add.mask() output, or get.roc.subgroup() output; add to existing
  # plot? (TRUE for yes, FALSE for no); plot legend? (T/F); vector of ROC curve colors
  #
  # Output: ROC graph.
  ######################################################################################
  
  args <- list(...)
  
  if(!is.null(getElement(roc.list, "group"))) {
    roc.list <- get.roc.subgroup(roc.list)
  }
  
  plot.roc(roc.list[[1]]$points, add = add, col = col[1], args)
  
  if(length(roc.list) > 1) {
    for(i in 2:length(roc.list)) {
      plot.roc(roc.list[[i]]$points, add = TRUE, col = col[i])
    }
  }
  
  if(legend == TRUE) {
    legend(x = "bottomright", legend = sub.elements(roc.list, "name"), fill = col)
  }
  
  return(list(names = unlist(sub.elements(roc.list, "name")), color = col))
}


load.matrices<-function(shared.drive, algorithm.name, track, split) {
  
  ######################################################################################
  # Description: Get multiple matrices
  #
  # Input: Shared drive folder, algorithm name, vector of tracks, and vector of splits
  #
  # Output: list of load.mtx() matrices
  ######################################################################################
  
  split <- as.character(split)
  output <- list()
  i <- 1
  
  for (trk in track) {
    for (spl in split) {
      
      if (algorithm.name == "COTS1") {
        
        output[[i]] <- load.mtx(sprintf("%s/CS0/benchmarks/COTS1/split%s/verify_%s_%s.mtx", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/benchmarks/COTS1/split%s/verify_%s_%s.mask", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_gal.csv", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_probe.csv", shared.drive, spl, spl, trk))
        i <- i + 1
        
      } else if (algorithm.name == "COTS2") {
        
        output[[i]] <- load.mtx(sprintf("%s/CS0/benchmarks/COTS2/split%s/verify_%s_%s.mtx", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/benchmarks/COTS2/split%s/verify_%s_%s.mask", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_gal.csv", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_probe.csv", shared.drive, spl, spl, trk))
        i <- i + 1
        
      } else if (algorithm.name == "stereo") {
        
        if (spl == "6" & trk == "B") {
          message("WARNING: at the time this code was written, stereo 6 B was not in the shared drive. Skipping it to avoid errors.")
        } else {
          output[[i]] <- load.mtx(sprintf("%s/stereo/stereo-results-%s-%s/scores.mtx", shared.drive, spl, tolower(trk)),
                                  sprintf("%s/stereo/stereo-results-%s-%s/mask.mask", shared.drive, spl, tolower(trk)),
                                  sprintf("%s/CS0/protocol/split%s/test_%s_%s_gal.csv", shared.drive, spl, spl, trk),
                                  sprintf("%s/CS0/protocol/split%s/test_%s_%s_probe.csv", shared.drive, spl, spl, trk))
          i <- i + 1
        }
                              
      } else if (algorithm.name == "umdfv") {
        
        output[[i]] <- load.mtx(sprintf("%s/umdfv/UMD-20150102-001-%s-%s.mtx", shared.drive, trk, spl),
                                sprintf("%s/umdfv/mask_%s%s.mask", shared.drive, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_gal.csv", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/protocol/split%s/test_%s_%s_probe.csv", shared.drive, spl, spl, trk))
        i <- i + 1
        
      } else {
        message("The algorithm you entered is not supported. Check your spelling/capitalization and try again.")
        message("Supported algorithms: COTS1, COTS2, stereo, umdfv")
        stop()
      }
      
      
      
    }
  }
  
  if(length(output) == 1) {
    output <- output[[1]]
  }
  
  return(output)
}


mean.roc <- function(data, n.points = 500) {
  
  ######################################################################################
  # Description: This function takes a list of load.mtx() matrices and gets their mean ROC
  #
  # Input: list of load.mtx() matrices
  #
  # Output: mean ROC points data frame
  ######################################################################################
  
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


add.mask <- function(matrix, meta.data, argument, mask.name, probe = TRUE) {

  #################################################################################################
  # Description: adds a sub-grouping mask to a load.mtx() matrix
  # 
  # Input: load.mtx() matrix, data frame from protocol csv, logical argument for subgrouping, name
  # of sub-group (for legend), protocol type (FALSE for gallery, TRUE for probe)
  #
  # Output: load.mtx() matrix + subgroup mask 
  #################################################################################################
  
  #generate mask from given argument
  output <- unique(meta.data$TEMPLATE_ID) %in% meta.data$TEMPLATE_ID[eval(parse(text = argument))]
  
  mask <- list(mask = output, name = mask.name, probe = probe)
  
  #if this 
  if (!is.null(getElement(matrix, "group"))) {
    matrix$group[[length(matrix$group) + 1]] <- mask
  } else {
    matrix$group <- list(mask)
  }
  
  return(matrix)
}


intersect <- function(matrix, mask.1, mask.2, name, argument = "and", clear.old = FALSE) {
  
  #################################################################################################
  # Description: combine masks to account for the complexity within templates
  # 
  # Input: load.mtx() matrix, mask 1 index, mask 2 index, name of logical operator to be used,
  # delete masks used to make new mask? (T/F)
  #
  # Output: logical (TRUE/FALSE) mask
  #################################################################################################
  
  mask.1.data <- matrix$group[[mask.1]]$mask
  mask.2.data <- matrix$group[[mask.2]]$mask
  probe <- matrix$group[[mask.1]]$probe
  
  if (matrix$group[[mask.1]]$probe != matrix$group[[mask.2]]$probe) {
    message("Gallery masks cannot be combined with probe masks!")
    stop()
  }
  
  if(argument == "and") {
    new.mask <- mask.1.data & mask.2.data
  } else if (argument == "or") {
    new.mask <- mask.1.data | mask.2.data
  } else if (argument == "xor") {
    new.mask <- xor(mask.1.data, mask.2.data)
  } else {
    message("Argument not recognized! Recognized arguments are: and, or, xor")
    return(NULL)
  }
  
  if (clear.old == TRUE) {
    matrix$group[[mask.1]] <- NULL
    if (mask.1 > mask.2) {
      matrix$group[[mask.2]] <- NULL
    } else {
      matrix$group[[mask.2 - 1]] <- NULL
    }
  }
  
  matrix$group[[length(matrix$group) + 1]] <- list(mask = new.mask,
                                               name = name,
                                               probe = probe)
  
  return(matrix)
}


sub.elements <- function(list, sub.element) {
  
  #################################################################################################
  # Description: get a list of sub-elements from a higher-level list. For instance, if you have
  # five load.mtx() matrices and you want to grab their $matrix data only, this function can extract
  # just those requested elements from the higher-level list.
  # 
  # Input: higher-level list (of anything)
  #
  # Output: lower-level elements
  #################################################################################################
  
  output <- list()
  for(i in 1:length(list)) {
    output[[i]] <- getElement(list[[i]], sub.element)
  }
  return(output)
}


add.img.frame.masks <- function(matrix, meta.data, probe = TRUE, clear.old = TRUE) {
  
  #################################################################################################
  # Description: Adds a mask for 1.) images only, 2.) video only, and 3.) both images and video.
  # These masks are important/complicated enough to warrant their own function.
  # 
  # Input: load.mtx() matrix, gallery or probe data frame from csv, is the meta-data probe or gallery? (T/F)
  #
  # Output: load.mtx() matrices + subgroup masks
  #################################################################################################
  
  #Get "image and video" group
  matrix <- add.mask(matrix, meta.data, "grepl('img', meta.data$FILE)", "Only video", probe = probe)
  matrix <- add.mask(matrix, meta.data, "grepl('frame', meta.data$FILE)", "Only images", probe = probe)
  matrix <- intersect(matrix, length(matrix$group), length(matrix$group) - 1, name = "Images and video", argument = "and", clear.old = clear.old)
  
  matrix$group[[length(matrix$group) - 2]]$mask <- !matrix$group[[length(matrix$group) - 2]]$mask
  matrix$group[[length(matrix$group) - 1]]$mask <- !matrix$group[[length(matrix$group) - 1]]$mask
  
  return(matrix)
}