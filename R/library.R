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
  # 2) Mask data (numeric vector: TRUE == match, FALSE == non-match, 0 == ignored)
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
  dim <- rev(as.numeric(strsplit(header[4], split = " ")[[1]][2:3]))
  
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
  if (!identical(rev(as.numeric(strsplit(header[4], split = " ")[[1]][2:3])), dim)) message("Warning: matrix/mask dimensions inconsistent")
  
  #read mask data
  mask.data <- as.numeric(readBin(to.read, logical(), n = (dim[1]*dim[2]), size = 1, endian = "little"))
  #convert mask data to logical vector
  mask.logical <- logical()
  mask.logical[mask.data == -1] <- TRUE
  mask.logical[mask.data == 127] <- FALSE
  mask.logical <- mask.logical[mask.data != 0]
  #convert vector to matrix
  mask.mat <- array(mask.logical, dim)
  #close connection to mask file
  close(to.read)
  
  #remove ignored values
  matrix.data <- matrix.data[mask.data != 0]
  #convert vector to matrix
  matrix.data <- array(matrix.data, dim)
  
  #if entered: load protocol CSVs
  if (!is.logical(gal)) {
    gal.data <- read.csv(gal)
    gal <- normalizePath(gal)
  } else gal.data <- gal
  if (!is.logical(probe)) {
    probe.data <- read.csv(probe)
    probe <- normalizePath(probe)
  } else probe.data <- probe
  
  if (!is.logical(gal) & !is.logical(probe)) {
    rownames(matrix.data) <- unique(gal.data$TEMPLATE_ID)
    colnames(matrix.data) <- unique(probe.data$TEMPLATE_ID)
    rownames(mask.mat) <- unique(gal.data$TEMPLATE_ID)
    colnames(mask.mat) <- unique(probe.data$TEMPLATE_ID)
  }
  
  output <- list(matrix = matrix.data,
                 mask = mask.mat,
                 dimensions = dim,
                 distance = distance,
                 gal = gal.data,
                 probe = probe.data,
                 n.failed = sum(matrix.data < -3.402*10^38),
                 filename = list(matrix = normalizePath(matrix),
                                 mask = normalizePath(mask),
                                 gal = gal,
                                 probe = probe))
  
#   if(check.templates(output, table = FALSE) > 0) {
#     message("Matrix and Protocol mismatch:")
#     message(sprintf("--Matrix %s", output$filename$matrix))
#     message(sprintf("--Probe %s", output$filename$probe))
#     
#     output$mismatch <- check.templates(output, table = TRUE)
#   }
  
  return(output)
}


get.roc.points <- function(data, n.points = 500, distance = FALSE) {
  
  ######################################################################################
  # Description: this function takes formatted data and plots the points of an ROC curve
  #
  # Input: load.mtx() output; or a data frame with a "match" column of "m" and "nm"
  # for "match" and "non-match", and a "resp" column of similarity scores.
  #
  # Output: roc points
  ######################################################################################
  
  #Allow load.mtx matrices to be input directly into this function
  if (!is.null(getElement(data, "matrix")) & !is.null(getElement(data, "mask")) & !is.null(getElement(data, "distance"))) {
    distance = data$distance
    match <- character()
    match[data$mask == TRUE] <- "m"
    match[data$mask == FALSE] <- "nm"
    
    data <- data.frame(match = match, resp = as.vector(data$matrix))
  }
  
  response <- data$resp
  
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


get.roc.subgroup <- function(data, n.points = 500, distance = FALSE) {
  
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
    
    if (sum(data$group[[i]]$mask) == 0) {
      output[[i]] <- list(name = sprintf("%s (no data)", data$group[[i]]$name), points = NULL)
      message(sprintf("No items in '%s' mask", data$group[[i]]$name))
    } else {
      
      #apply the mask to the data
      group.matrix <- data$matrix[data$group[[i]]$mask]
      group.mask <- data$mask[data$group[[i]]$mask]
      
      #convert masked data to get.roc.points() format (data frame with "match" and "resp" columns)
      match <- character()
      match[group.mask == TRUE] <- "m"
      match[group.mask == FALSE] <- "nm"
      group.data <- data.frame(match = match, resp = as.vector(group.matrix))
      
      #output a list with the provided subgroup title and roc points (this repeats for each mask in data$group)
      output[[i]] <- list(name = data$group[[i]]$name,
                          points = get.roc.points(group.data, n.points = n.points,
                                                  distance = data$distance),
                          n = data$group[[i]]$n)
    }
  } #end of loop
  
  return(output)
}


plot.roc <- function(data, add = FALSE, n.points = 500, col = NULL, legend = TRUE, name = "data", ...) {
  
  ######################################################################################
  # Description: this function uses the plot function with defaults optimized for ROC
  # curve plotting. Set "add = TRUE" to add a new line to an existing graph.
  #
  # Input: load.br.matrix() or get.roc.points() output
  #
  # Output: roc graph
  ######################################################################################
  
  if(!is.null(getElement(data, "group"))) {
    data <- get.roc.subgroup(data, n.points = n.points, distance = data$distance)
  } else if (!is.null(getElement(data, "matrix")) & !is.null(getElement(data, "mask")) & !is.null(getElement(data, "distance"))) {
    data <- list(list(name = name, points = get.roc.points(data, n.points = n.points, distance = data$distance)))
  } else if (class(data) == "data.frame") {
    data <- list(list(name = name, points = data))
  }
  
  args <- list(...)
  
  if (is.null(getElement(args, "type"))) args$type <- "l"
  if (is.null(col)) col <- rainbow(length(data))
  
  for (i in 1:length(data)) {
    
    if (!is.null(data[[i]]$points)) {
      args$x = data[[i]]$points
      args$col <- col[i]
      
      if (add == TRUE) do.call(lines, args)
      else do.call(plot, args)
      
      if (!is.null(getElement(args, "log"))) args$log <- NULL
      add <- TRUE
    }
  }
  
  if(legend == TRUE) {
    legend(x = "bottomright", legend = sub.elements(data, "name"), fill = col)
  }
  
  return(list(names = unlist(sub.elements(data, "name")), color = col, n = unlist(sub.elements(data, "n"))))
}


load.matrices<-function(shared.drive, algorithm.name, track, split, protocol.folder = sprintf("%s/CS0/protocol/",shared.drive)) {
  
  ######################################################################################
  # Description: Get multiple matrices
  #
  # Input: Shared drive folder, algorithm name, vector of tracks, and vector of splits
  #
  # Output: list of load.mtx() matrices
  ######################################################################################
  
  output <- list()
  i <- 1
  
  for (trk in toupper(track)) {
    for (spl in as.character(split)) {
      
      if (algorithm.name == "COTS1") {
        
        if (protocol.folder == sprintf("%s/CS0/protocol/",shared.drive)) {
          gal <- FALSE
          probe <- FALSE
        } else {
          gal <- sprintf("%s/split%s/test_%s_%s_gal.csv", protocol.folder, spl, spl, trk)
          probe <- sprintf("%s/split%s/test_%s_%s_probe.csv", protocol.folder, spl, spl, trk)
        }
        
        output[[i]] <- load.mtx(sprintf("%s/CS0/benchmarks/COTS1/split%s/verify_%s_%s.mtx", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/benchmarks/COTS1/split%s/verify_%s_%s.mask", shared.drive, spl, spl, trk),
                                gal,
                                probe)
        i <- i + 1
        
      } else if (algorithm.name == "COTS2") {
        
        if (protocol.folder == sprintf("%s/CS0/protocol/",shared.drive)) {
          gal <- FALSE
          probe <- FALSE
        } else {
          gal <- sprintf("%s/split%s/test_%s_%s_gal.csv", protocol.folder, spl, spl, trk)
          probe <- sprintf("%s/split%s/test_%s_%s_probe.csv", protocol.folder, spl, spl, trk)
        }
        
        output[[i]] <- load.mtx(sprintf("%s/CS0/benchmarks/COTS2/split%s/verify_%s_%s.mtx", shared.drive, spl, spl, trk),
                                sprintf("%s/CS0/benchmarks/COTS2/split%s/verify_%s_%s.mask", shared.drive, spl, spl, trk),
                                gal,
                                probe)
        i <- i + 1
        
      } else if (algorithm.name == "stereo") {
        
        output[[i]] <- load.mtx(sprintf("%s/stereo/stereo-results-%s-%s/scores.mtx", shared.drive, spl, tolower(trk)),
                                sprintf("%s/stereo/stereo-results-%s-%s/mask.mask", shared.drive, spl, tolower(trk)),
                                sprintf("%s/split%s/test_%s_%s_gal.csv", protocol.folder, spl, spl, trk),
                                sprintf("%s/split%s/test_%s_%s_probe.csv", protocol.folder, spl, spl, trk))
        i <- i + 1
        
      } else if (algorithm.name == "umdfv") {
        
        output[[i]] <- load.mtx(sprintf("%s/umdfv/UMD-20150102-001-%s-%s.mtx", shared.drive, trk, spl),
                                sprintf("%s/umdfv/mask_%s%s.mask", shared.drive, spl, trk),
                                sprintf("%s/split%s/test_%s_%s_gal.csv", protocol.folder, spl, spl, trk),
                                sprintf("%s/split%s/test_%s_%s_probe.csv", protocol.folder, spl, spl, trk))
        i <- i + 1
        
      } else if (algorithm.name == "umdfv2") {
        
        output[[i]] <- load.mtx(sprintf("%s/umdfv/jan-2-2015-fv-filtered-nolfwOverlap-deliverable/UMD-20150102-002-%s-%s.mtx", shared.drive, trk, spl),
                                sprintf("%s/umdfv/jan-2-2015-fv-filtered-nolfwOverlap-deliverable/UMD-20150102-002-%s-%s.mask", shared.drive, trk, spl),
                                sprintf("%s/split%s/test_%s_%s_gal.csv", protocol.folder, spl, spl, trk),
                                sprintf("%s/split%s/test_%s_%s_probe.csv", protocol.folder, spl, spl, trk))
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
  
  #If there are subgroups in the matrices
  if(!is.null(getElement(data[[1]], "group"))) {
    
    #Convert raw scores to ROC points
    points <- lapply(data, get.roc.subgroup, n.points = n.points)
    
    n.subgroups <- length(points[[1]])
    n.matrices <- length(points)
    
    FAR <- list() ; HR <- list()
    #for every subgroup
    for (i in 1:n.subgroups) {
      FAR[[i]] <- numeric()
      HR[[i]] <- numeric()
    }
    
    #for every matrix
    for (i in 1:n.matrices) {
      #for every subgroup
      for (j in 1:n.subgroups) {
        FAR[[j]] <- c(FAR[[j]], as.vector(points[[i]][[j]]$points$False.Alarms))
        HR[[j]] <- c(HR[[j]], as.vector(points[[i]][[j]]$points$Hits))
      }
    }
    
    #for every subgroup
    for (i in 1:n.subgroups) {
      #turn the vector into a matrix
      
      ###THIS LINE WORKS IN THE CONSOLE###
      #array(test[[1]][[1]], c(23,5))
      
      FAR[[i]] <- rowMeans(array(FAR[[i]], c(n.points+3, n.matrices)))
      HR[[i]] <- rowMeans(array(HR[[i]], c(n.points+3, n.matrices)))
    }
    
    #return(list(False.Alarms = data.frame(), FAR = FAR))
    
    output <- list()
    
    for (i in 1:n.subgroups) {
      output[[i]] <- list(name = points[[1]][[i]]$name,
                          points = data.frame(False.Alarms = FAR[[i]],
                                              Hits = HR[[i]]))
    }
    
    #if there are not subgroups in the matrices
  } else {
    #Convert raw scores to ROC points
    data <- lapply(data, get.roc.points, n.points = n.points)
    
    
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
  }
  
  return(output)
}


add.mask <- function(matrix, argument, mask.name, probe = TRUE, meta.data = NULL, type = "containing") {
  
  #################################################################################################
  # Description: adds a sub-grouping mask to a load.mtx() matrix
  # 
  # Input: load.mtx() matrix, data frame from protocol csv, logical argument for subgrouping, name
  # of sub-group (for legend), protocol type (FALSE for gallery, TRUE for probe)
  #
  # Output: load.mtx() matrix + subgroup mask 
  #################################################################################################
  
  if (is.null(meta.data) & probe == TRUE) {
    meta.data <- matrix$probe
  } else if (is.null(meta.data) & probe == FALSE) {
    meta.data <- matrix$gal
  }
  
  #generate mask from given argument
  if (type == "only") {
    if (probe == TRUE) {
      output <- colnames(matrix$matrix) %in% meta.data$TEMPLATE_ID[eval(parse(text = argument))]
      output <- output & !(colnames(matrix$matrix) %in% meta.data$TEMPLATE_ID[!eval(parse(text = argument))])
    } else {
      output <- rownames(matrix$matrix) %in% meta.data$TEMPLATE_ID[eval(parse(text = argument))]
      output <- output & !(rownames(matrix$matrix) %in% meta.data$TEMPLATE_ID[!eval(parse(text = argument))])
    }
  } else {
    if (probe == TRUE) {
      output <- colnames(matrix$matrix) %in% meta.data$TEMPLATE_ID[eval(parse(text = argument))]
    } else {
      output <- rownames(matrix$matrix) %in% meta.data$TEMPLATE_ID[eval(parse(text = argument))]
    }
  }
  
  #names(output) <- unique(meta.data$TEMPLATE_ID)
  if (probe == TRUE) {
    output <- array(rep(output, each = matrix$dimensions[1]), c(matrix$dimensions[1],matrix$dimensions[2]))
  } else {
    output <- array(rep(output, times = matrix$dimensions[1]), c(matrix$dimensions[1],matrix$dimensions[2]))
  }
  
  mask <- list(mask = output, name = mask.name, probe = probe, n = sum(output))
  
  if (mask$n == 0) message(sprintf("WARNING: mask '%s' contains no values", mask$name))
  
  if (!is.null(getElement(matrix, "group"))) {
    matrix$group[[length(matrix$group) + 1]] <- mask
  } else {
    matrix$group <- list(mask)
  }
  
  return(matrix)
}


batch.add.mask <- function(matrices, argument, mask.name, probe = TRUE, containing = TRUE) {
  for (i in 1:length(matrices)) {
    output <- add.mask(matrices[[i]], argument = argument, mask.name = mask.name, probe = probe, containing = containing)
  }
}


combine.masks <- function(matrix, mask.1, mask.2, name, argument = "and", clear.old = FALSE) {
  
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
  
  if(argument == "and") {
    new.mask <- mask.1.data & mask.2.data
  } else if (argument == "or") {
    new.mask <- mask.1.data | mask.2.data
  } else if (argument == "xor") {
    new.mask <- xor(mask.1.data, mask.2.data)
  } else if (argument == "nor") {
    new.mask <- !(mask.1.data | mask.2.data)
  } else {
    message("Argument not recognized! Recognized arguments are: and, or, xor, nor")
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
                                                   probe = probe,
                                                   n = sum(new.mask))
  
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


add.img.frame.masks <- function(matrix, probe = TRUE, meta.data = NULL, img = TRUE, frame = TRUE, both = TRUE) {
  
  #################################################################################################
  # Description: Adds a mask for 1.) images only, 2.) video only, and 3.) both images and video.
  # These masks are important/complicated enough to warrant their own function.
  # 
  # Input: load.mtx() matrix, gallery or probe data frame from csv, is the meta-data probe or gallery? (T/F)
  #
  # Output: load.mtx() matrices + subgroup masks
  #################################################################################################
  
  if (is.null(meta.data) & probe == TRUE) {
    meta.data <- matrix$probe
  } else if (is.null(meta.data) & probe == FALSE) {
    meta.data <- matrix$gal
  }
  
  #Get "image and video" group
  if (frame == TRUE) matrix <- add.mask(matrix, "grepl('frame', meta.data$FILE)", "Only video", probe = probe, type = "only")
  if (img == TRUE) matrix <- add.mask(matrix, "grepl('img', meta.data$FILE)", "Only images", probe = probe, type = "only")
  if (both == TRUE) matrix <- combine.masks(matrix, length(matrix$group), length(matrix$group) - 1, name = "Images and video",
                                            argument = "nor", clear.old = FALSE)
  
  return(matrix)
}


check.templates <- function(matrix, table = TRUE) {
  
  mat <- matrix
  
  matches.per.row <- rowSums(matrix(mat$mask, mat$dimensions[2], mat$dimensions[1]) == TRUE)
  
  probe.temp.per.sub <- mat$probe[,1:2]
  probe.temp.per.sub <- probe.temp.per.sub[!duplicated(probe.temp.per.sub[,1]),]
  probe.temp.per.sub <- as.vector(table(factor(probe.temp.per.sub[,2], levels=unique(probe.temp.per.sub[,2]))))
  
  print(matches.per.row == probe.temp.per.sub)
  print(matches.per.row)
  print(probe.temp.per.sub)
  
  if (table == TRUE) {
    return(data.frame(from.mtx = matches.per.row, from.csv = probe.temp.per.sub, equal = matches.per.row == probe.temp.per.sub))
  } else {
    return(sum(matches.per.row != probe.temp.per.sub))
  }
}


batch.img.frame.masks <- function(matrices, probe = TRUE, meta.data = NULL) {
  
  for (i in 1:length(matrices)) {
    
    if (is.null(meta.data) & probe == TRUE) {
      protocol <- matrices[[i]]$probe
    } else if (is.null(meta.data) & probe == FALSE) {
      protocol <- matrices[[i]]$gal
    } else {
      protocol <- meta.data
    }
    
    matrices[[i]] <- add.img.frame.masks(matrices[[i]], meta.data = protocol, probe = probe)
  }
  
  return(matrices)
}

visualize <- function(data, square = TRUE, ...) {
  
  #################################################################################################
  # Description: Visualize data in a matrix table without the weird defaults of heatmap() and image()
  # 
  # Input: Matrix (as in an actual n x m matrix of data)
  #
  # Output: "heat map" of values in the matrix
  #################################################################################################
  
  args <- list(...)
  if (is.null(getElement(args, "asp")) & square == TRUE) args$asp <- nrow(data)/ncol(data)
  data <- t(data[nrow(data):1,])
  args$x <- data
  do.call(image, args)
}


plot.cmc <- function(matrix, n.points = nrow(matrix$matrix), plot = TRUE, add = FALSE, cutoff = NULL, cutoff.col = "blue", ...) {
  
  #################################################################################################
  # Description: Plots a Cumulative Match Characteristic (CMC) curve from a matrix or precalculated points
  # 
  # Input: load.mtx() matrix; maximum rank to plot; plot graph, or save points?; add to existing plot?
  #
  # Output: CMC graph, or CMC points, depending on value of "plot" argument
  #################################################################################################
  
  #get extra user input
  args <- list(...)
  
  #if the input is a data frame, assume it contains cmc points and skip calculating them
  if(class(matrix) == "data.frame") cmc <- matrix
  else {
    
    #initialize variables for loop
    match.ranks <- numeric()
    rank <- numeric()
    rr <- numeric()
    
    #for every column, get the rank (in descending order) of that column's matched pair
    for (i in 1:ncol(matrix$matrix)) {
      #If the scores are similarity scores, flip the sign to reverse the rank order
      if (matrix$distance == FALSE) col.values <- -matrix$matrix[,i]
      else col.values <- matrix$matrix[,i]
      #get the rank of this column's matched pair
      match.ranks[i] <- rank(col.values)[matrix$mask[,i] == TRUE]
    }
    #for every point to be plotted (defaults to number of rows -- the highest possible rank)
    for (i in 1:n.points) {
      rank[i] <- i
      rr[i] <- sum(match.ranks <= i)
    }
    rr <- rr/max(rr)
    cmc <- data.frame(Rank = rank, Retrieval.Rate = rr)
  }
  
  if (plot == TRUE) {
    if (is.null(getElement(args, "type"))) args$type <- "l"
    args$x <- cmc
    
    if (add == FALSE) do.call("plot", args)
    else do.call("lines", args)
    
    if (!is.null(cutoff)) {
      for (i in 1:length(cutoff)) {
        
        lines(cbind(c(cutoff[i],cutoff[i]), c(0, cmc[cutoff[i],2])), lty = 2, col = cutoff.col)
        lines(cbind(c(-10,cutoff[i]), c(cmc[cutoff[i],2], cmc[cutoff[i],2])), lty = 2, col = cutoff.col)
      }
    }
  }
  else return(cmc)
}

mask.for.MEDIA_IDs <- function(old.data, split, components, operator){
  new.data = old.data
  
  full.data <- numeric()
  
  for(spl in split){
    media = paste(new.data[[spl]]$probe$TEMPLATE_ID, new.data[[spl]]$probe$MEDIA_ID, sep = ".")
    probe.medid.per.temp = data.frame(media = media, template = new.data[[spl]]$probe$TEMPLATE_ID)
    probe.medid.per.temp <- probe.medid.per.temp[!duplicated(probe.medid.per.temp[,1]),]
    probe.medid.per.temp <- table(factor(probe.medid.per.temp[,2], levels=unique(probe.medid.per.temp[,2])))
    
    full.data = c(full.data, probe.medid.per.temp)
    
    probe.medid.per.temp.matrix = probe.medid.per.temp
    
    for(i in 1:49){
      probe.medid.per.temp.matrix = rbind(probe.medid.per.temp.matrix, probe.medid.per.temp)
    }
    
    for(com in as.character(components)){
      argument = sprintf('probe.medid.per.temp.matrix %s %s', operator, com)
      mask.name = sprintf('%s%s MEDIA_IDs', operator, com)
      new.data[[spl]] = add.mask(new.data[[spl]], argument, mask.name)
    }
  }
  return(new.data)
}