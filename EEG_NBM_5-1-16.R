# Required packages.
library("ggplot2")
library("stats")
library("reshape2")
library("plyr")
#library("plotly")

#
# Functions for data extraction
#

# Apply the fast Fourier Transform to each electrode for each tester (subject).
ApplyFFT <- function(df, subj.ids) {
  # Build the initial data frame...
  # fft.frame stores all FFT data in a 256 x (34*25) data frame.
  fft.frame <- data.frame(id = matrix(nrow = 0, ncol = 1),
                          tp = matrix(nrow = 0, ncol = 1),
                          matrix(nrow = 0, ncol = 25))
  
  # Loop through all the data in combined table.
  for (i in 1:length(subj.ids)) {
    # Declare an empty matrix to store FFT results.
    fft.results <- matrix(nrow = 256, ncol = 0)
    
    # Now, loop through the electrodes contained in columns [4,28].
    for (j in 4:length(df)) {
      # Create a matrix of FFT values for 1 electrode, 1 subject.
      # Get the results of FFT.
      fft.results <- cbind(fft.results, fft(df[df[,2] == subject.ids[i],j]))
    }
    
    # Populate the table.
    fft.table <- data.frame(id = rep(subj.ids[i],256),
                            tp = 1:256,
                            fft.results)
    # Correct the electrode column names.
    colnames(fft.table) <- c("id", "tp", seq(1,25))
    
    # Store the frequencies in the fft.frame.
    fft.frame <- rbind(fft.frame, fft.table)
  }
  
  fft.frame # Return the fft.frame.
}

# Calculate the magnitude for all members of an fft dataset.
ApplyMagnitude <- function(df){
  data.frame(df[,1], abs(df[,2:length(df)]))
}

# Extract data for specific subject ID sets.  Returns a data.frame 
# containing only data with for the list of subjects passed to the 
# function.
ExtractSubjectData <- function(df, subj.ids){
  df.subset <- data.frame(matrix(nrow = 0, ncol = 27))
  for(i in 1:length(subj.ids)){
    df.tmp <- data.frame(df[df[,1] == subj.ids[i],])
    df.subset <- rbind(df.subset,df.tmp)
  }
  df.subset
}

# Obtain mean data for each electrode in the list and data.frame passed 
# to the function. Returns a n x 26 data.frame, with subject IDs in 
# column 1.
ExtractMeans <- function(df) {
  df.subset <- data.frame(matrix(NA, nrow = 1, ncol = 25))
  #print(df.subset)
  df.subset[1:25] <- colMeans(df[,3:27])
  y <- data.frame(el=seq(1:25))
  y <- cbind(y, mag.mean=matrix(df.subset[1,1:25], nrow=25, ncol=1))
  y
}

# Reunite the delta values into a single df. First frame should be "happy", 2nd "neutral"
UniteFrame <- function(h.df, n.df, title) {
  cond <- c("")
  df.united <- data.frame(Condition=rep(c("Happy","Neutral"), times=1, each=25))
  df.united <- cbind(df.united, rbind(h.df,n.df))
  df.united <- cbind(df.united,h.df)
  df.united <- cbind(df.united, g.mag.mean=unlist(n.df[,2]))
  df.united
}

# Accept a data.frame with the layout creawted by UniteFrame() and 
# return a list of plots.
CreatePlots <- function(df, title){
  # Unlist the magnitude means
  df$mag.mean <- unlist(df$mag.mean)
  
  # Set names
  names(df)[2:3] <- c("Electrode", "Magnitude")
  
  # Create the plot
  tmp.list <- ggplot(data=df, 
                     aes(x=Electrode, 
                         y=Magnitude, 
                         group = Condition, 
                         colour = Condition)) +
    geom_line() + 
    ggtitle(title)
  tmp.list
}

# Spectral density function. Takes a subject's magnitude frame and a frequency
# band (frequency range, example: c(1, 6.5), which would indicate 1-6.5 Hz) and 
# returns a subject data frame populatd with the spectral densities for each 
# electrode.
ObtainDensity <- function(df, freq_band) {
  # Construct empty data frame.
  density.frame <- data.frame(matrix(nrow = 1, ncol = 26))
  density.frame[1] <- df[1,1] 
  colnames(density.frame) <- c("id", seq(1,25))
  
  # Calculate total length of frequency band for dividing integral by later.
  a <- freq_band[1]
  b <- freq_band[2]
  length <- b - a
  
  # Loop through the discrete time range.
  for (i in 3:27) {
    value <- 0
    
    # Loop through all electrodes.
    for (j in a:b) {
      value <- value + df[j,i]**2
    }
    
    # Finish calculation of spectral density.
    value <- value / length
    
    # Ultimately, new subject data frame will contain a total of:
    # 26 columns (first is id, 2-26 are electrodes 1-25).
    # 1 appended row of the specified frequency band.
    density.frame[1,i-1] <- value
  }
  
  # Return...
  density.frame
}

FindMinMax <- function(df) {
  # Build empty matrix.
  values <- matrix(nrow = 2, ncol = 4)
  colnames(values) <- c("delta", "theta", "alpha", "beta")
  rownames(values) <- c("min", "max")
  
  # Assign min values for each band.
  values[1,1] <- min(df[df[,1] == "delta", 3:27])  
  values[1,2] <- min(df[df[,1] == "theta", 3:27])  
  values[1,3] <- min(df[df[,1] == "alpha", 3:27])  
  values[1,4] <- min(df[df[,1] == "beta", 3:27])  
  
  # Assign max values for each band.
  values[2,1] <- max(df[df[,1] == "delta", 3:27])  
  values[2,2] <- max(df[df[,1] == "theta", 3:27])  
  values[2,3] <- max(df[df[,1] == "alpha", 3:27])  
  values[2,4] <- max(df[df[,1] == "beta", 3:27])  
  
  # Return the matrix.
  values
}

# GetAllDensities function accepts a data.frame and list of subject IDs.
# The function calls ObtainDensity function for each frequency band and
# consolidates them in a single data.frame.  The data.frame is returned.
GetAllDensities <- function(df, subj.ids){
  # ApplyMagnitude to df.
  df <- ApplyMagnitude(df)
  
  # Obtain the spectral density for the delta/theta/alpha/beta
  spec.den <- data.frame(matrix(nrow = 0, ncol = 27))
  
  # Enter the names for each band to be requested
  colnames(spec.den) <- (c("freq.band", "id", seq(1,25)))
  
  for(i in 1:length(subj.ids)){
    # Temporary frames to collect densities for 5 frequency bands.
    df.tmp <- data.frame(df[df[,1] == subj.ids[i],])
    d.block <- data.frame(matrix(nrow = 4, ncol = 1))
    d.block[,1] <- c("delta", "theta", "alpha", "beta")
    colnames(d.block) <- c("freq.band")
    d.vals <- data.frame(matrix(nrow = 0, ncol = 26))
    
    # Delta
    tmp <- ObtainDensity(df.tmp, c(1,7))
    d.vals <- rbind(d.vals, tmp)
    # Theta
    tmp <- ObtainDensity(df.tmp, c(8,9))
    d.vals <- rbind(d.vals, tmp)
    # Alpha
    tmp <- ObtainDensity(df.tmp, c(10, 12))
    d.vals <- rbind(d.vals, tmp)
    # Beta
    tmp <- ObtainDensity(df.tmp, c(13, 30))
    d.vals <- rbind(d.vals, tmp)
    
    # Add to spec.den
    d.block <- cbind(d.block, d.vals)
    spec.den <- rbind(spec.den, d.block)
  }
  
  # Return...
  spec.den
}

# Create bins (100 bins per band, per subject type)
CreateBins <- function(df.min.max, df.spec.dens, bin.total){
  # Create a data.frame to contain all the bins
  bin.frame <- data.frame(matrix(0, nrow = 4, ncol = bin.total+1))
  colnames(bin.frame) <- c("freq.band", seq(1,100))
  bin.frame[1,1] <- "delta"
  bin.frame[2,1] <- "theta"
  bin.frame[3,1] <- "alpha"
  bin.frame[4,1] <- "beta"
  
  # Calculate partition widths.
  delta.max <- df.min.max[2,1]
  delta.min <- df.min.max[1,1]
  delta.width <- (delta.max - delta.min) / bin.total
  
  theta.max <- df.min.max[2,2]
  theta.min <- df.min.max[1,2]
  theta.width <- (theta.max - theta.min) / bin.total
  
  alpha.max <- df.min.max[2,3]
  alpha.min <- df.min.max[1,3]
  alpha.width <- (alpha.max - alpha.min) / bin.total
  
  beta.max <- df.min.max[2,4]
  beta.min <- df.min.max[1,4]
  beta.width <- (beta.max - beta.min) / bin.total
  
  # Bin limit.
  bin.limit <- data.frame(matrix(0, nrow = 4, ncol = bin.total))
  
  # Iterate through bin.total to find bin frequencies.
  for (i in 1:bin.total) {
    # Iterate through df.spec.dens. 
    for (j in 1:nrow(df.spec.dens)) {
      # Calculate current bin limit. Current bin is i.
      if (df.spec.dens[j,1] == "delta") {
        a <- delta.min + ((i-1) * delta.width)
        b <- delta.min + (i * delta.width)
        band <- "delta"
        bin.limit[1,i] <- b
      }
      else if (df.spec.dens[j,1] == "theta") {
        a <- theta.min + ((i-1) * theta.width)
        b <- theta.min + (i * theta.width)
        band <- "theta"
        bin.limit[2,i] <- b
      }
      else if (df.spec.dens[j,1] == "alpha") {
        a <- alpha.min + ((i-1) * alpha.width)
        b <- alpha.min + (i * alpha.width)
        band <- "alpha"
        bin.limit[2,i] <- b
      }    
      else if (df.spec.dens[j,1] == "beta") {
        a <- beta.min + ((i-1) * beta.width)
        b <- beta.min + (i * beta.width)
        band <- "beta"
        bin.limit[2,i] <- b
      }
      
      # Perform counting in current bin.
      for (k in 3:length(df.spec.dens)) {
        if (df.spec.dens[j,k] >= a && df.spec.dens[j,k] < b) {
          bin.frame[bin.frame[1] == band,i+1] <- bin.frame[bin.frame[1] == band,i+1] + 1
        }
        else if (i == 100) {
          if (df.spec.dens[j,k] >= a && df.spec.dens[j,k] <= b) {
            bin.frame[bin.frame[1] == band,i+1] <- bin.frame[bin.frame[1] == band,i+1] + 1
          }
        }
      }
    }
  }
  
  # bin.frame returned contains 4 rows (one for each freq band), and 101 columns
  # (1 for freq band label, and 100 for bins).
	bin.limit <- cbind(bin.frame[,1], bin.limit)
  list(bin.frame, bin.limit)
}

GetProbabilites <- function(t.avg, f.log, g.log){
  # Create a data.frame to store the probabilites
  t.log <- cbind(t.avg, matrix(NA, nrow = nrow(t.avg), ncol = 2))
	
  # Iterate through the bins to find the probabilities for each subject
	# Iterate through each subject
	# i=number of subjects, j=frequency band, k=f.log bin
  for (i in 1:nrow(t.log)){
		# Iterate through delta through beta	
    for (j in 1:nrow(f.log[[1]])){  
			# Iterate through all bins of f.log
      for (k in 2:(length(f.log[[1]]))){
        if (is.na(t.log[i,4]) && (t.log[i,3] < f.log[[2]][j,k])){
          t.log[i,4] <- f.log[[1]][j,k]
        }
        if (is.na(t.log[i,5]) && (t.log[i,3] < g.log[[2]][j,k])){
          t.log[i,5] <- g.log[[1]][j,k]
        }
      }
    }
  }
	colnames(t.log)[4:5] <- c("Happy", "Neutral")
  t.log
}

# Show the happy and neutral test results.
TestSummary <- function(df){

	# Extract the subject IDs
	df.subs <- levels(unique(df[,2])[drop = TRUE])

	# Summary data.frame will contain the sum of log probabilities
	# and report the final result ("Happy" or "Neutral").
	summary <- data.frame(matrix(NA, nrow = length(df.subs), ncol = 4))
	colnames(summary) <- c("id", "happy", "neutral", "Result")

	# Cycle through the subjects, adding all log probabilities.
	for(i in 1:length(df.subs)){
		print(df.subs[i])
		summary[i,1] <- df.subs[i]
		summary[i,2] <- sum(df[df[,2]==df.subs[i],4])
		summary[i,3] <- sum(df[df[,2]==df.subs[i],5])	
	
		# Decide if the result is "Happy" or "Neutral".
		if (summary[i,2] < summary[i,3]){
			summary[i,4] <- c("Neutral")
		}
		else{
			summary[i,4] <- c("Happy")
		}
	}
	summary
}

#
# END OF FUNCTIONS
#


# Import the raw data.
combined.data <- read.csv("Flower_Gray_Combined_Data.csv", header=FALSE)
#library("plotly")

# Rename the columns in combined.data.
colnames(combined.data) <- c("type", "id", "tp", paste0(seq(1:25)))

# Extract all unique subject IDs.
subject.ids <- levels(combined.data[,2])

# Extract flower subject IDs
f.subject.ids <- combined.data[,1]=="flower"
f.subject.ids <- combined.data[f.subject.ids,2]
f.subject.ids <- factor(f.subject.ids)
f.subject.ids <- levels(f.subject.ids)

#Extract gray subject IDs
g.subject.ids <- combined.data[,1]=="gray"
g.subject.ids <- combined.data[g.subject.ids,2]
g.subject.ids <- factor(g.subject.ids)
g.subject.ids <- levels(g.subject.ids)

# Isolate the flower data
flower.rows <- combined.data[,1] == "flower"
flower.data <- combined.data[flower.rows, ]

# Isolate the gray data
gray.rows <- combined.data[,1] == "gray"
gray.data <- combined.data[gray.rows, ]


# Calculate the magnitudes, then create frames for the delta, theta, alpha, and
# beta frequencies.
fft.data <- ApplyFFT(combined.data, subject.ids)
magnitude.frame <- ApplyMagnitude(fft.data)

# Rename electrode columns back to 1..25 instead of X1..X25.
colnames(magnitude.frame) <- c("id", "tp", seq(1,25))

# Build one data frame for each significant frequency spectrum and append each
# subject's mean values for their respective electrode readings.
delta.frame <- data.frame(matrix(nrow = 0, ncol = 27))
theta.frame <- data.frame(matrix(nrow = 0, ncol = 27))
alpha.frame <- data.frame(matrix(nrow = 0, ncol = 27))
beta.frame <- data.frame(matrix(nrow = 0, ncol = 27))

# Rename initial frame columns.
colnames(delta.frame) <- c("id", "tp", seq(1,25))
colnames(theta.frame) <- c("id", "tp", seq(1,25))
colnames(alpha.frame) <- c("id", "tp", seq(1,25))
colnames(beta.frame) <- c("id", "tp", seq(1,25))

# Obtain average electrode magnitudes for every subject and every electrode.
j <- 1
for (i in seq(1,length(subject.ids))) {
  # Bind the rows to the ends of each frame.
  delta.frame <- rbind(delta.frame, magnitude.frame[(j):(j+6), ])
  theta.frame <- rbind(theta.frame, magnitude.frame[(j+7):(j+8), ])
  alpha.frame <- rbind(alpha.frame, magnitude.frame[(j+9):(j+12), ])
  beta.frame <- rbind(beta.frame, magnitude.frame[(j+13):(j+30), ])
  
  # Create a matrix for each list of electrode means.
  delta.means <- cbind(id = subject.ids[i], 
                       electrode = seq(1,25), 
                       mean = colMeans(delta.frame[3:27]))
  theta.means <- cbind(id = subject.ids[i], 
                       electrode = seq(1,25), 
                       mean = colMeans(theta.frame[3:27]))
  alpha.means <- cbind(id = subject.ids[i], 
                       electrode = seq(1,25), 
                       mean = colMeans(alpha.frame[3:27]))
  beta.means <- cbind(id = subject.ids[i], 
                      electrode = seq(1,25), 
                      mean = colMeans(beta.frame[3:27]))
  
  # Convert each means to a data frame.
  delta.means <- data.frame(delta.means)
  theta.means <- data.frame(theta.means)
  alpha.means <- data.frame(alpha.means)
  beta.means <- data.frame(beta.means)
  
  # Create a list of all the frames
  wave.means <- list(delta.means, theta.means, alpha.means, beta.means)
  
  # Output a plot for each subject.
  
  # Increment j.
  if (i < 34) {
    j <- (j+256)
  }
}

delta.plot <- ggplot(data=delta.means,
                     aes(x=electrode,
                         y=factor(mean),
                         fill=factor(mean))) +
  geom_bar(stat="identity") + 
  guides(fill=FALSE) + 
  ggtitle(subject.ids[1]) + 
  scale_x_discrete(limits=seq(1,25))

# Generate all Mean Frames for delta, theta, etc.
f.delta.frame <- ExtractSubjectData(delta.frame, f.subject.ids)
f.delta.mean <- ExtractMeans(f.delta.frame)
g.delta.frame <- ExtractSubjectData(delta.frame, g.subject.ids)
g.delta.mean <- ExtractMeans(g.delta.frame)
u.delta.mean <- UniteFrame(f.delta.mean, g.delta.mean)

f.theta.frame <- ExtractSubjectData(theta.frame, f.subject.ids)
f.theta.mean <- ExtractMeans(f.theta.frame)
g.theta.frame <- ExtractSubjectData(theta.frame, g.subject.ids)
g.theta.mean <- ExtractMeans(g.theta.frame)
u.theta.mean <- UniteFrame(f.theta.mean, g.theta.mean)

f.alpha.frame <- ExtractSubjectData(alpha.frame, f.subject.ids)
f.alpha.mean <- ExtractMeans(f.alpha.frame)
g.alpha.frame <- ExtractSubjectData(alpha.frame, g.subject.ids)
g.alpha.mean <- ExtractMeans(g.alpha.frame)
u.alpha.mean <- UniteFrame(f.alpha.mean, g.theta.mean)

f.beta.frame <- ExtractSubjectData(beta.frame, f.subject.ids)
f.beta.mean <- ExtractMeans(f.beta.frame)
g.beta.frame <- ExtractSubjectData(beta.frame, g.subject.ids)
g.beta.mean <- ExtractMeans(g.beta.frame)
u.beta.mean <- UniteFrame(f.beta.mean, g.beta.mean)

# Obtain all the spectral densities for subjects
f.spectral.densities <- GetAllDensities(fft.data, f.subject.ids[3:length(f.subject.ids)])
g.spectral.densities <- GetAllDensities(fft.data, g.subject.ids[3:length(g.subject.ids)])

# Obtain min/max values.
f.min.max <- FindMinMax(f.spectral.densities)
g.min.max <- FindMinMax(g.spectral.densities)

# Obtain bin frequencies
f.freqs <- CreateBins(f.min.max, f.spectral.densities, 100)
g.freqs <- CreateBins(g.min.max, g.spectral.densities, 100)

# Calculate bin probabilities
f.probs <- f.freqs
f.probs[[1]][,2:101] <- f.freqs[[1]][,2:101] / sum(f.freqs[[1]][1,2:101])
print(dim(f.probs[[1]]))
f.logs <- f.probs
f.logs[[1]][,2:101] <- log(f.logs[[1]][,2:101])
g.probs <- g.freqs
g.probs[[1]][,2:101] <- g.freqs[[1]][,2:101] / sum(g.freqs[[1]][1,2:101])
g.logs <- g.probs
g.logs[[1]][,2:101] <- log(g.probs[[1]][,2:101])

# Prep and test the Test Subjects
# Obtain the spectral densities for the remaining subjects
f.test.densities <- GetAllDensities(fft.data, f.subject.ids[1:2])
f.average.densities <- cbind(f.test.densities[,1:2],
                             rowMeans(f.test.densities[,3:27]))
colnames(f.average.densities)[3] <- c("avg.dens")

g.test.densities <- GetAllDensities(fft.data, g.subject.ids[1:2])
g.average.densities <- cbind(g.test.densities[,1:2], 
                             rowMeans(g.test.densities[,3:27]))
colnames(g.average.densities)[3] <- c("avg.dens")

# average.densities contains all average spectral densities, for all test subjects.
test.average.densities <- rbind(f.average.densities, g.average.densities)

# Find the log probabilities for a subject being happy and neutral for each band.
test.log <- GetProbabilites(test.average.densities, f.logs, g.logs)
test.summary <- TestSummary(test.log)
print(test.summary)

# Create plots for all mean values
delta.plot <- CreatePlots(u.delta.mean, c("Mean Delta Magnitudes"))
theta.plot <- CreatePlots(u.theta.mean, c("Mean Theta Magnitudes"))
alpha.plot <- CreatePlots(u.alpha.mean, c("Mean Alpha Magnitudes"))
beta.plot <- CreatePlots(u.beta.mean, c("Mean Beta Magnitudes"))

super.plot <- alpha.plot + 
  geom_line(data=u.beta.mean,
            aes(x=el, y=unlist(mag.mean), 
                group=Condition)) +
  geom_line(data=u.delta.mean, 
            aes(x=el, y=unlist(mag.mean), 
                group=Condition)) +
  geom_line(data=u.theta.mean, 
            aes(x=el, y=unlist(mag.mean), 
                group=Condition)) +
  ggtitle("Mean Magnitudes, All Bands")

# Single subject plot (happy), all 256 timpepoints and 25 electrodes
f.oneSub <- combined.data[combined.data[,2]==f.subject.ids[1],3:28]
f.mdat <- melt(f.oneSub)
f.mdat$tp <- seq(1:256)
f.oneSub.plot <- ggplot(data=f.mdat, 
                        aes(x=factor(tp), 
                            y=value, 
                            group=factor(variable), 
                            colour=factor(variable))) + 
  geom_line() + 
  ggtitle("One Subject, Happy, All Electrodes") + 
  xlab("Timepoint (1:256)") + 
  scale_y_continuous(limits=c(-15, 30))

# Single subject plot (neutral), all 256 timpepoints and 25 electrodes
g.oneSub <- combined.data[combined.data[,2]==g.subject.ids[1],3:28]
g.mdat <- melt(g.oneSub)
g.mdat$tp <- seq(1:256)
g.oneSub.plot <- ggplot(data=g.mdat, 
                        aes(x=factor(tp), y=value, group=factor(variable), 
                            colour=factor(variable))) + 
  geom_line() + 
  ggtitle("One Subject, Neutral, All Electrodes") + 
  xlab("Timepoint (1:256)") + 
  scale_y_continuous(limits=c(-15, 30))

# One subject (happy) transformed data, all 256 timepoints, and 25 electrodes
f.oneSub.fft <- magnitude.frame[magnitude.frame[,1]==f.subject.ids[1],2:27]
f.oneSub.fft <- f.oneSub.fft[1:35,]
f.oneSub.fft[,1] <- factor(seq(1:35))
f.mdat <- melt(f.oneSub.fft)
f.oneSub.fft.plot <- ggplot(data=f.mdat, 
                            aes(x=factor(tp), 
                                y=value, 
                                group=factor(variable), 
                                colour=factor(variable))) + 
  geom_line() + 
  ggtitle("Magnitudes, Happy, All Electrodes") + 
  xlab("Timepoint (1:35)") +
  scale_y_continuous(limit=c(0,1000))

# One subject (happy) transformed data, all 256 timepoints, and 25 electrodes
g.oneSub.fft <- magnitude.frame[magnitude.frame[,1]==g.subject.ids[1],2:27]
g.oneSub.fft <- g.oneSub.fft[1:35,]
#g.oneSub.fft[,2:26] <- rowMeans(g.oneSub.fft[,2:26])
g.oneSub.fft[,1] <- factor(seq(1:35))
g.mdat <- melt(g.oneSub.fft)
g.oneSub.fft.plot <- ggplot(data=g.mdat, 
                            aes(x=factor(tp), 
                                y=value, 
                                group=factor(variable), 
                                colour=factor(variable))) + 
  geom_line() + 
  ggtitle("Magnitudes, Neutral, All Electrodes") + 
  xlab("Timepoint (1:35)") +
  scale_y_continuous(limit=c(0,1000))

# All plots in one convenient location! Output after each.
g.oneSub.plot 
f.oneSub.plot 
g.oneSub.fft.plot
f.oneSub.fft.plot
super.plot # Needs more colour to differentiate freq. bands.
alpha.plot
beta.plot
delta.plot
theta.plot
