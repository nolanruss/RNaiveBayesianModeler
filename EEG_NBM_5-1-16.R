# Required packages.
library("ggplot2")
library("stats")
library("reshape2")
library("plyr")
library("datasets")
library("graphics")
library("methods")
library("stats")
library("utils")

# Variables
bin.tot <- 100 # Number of frequency bins to use with the model
subj.tot <- 8 # Number of subjects to test against the model

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

# Calculate the Amplitude for all members of an fft dataset.
ApplyAmplitude <- function(df){
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
  y <- data.frame(Electrode=seq(1:25))
  y <- cbind(y, Amplitude=matrix(df.subset[1,1:25], nrow=25, ncol=1))
  y
}

# Reunite the delta values into a single df. First frame should be "happy", 2nd "neutral"
UniteFrame <- function(h.df, n.df, df.names){
  df.united <- data.frame(Condition=rep(df.names, times=1, each=25))
  df.united <- cbind(df.united, rbind(h.df, n.df))
  df.united
}

# Accept a data.frame with the layout created by UniteFrame() and 
# return a list of plots.
CreatePlots <- function(df, title){
  # Unlist the Amplitude means and conditions
  df$Amplitude <- unlist(df$Amplitude)
  df$Condition <- unlist(df$Condition)

  # Create the plot
  tmp.plot <- ggplot(data=df, 
                     aes(x=Electrode, 
                         y=Amplitude, 
                         group = Condition, 
                         colour = Condition)) + 
                        geom_line() + 
                        ggtitle(title)
  tmp.plot
}

# Spectral density function. Takes a subject's Amplitude frame and a frequency
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
  # ApplyAmplitude to df.
  df <- ApplyAmplitude(df)
  
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
  colnames(bin.frame) <- c("freq.band", seq(1,bin.total))
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
        bin.limit[3,i] <- b
      }    
      else if (df.spec.dens[j,1] == "beta") {
        a <- beta.min + ((i-1) * beta.width)
        b <- beta.min + (i * beta.width)
        band <- "beta"
        bin.limit[4,i] <- b
      }
      
      # Perform counting in current bin.
      for (k in 3:length(df.spec.dens)) {
        if (df.spec.dens[j,k] >= a && df.spec.dens[j,k] < b) {
          bin.frame[bin.frame[,1] == band,i+1] <- bin.frame[bin.frame[,1] == band,i+1] + 1
        }
        else if (i == bin.total) {
          if (df.spec.dens[j,k] >= a && df.spec.dens[j,k] <= b) {
            bin.frame[bin.frame[,1] == band,i+1] <- bin.frame[bin.frame[,1] == band,i+1] + 1
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
  t.log <- cbind(t.avg, matrix(-Inf, nrow = nrow(t.avg), ncol = 2))
  
  # Iterate through the bins to find the probabilities for each subject
	# Iterate through each subject
	# i=number of subjects, j=frequency band, k=f.log bin
  for (i in 1:nrow(t.log)){
		# Iterate through delta through beta	
    for (j in 1:nrow(f.log[[1]])){  
			# Iterate through all bins of f.log and g.log
      for (k in 2:(length(f.log[[1]]))){
        if (is.infinite(t.log[i,4]) && (t.log[i,3] < f.log[[2]][j,k])){
          t.log[i,4] <- f.log[[1]][j,k]
        }
        if (is.infinite(t.log[i,5]) && (t.log[i,3] < g.log[[2]][j,k])){
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
		summary[i,1] <- df.subs[i]
		summary[i,2] <- sum(df[df[,2]==df.subs[i],4])
		summary[i,3] <- sum(df[df[,2]==df.subs[i],5])	
	
		# Determine if the result is "Happy" or "Neutral".
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


# Calculate the Amplitudes, then create frames for the delta, theta, alpha, and
# beta frequencies.
fft.data <- ApplyFFT(combined.data, subject.ids)
Amplitude.frame <- ApplyAmplitude(fft.data)

# Rename electrode columns back to 1..25 instead of X1..X25.
colnames(Amplitude.frame) <- c("id", "tp", seq(1,25))

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

# Obtain average electrode Amplitudes for every subject and every electrode.
j <- 1
for (i in seq(1,length(subject.ids))) {
  # Bind the rows to the ends of each frame.
  delta.frame <- rbind(delta.frame, Amplitude.frame[(j):(j+6), ])
  theta.frame <- rbind(theta.frame, Amplitude.frame[(j+7):(j+8), ])
  alpha.frame <- rbind(alpha.frame, Amplitude.frame[(j+9):(j+12), ])
  beta.frame <- rbind(beta.frame, Amplitude.frame[(j+13):(j+30), ])
  
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

# Generate all Mean Frames for delta, theta, etc.
f.delta.frame <- ExtractSubjectData(delta.frame, f.subject.ids)
f.delta.mean <- ExtractMeans(f.delta.frame)
g.delta.frame <- ExtractSubjectData(delta.frame, g.subject.ids)
g.delta.mean <- ExtractMeans(g.delta.frame)
u.delta.mean <- UniteFrame(f.delta.mean, g.delta.mean, c("Delta (Happy)", "Delta (Neutral)"))

f.theta.frame <- ExtractSubjectData(theta.frame, f.subject.ids)
f.theta.mean <- ExtractMeans(f.theta.frame)
g.theta.frame <- ExtractSubjectData(theta.frame, g.subject.ids)
g.theta.mean <- ExtractMeans(g.theta.frame)
u.theta.mean <- UniteFrame(f.theta.mean, g.theta.mean, c("Theta (Happy)", "Theta (Neutral)"))

f.alpha.frame <- ExtractSubjectData(alpha.frame, f.subject.ids)
f.alpha.mean <- ExtractMeans(f.alpha.frame)
g.alpha.frame <- ExtractSubjectData(alpha.frame, g.subject.ids)
g.alpha.mean <- ExtractMeans(g.alpha.frame)
u.alpha.mean <- UniteFrame(f.alpha.mean, g.alpha.mean, c("Alpha (Happy)", "Alpha (Neutral)"))

f.beta.frame <- ExtractSubjectData(beta.frame, f.subject.ids)
f.beta.mean <- ExtractMeans(f.beta.frame)
g.beta.frame <- ExtractSubjectData(beta.frame, g.subject.ids)
g.beta.mean <- ExtractMeans(g.beta.frame)
u.beta.mean <- UniteFrame(f.beta.mean, g.beta.mean, c("Beta (Happy)", "Beta (Neutral)"))

# Obtain all the spectral densities for subjects
f.spectral.densities <- GetAllDensities(fft.data, f.subject.ids[1:length(f.subject.ids)])
g.spectral.densities <- GetAllDensities(fft.data, g.subject.ids[1:length(g.subject.ids)])

# Obtain min/max values.
u.spectral.densities <- rbind(f.spectral.densities, g.spectral.densities)
u.min.max <- FindMinMax(u.spectral.densities)
rownames(u.min.max) <- c("Min", "Max")

f.min.max <- FindMinMax(f.spectral.densities)
g.min.max <- FindMinMax(g.spectral.densities)
rownames(f.min.max) <- c("Happy-min", "Happy-max")
rownames(g.min.max) <- c("Neutral-min", "Neutral-max")
mmax <- rbind(f.min.max, g.min.max)

# Obtain bin frequencies
u.f.freqs <- CreateBins(u.min.max, f.spectral.densities, bin.tot)
u.g.freqs <- CreateBins(u.min.max, g.spectral.densities, bin.tot)
f.freqs <- CreateBins(f.min.max, f.spectral.densities, bin.tot)
g.freqs <- CreateBins(g.min.max, g.spectral.densities, bin.tot)

# Calculate bin probabilities
bin.lim <- bin.tot + 1
f.probs <- f.freqs
f.probs[[1]][,2:bin.lim] <- f.freqs[[1]][,2:bin.lim] / sum(f.freqs[[1]][1,2:bin.lim])
print(dim(f.probs[[1]]))
f.logs <- f.probs
f.logs[[1]][,2:bin.lim] <- log(f.logs[[1]][,2:bin.lim])
g.probs <- g.freqs
g.probs[[1]][,2:bin.lim] <- g.freqs[[1]][,2:bin.lim] / sum(g.freqs[[1]][1,2:bin.lim])
g.logs <- g.probs
g.logs[[1]][,2:bin.lim] <- log(g.probs[[1]][,2:bin.lim])

# U. Tryout
u.f.probs <- u.f.freqs
u.f.probs[[1]][,2:bin.lim] <- u.f.freqs[[1]][,2:bin.lim] / sum(u.f.freqs[[1]][1,2:bin.lim])
print(dim(u.f.probs[[1]]))
u.f.logs <- u.f.probs
u.f.logs[[1]][,2:bin.lim] <- log(u.f.logs[[1]][,2:bin.lim])
u.g.probs <- u.g.freqs
u.g.probs[[1]][,2:bin.lim] <- u.g.freqs[[1]][,2:bin.lim] / sum(u.g.freqs[[1]][1,2:bin.lim])
u.g.logs <- u.g.probs
u.g.logs[[1]][,2:bin.lim] <- log(u.g.probs[[1]][,2:bin.lim])


# Test all electrodes against the model
ElectrodeSummary <- data.frame(matrix(NA, nrow = 0, ncol = 4))
electrodes <-25
for(k in 1:subj.tot){
  j.range <- 9 + 2
#  for(j in 3:j.range){
    f.oneEl <- fft.data[,1:2]
    for(i in 1:electrodes){
      f.oneEl <- cbind(f.oneEl, fft.data$`9`)
    }
    g.oneEl <- fft.data[,1:2]
    for(i in 1:electrodes){
      g.oneEl <- cbind(g.oneEl, fft.data$`9`)
    }
    
    f.oneEl.densities <- GetAllDensities(f.oneEl, f.subject.ids[k])
    g.oneEl.densities <- GetAllDensities(g.oneEl, g.subject.ids[k])
#    print(head(f.oneEl.densities))
    oneEl.densities <- rbind(f.oneEl.densities, g.oneEl.densities)
#    print(oneEl.densities)
    oneEl.avg <- oneEl.densities[,1:3]
    oneEl.oneB <- oneEl.avg[oneEl.avg[,1]=="delta",]
    for(n in 1:3){
      oneEl.oneB <- rbind(oneEl.oneB, oneEl.oneB[k,])
      print(oneEl.oneB)
    }
    print(oneEl.oneB)
 #   print(oneEl.avg)
    oneEl.log <- GetProbabilites(oneEl.avg, u.f.logs, u.g.logs)
    oneEl.summary <- TestSummary(oneEl.log)
    rownames(oneEl.summary) <- list(c("H-",j.range-2), c("N-",(j.range-2)))
    ElectrodeSummary <- rbind(ElectrodeSummary, oneEl.summary)
#  }
}
View(ElectrodeSummary)


# Prep and test the Test Subjects
# Obtain the spectral densities for the remaining subjects
f.test.densities <- GetAllDensities(fft.data, f.subject.ids[1:subj.tot])
f.average.densities <- cbind(f.test.densities[,1:2],
                             rowMeans(f.test.densities[,3:27]))
colnames(f.average.densities)[3] <- c("avg.dens")

g.test.densities <- GetAllDensities(fft.data, g.subject.ids[1:subj.tot])
g.average.densities <- cbind(g.test.densities[,1:2], 
                             rowMeans(g.test.densities[,3:27]))
colnames(g.average.densities)[3] <- c("avg.dens")

# average.densities contains all average spectral densities, for all test subjects.
test.average.densities <- rbind(f.average.densities, g.average.densities)

# Find the log probabilities for a subject being happy and neutral for each band.
test.log <- GetProbabilites(test.average.densities, f.logs, g.logs)
test.summary <- TestSummary(test.log)
View(test.summary)

# Test using uniform bin sizes and widths for happy and neutral (prob given either happy or neutral)
u.test.log <- GetProbabilites(test.average.densities, u.f.logs, u.g.logs)
u.test.summary <- TestSummary(u.test.log)
View(u.test.summary)

# Create plots for all mean values
delta.plot <- CreatePlots(u.delta.mean, c("Mean Delta Amplitudes"))
theta.plot <- CreatePlots(u.theta.mean, c("Mean Theta Amplitudes"))
alpha.plot <- CreatePlots(u.alpha.mean, c("Mean Alpha Amplitudes"))
beta.plot <- CreatePlots(u.beta.mean, c("Mean Beta Amplitudes"))

# Select a better color pallete, and create the plot of all Amplitudes.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#D55E00", "#0072B2")
super.plot <- alpha.plot + 
              geom_line(data=u.beta.mean, 
                        aes(x=Electrode,
                            y=unlist(Amplitude), 
                            group=unlist(Condition))) +
              geom_line(data=u.delta.mean, 
                        aes(x=Electrode, 
                            y=unlist(Amplitude), 
                            group=unlist(Condition))) +
              geom_line(data=u.theta.mean, 
                        aes(x=Electrode, 
                            y=unlist(Amplitude), 
                            group=unlist(Condition))) +
              ggtitle("Mean Amplitudes, All Bands") +
              scale_color_manual(values=cbbPalette)

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
                        aes(x=factor(tp), 
                            y=value, 
                            group=factor(variable), 
                            colour=factor(variable))) + 
                        geom_line() + 
                        ggtitle("One Subject, Neutral, All Electrodes") + 
                        xlab("Timepoint (1:256)") + 
                        scale_y_continuous(limits=c(-15, 30))

# One subject (happy) transformed data, all 256 timepoints, and 25 electrodes
f.oneSub.fft <- Amplitude.frame[Amplitude.frame[,1]==f.subject.ids[1],2:27]
f.oneSub.fft <- f.oneSub.fft[1:35,]
f.oneSub.fft[,1] <- factor(seq(1:35))
f.mdat <- melt(f.oneSub.fft)
f.oneSub.fft.plot <- ggplot(data=f.mdat, 
                            aes(x=factor(tp), 
                                y=value, 
                                group=factor(variable), 
                                colour=factor(variable))) + 
                            geom_line() + 
                            ggtitle("Amplitudes, Happy, All Electrodes") + 
                            xlab("Timepoint (1:35)") +
                            scale_y_continuous(limit=c(0,1000))

# One subject (happy) transformed data, all 256 timepoints, and 25 electrodes
g.oneSub.fft <- Amplitude.frame[Amplitude.frame[,1]==g.subject.ids[1],2:27]
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
                            ggtitle("Amplitudes, Neutral, All Electrodes") + 
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


# Generate plots for the log probabilites of both subject types, for each subject type
f.d.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(f.d.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  f.d.log[i-1, 1] <- c("Delta-Happy")
  f.d.log[i-1, 2] <- u.f.logs[[2]][1,i]
  f.d.log[i-1, 3] <- u.f.logs[[1]][1,i]
}

g.d.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(g.d.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  g.d.log[i-1, 1] <- c("Delta-Neutral")
  g.d.log[i-1, 2] <- u.g.logs[[2]][1,i]
  g.d.log[i-1, 3] <- u.g.logs[[1]][1,i]
}

# Theta Frame
f.t.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(f.t.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  f.t.log[i-1, 1] <- c("Theta-Happy")
  f.t.log[i-1, 2] <- u.f.logs[[2]][2,i]
  f.t.log[i-1, 3] <- u.f.logs[[1]][2,i]
}

g.t.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(g.t.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  g.t.log[i-1, 1] <- c("Theta-Neutral")
  g.t.log[i-1, 2] <- u.g.logs[[2]][2,i]
  g.t.log[i-1, 3] <- u.g.logs[[1]][2,i]
}

# Alpha Frame
f.a.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(f.a.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  f.a.log[i-1, 1] <- c("Alpha-Happy")
  f.a.log[i-1, 2] <- u.f.logs[[2]][3,i]
  f.a.log[i-1, 3] <- u.f.logs[[1]][3,i]
}

g.a.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(g.a.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  g.a.log[i-1, 1] <- c("Alpha-Neutral")
  g.a.log[i-1, 2] <- u.g.logs[[2]][3,i]
  g.a.log[i-1, 3] <- u.g.logs[[1]][3,i]
}

# Beta frame
f.b.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(f.b.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  f.b.log[i-1, 1] <- c("Beta-Happy")
  f.b.log[i-1, 2] <- u.f.logs[[2]][4,i]
  f.b.log[i-1, 3] <- u.f.logs[[1]][4,i]
}

g.b.log <- data.frame(NA, nrow = (length(u.f.logs[[1]])-1), ncol = 3)
colnames(g.b.log) <- c("type", "bin.val", "bin.freq")

for(i in 2:length(u.f.logs[[1]])){
  g.b.log[i-1, 1] <- c("Beta-Neutral")
  g.b.log[i-1, 2] <- u.g.logs[[2]][4,i]
  g.b.log[i-1, 3] <- u.g.logs[[1]][4,i]
}

d.groupFrame <- rbind(f.d.log, g.d.log)
t.groupFrame <- rbind(f.t.log, g.t.log)
a.groupFrame <- rbind(f.a.log, g.a.log)
b.groupFrame <- rbind(f.b.log, g.b.log)

d.log.plot <- ggplot(data=d.groupFrame, 
                    aes(x=factor(bin.val), 
                        y=bin.freq, 
                        group=factor(type), 
                        colour=factor(type))) + 
  geom_line() + 
  ggtitle("Delta Log Probabilites vs. Bin") + 
  xlab("Spectral Density Bin Limits") +
  ylab("Log Probability")

t.log.plot <- ggplot(data=t.groupFrame, 
                     aes(x=factor(bin.val), 
                         y=bin.freq, 
                         group=factor(type), 
                         colour=factor(type))) + 
  geom_line() + 
  ggtitle("Theta Log Probabilites vs. Bin") + 
  xlab("Spectral Density Bin Limits") +
  ylab("Log Probability")

a.log.plot <- ggplot(data=a.groupFrame, 
                     aes(x=factor(bin.val), 
                         y=bin.freq, 
                         group=factor(type), 
                         colour=factor(type))) + 
  geom_line() + 
  ggtitle("Alpha Log Probabilites vs. Bin") + 
  xlab("Spectral Density Bin Limits") +
  ylab("Log Probability")

b.width <- (b.groupFrame[2,2]-b.groupFrame[1,2])
#b.groupFrame[,2] <- b.groupFrame[,2]-(b.width/2)
b.log.plot <- ggplot(data=b.groupFrame, 
                     aes(x=factor(bin.val), 
                         y=bin.freq, 
                         group=factor(type), 
                         colour=factor(type))) + 
  ggtitle("Beta Log Probabilites vs. Bin") + 
  xlab("Spectral Density Bin Limits") +
  ylab("Log Probability")

g.log.plot <- ggplot(data=b.groupFrame, 
                     aes(x=bin.freq, 
                         group=factor(type), 
                         colour=factor(type))) + 
  geom_ + 
  ggtitle("Beta Log Probabilites vs. Bin") + 
  xlab("Spectral Density Bin Limits") +
  ylab("Log Probability")

#head(b.groupFrame)
#x <- b.groupFrame[2,2]-b.groupFrame[1,2]
