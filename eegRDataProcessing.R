# Import the raw data. Class: data.frame
Flower_Gray_Combined_Data <- read.csv("~/Documents/School/CSC_511_AI/Research_Proj/repository of eeg/ERPPaired/Flower_Gray_Combined_Data.csv", header=FALSE)
colnames(Flower_Gray_Combined_Data) <- c("Type", "Subject", "Timepoint", paste0("El", seq(1:25))) # Set column names.

# Order the data by subject, where each subject contains a data.frame of timepoint x electrode matrices.
combinedDataBySubject <- data.frame(x=1:34)
combinedDataBySubject <- combinedDataBySubject[,FALSE]

# Extract all unique subject IDs.  Class: character array.
subjectIDs <- levels(Flower_Gray_Combined_Data[,2])

# binCreator is a function that builds a matrix of subject IDs, max value, and min value.
binCreator <- function(df, subjMat){
  binTotal <- 32
  
  # Create a dataframe to store the bin values in order from min to max.
  initialFrame <- data.frame(x=1:33)
  binFrame <- initialFrame[,FALSE]
  print(binFrame)
  # Loop through each subject in the data.frame.
  for(i in 1:length(subjMat)){
    # Loop through each set of electrode potentials.
    for(j in 4:length(df)){
      sortedPotentials <- sort(df[subjMat[i]==df[,2],j], decreasing = FALSE) # Sort array of electrode potentials for subject i.
      maxPotential <- sortedPotentials[256] # Max electrode EEG value.
      minPotential <- sortedPotentials[1] # Min electrode EEG value.

      # Calculate the deltaV (max - min)/32 for each subject set.
      deltaV <- (maxPotential-minPotential)/binTotal
      binArray <- matrix(0, nrow = 33, ncol = 1) # Matrix ore all bin values.
      colnames(binArray) <- paste0(subjMat[i],"-El",(j-3)) # Name the matrix columns.
      binValue <-  minPotential # Set the initial bin variable value.
      
      # Build an array of the BIN values
      for(k in 1:(binTotal+1)){
        #print(binValue)
        binArray[k] <- binValue # Add a new bin value to the binArray.
        binValue <- binValue + deltaV # Increment the next bin.
      }
      # Add the sorted data as a new column to binFrame.
      
      binFrame <- cbind(binFrame, binArray) # Add the binArray as a new column to the Frame.
    }
  }
  row.names(binFrame) = paste0("BIN",seq(1:33))
  binFrame # Return the binFrame
}


binFrequency <- function(df, binF, subjMat){
  initialFrame <- data.frame(x=1:33)
  frequencyFrame <- initialFrame[,FALSE] # Create the frequency frame
  # Loop through all subjects.
  for(i in 2:length(subjMat)){
    binArray <- matrix(0, nrow = 33, ncol = 1) # Matrix ore all bin values.
    
    # Loop through all electrodes.
    for(j in 1:25){
       # Loop through all bins, add the count of values < the bin value to the frequencyFrame.
       for(k in 1:33){
         binArray[k] <- count()
       }
    }
  }
}


# Isolate the flower data
flowerRows <- Flower_Gray_Combined_Data[,1] == "flower"
flowerOnlyData <- Flower_Gray_Combined_Data[flowerRows,]
# Isolate the gray data
grayRows <- Flower_Gray_Combined_Data[,1] == "gray"
grayOnlyData <- Flower_Gray_Combined_Data[grayRows,]

# Performs Full Fourier Transform on each set of electrodes for each subject.
# The function accepts a matrix of subject IDs (subjIDs), and data.frame (df).
EegFFT <- function(df, subjIDs){
    initialFrame <- data.frame(x=1:256)
    
    # frequencyFrame stores all fft data in a 256 x (34*25) data frame.
    fftFrame <- initialFrame[,FALSE] # Create the frequency frame
  
    for(i in 1:length(subjIDs)){
        for(j in 4:28){
            # create a matrix of fft values for 1 electrode, 1 subject.
            fftArray <- matrix(fft(df[df[,2]==subjIDs[i],j]), nrow = 256, ncol = 1)
            colnames(fftArray) <- paste0(subjIDs[i],"-El",(j-3)) # Name the matrix columns.
            # Store the frequencies in the frequencyFrame.
            fftFrame <- cbind(fftFrame, fftArray)
        }
    }
    row.names(fftFrame) = paste0(1:256)
    fftFrame # Return the binFrame
}

# To write the fft file to a .csv file
# write.csv(x, file = "~/Documents/School/CSC_511_AI/Research_Proj/repository of eeg/ERPPaired/FFT_Flower_Gray_Combined.csv")

# Calculate the magnitudes, then create frames for the delta, theta, alpha, and beta frequencies
fft_data <- EegFFT(Flower_Gray_Combined_Data, subjectIDs)
magnitudeFrame <- abs(fft_data[,2:length(fft_data)])
deltaFrame <- magnitudeFrame[1:7,]
thetaFrame <- magnitudeFrame[8:9,]
alphaFrame <- magnitudeFrame[10:13,]
betaFrame <- magnitudeFrame[14:31,]

# Obtain the mean of magnitudes for each electrode from magnitudes stored in data.frames (deltaFrame, etc.).
delta <- colMeans(deltaFrame)
theta <- colMeans(thetaFrame)
alpha <- colMeans(alphaFrame)
beta <- colMeans(betaFrame)
