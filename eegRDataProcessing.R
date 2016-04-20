# Import the raw data. Class: data.frame
Flower_Gray_Combined_Data <- read.csv("~/Documents/School/CSC_511_AI/Research_Proj/repository of eeg/ERPPaired/Flower_Gray_Combined_Data.csv", header=FALSE)
colnames(Flower_Gray_Combined_Data) <- c("Type", "Subject", "Timepoint", paste0("El", seq(1:25))) # Set column names.

# Order the data by subject, where each subject contains a data.frame of timepoint x electrode matrices.
combinedDataBySubject <- data.frame(x=1:34)
combinedDataBySubject <- combinedDataBySubject[,FALSE]

# Extract all unique subject IDs.  Class: character array.
subjectIDs <- levels(Flower_Gray_Combined_Data[,2])
fSubjectIDs <- levels(Flower_Gray_Combined_Data["flower",2])
gSubjectIDs <- levels(Flower_Gray_Combined_Data["gray",2])

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
            tmpDf <- df[df[,2]==subjIDs[i],j]
            if(length(tmpDf != 0)){
                fftArray <- matrix(fft(df[df[,2]==subjIDs[i],j]), nrow = 256, ncol = 1)
                colnames(fftArray) <- paste0(subjIDs[i],"-El",(j-3)) # Name the matrix columns.
                # Store the frequencies in the frequencyFrame.
                fftFrame <- cbind(fftFrame, fftArray)
            }
        }
    }
    row.names(fftFrame) = paste0(1:256)
    fftFrame # Return the binFrame
}

# To write the fft file to a .csv file
# write.csv(x, file = "~/Documents/School/CSC_511_AI/Research_Proj/repository of eeg/ERPPaired/FFT_Flower_Gray_Combined.csv")

# For flower, calculate the alpha, beta, delta, theta values.
# Calculate the magnitudes, then create frames for the delta, theta, alpha, and beta frequencies
fFFTData <- EegFFT(flowerOnlyData, fSubjectIDs)
fMagnitudeFrame <- abs(fFFTData[,1:length(fFFTData)])
fDeltaFrame <- fMagnitudeFrame[1:7,]
fThetaFrame <- fMagnitudeFrame[8:9,]
fAlphaFrame <- fMagnitudeFrame[10:13,]
fBetaFrame <- fMagnitudeFrame[14:31,]

# Obtain the magnitudes for each electrode from magnitudes stored in data.frames (deltaFrame, etc.).
fDelta <- colMeans(fDeltaFrame)
fTheta <- colMeans(fThetaFrame)
fAlpha <- colMeans(fAlphaFrame)
fBeta <- colMeans(fBetaFrame)

# For gray, calculate the alpha, beta, delta, theta values.
# Calculate the magnitudes, then create frames for the delta, theta, alpha, and beta frequencies
gFFTData <- EegFFT(grayOnlyData, fSubjectIDs)
gMagnitudeFrame <- abs(gFFTData[,1:length(gFFTData)])
gDeltaFrame <- gMagnitudeFrame[1:7,]
gThetaFrame <- gMagnitudeFrame[8:9,]
gAlphaFrame <- gMagnitudeFrame[10:13,]
gBetaFrame <- gMagnitudeFrame[14:31,]

# Obtain the magnitudes for each electrode from magnitudes stored in data.frames (deltaFrame, etc.).
gDelta <- colMeans(gDeltaFrame)
gTheta <- colMeans(gThetaFrame)
gAlpha <- colMeans(gAlphaFrame)
gBeta <- colMeans(gBetaFrame)


# Plot the data
xAxis <- rep(1:25, each = 1, times = 17) # set the axis values 1-25, repeated for each subject.
fDeltaMean <- data.frame(fDelta) # Place all data in a data.frame.
fDeltaPlot <- ggplot(fDeltaMean, aes(xAxis, y=fDelta))+geom_point() # Make a scatter plot.
gDeltaMean <- data.frame(gDelta)
gDeltaPlot <- ggplot(gDeltaMean, aes(xAxis, y=gDelta))+geom_point()

fThetaMean <- data.frame(fTheta)
fThetaPlot <- ggplot(fThetaMean, aes(xAxis, y=fTheta))+geom_point()

gThetaMean <- data.frame(gTheta)
gThetaPlot <- ggplot(gThetaMean, aes(xAxis, y=gTheta))+geom_point()

gAlphaMean <- data.frame(gAlpha)
gAlphaPlot <- ggplot(gAlphaMean, aes(xAxis, y=gAlpha))+geom_point()

fAlphaMean <- data.frame(fAlpha)
fAlphaPlot <- ggplot(fAlphaMean, aes(xAxis, y=fAlpha))+geom_point()

gBetaMean <- data.frame(gBeta)
gBetaPlot <- ggplot(gBetaMean, aes(xAxis, y=gBeta))+geom_point()

fBetaMean <- data.frame(fBeta)
fBetaPlot <- ggplot(fBetaMean, aes(xAxis, y=fBeta))+geom_point()

gDeltaPlot
fDeltaPlot
gThetaPlot
fThetaPlot
gAlphaPlot
fAlphaPlot
gBetaPlot
fBetaPlot

dSubjectMeans <- data.frame(matrix(NA, nrow=25, ncol=17))
for(i in 1:length(fDeltaMean)){
    for(j in 1:17){
        dSubjectMeans[]
    }
}

