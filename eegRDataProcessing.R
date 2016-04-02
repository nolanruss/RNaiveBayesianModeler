# Import the raw data. Class: data.frame
Flower_Gray_Combined_Data <- read.csv("~/Documents/School/CSC_511_AI/Research_Proj/repository of eeg/ERPPaired/Flower_Gray_Combined_Data.csv", header=FALSE)
colnames(Flower_Gray_Combined_Data) <- c("Type", "Subject", "Timepoint", paste0("El", seq(1:25))) # Set column names.

# Extract all unique subject IDs.  Class: character array.
subjectIDs <- levels(Flower_Gray_Combined_Data[,2])

# binCreator is a function that builds a matrix of subject IDs, max value, and min value.
binCreator <- function(df, subjMat){
  binTotal <- 32
  
  # Create a dataframe to store the bin values in order from min to max.
  binFrame <- data.frame(paste0("BIN",seq(1:33)), row.names = )
  
  # Loop through each subject in the data.frame.
  for(i in 1:length(subjMat)){
    # Loop through each set of electrode potentials.
    for(j in 4){
    #for(j in 4:length(df)){
      sortedPotentials <- sort(df[subjMat[i]==df[,2],j], decreasing = FALSE) # Sort array of electrode potentials for subject i.
      maxPotential <- sortedPotentials[256] # Max electrode EEG value.
      minPotential <- sortedPotentials[1] # Min electrode EEG value.
 #     print(sortedPotentials)
      # Calculate the deltaV (max - min)/32 for each subject set.
      deltaV <- (maxPotential-minPotential)/binTotal
      binArray <- matrix(0, nrow = 0, ncol = 1) # Stores all bin values.
      colnames(binArray) <- paste0(subjMat[i],"-El",(j-3))
      binValue <-  minPotential # Set the initial bin variable value.
      
      # Build an array of the BIN values
      #print(sortedPotentials)
      for(k in 1:(binTotal+1)){
        #print(binValue)
        binArray <- rbind(binArray, binValue) # Add a new bin value to the binArray.
        binValue <- binValue + deltaV # Increment the next bin.
      }
      print(binArray)
      binFrame <- cbind(binFrame[1,1:i], binArray, binFrame[1,-(1:i)]) # Add the binArray as a new column to the Frame.
      #print(binFrame)
    }
  }
  binFrame # Return the binFrame
}


# Isolate the flower data
flowerRows <- Flower_Gray_Combined_Data[,1] == "flower"
flowerOnlyData <- Flower_Gray_Combined_Data[flowerRows,]
# Isolate the gray data
grayRows <- Flower_Gray_Combined_Data[,1] == "gray"
grayOnlyData <- Flower_Gray_Combined_Data[grayRows,]

