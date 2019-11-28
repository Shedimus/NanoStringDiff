
NanoStringDataNormalization <- function(path=path, header=TRUE, designs){
  # ------------------------------------ #
  # path: points out the directory in which your csv.file is located #
  # header: a logical value(TRUE or FALSE) indicating whether the file contains the names of the variables as its first line. #
  # designs: data frame for pheno type data storage #
  # ------------------------------------ #
  NanoStringData <- createNanoStringSetFromCsv(path, header, designs) # Create a NanoStringSet needed in NanoStringDiff from csv file #
  NanoStringData <- estNormalizationFactors(NanoStringData) # get normalization factors #
  
  PosFactor <- positiveFactor(NanoStringData) # get positive size factor #
  HouseFactor <- housekeepingFactor(NanoStringData) # get housekeeping size factor #
  NegFactor <- negativeFactor(NanoStringData) # get background niose #
  RawData <- exprs(NanoStringData) # get raw data #
  
  NormalizedData <- round((RawData - NegFactor)/(PosFactor * HouseFactor)) # normalization #
  NormalizedData[NormalizedData < 0] <- 0 # set negative value as 0 #
  
  ColName <- colnames(NormalizedData)
  colnames(NormalizedData) <- paste("Normalized", ColName, sep = "") # change column for normalized data #
  
  return(list(RawData = RawData, NormalizedData = NormalizedData))
} # end function #
