## A script for creating the dataset for machine learning --------------
# Modified dataset created as a part of the Capstone Project for Harvard X professional data science certificate
# Author: Tero Jalkanen
#####

## Packages -------------------------

# set options for downloading old Bioconductor packages
#options(repos = c(CRAN = "https://mran.microsoft.com/snapshot/2020-10-28"))
#BiocManager::install("ChemmineR", force = TRUE)
#BiocManager::install("Rcpi", force = TRUE)
#BiocManager::install("ChemmineOB", force = TRUE)


#Packages for calculating small molecule properties
require(ChemmineR)
require(Rcpi)

## Load original dataset ------------
# Original data from https://doi.org/10.6084/m9.figshare.8038913
# By Xavier Domingo-Almenara et al.

# Data created by Domingo-Almenara et al.
df <- read.csv(file = "data/SMRT_dataset.csv", header = TRUE, sep = ";", dec = ".")

# Let's create molecular descriptors for the molecules in the original dataset

## Take chunks of the dataset to create the descriptors --------

chunk_size <- 100

# set-up indices
low_index <- 1
high_index <- chunk_size

# new data.frame with molecular descriptors
new_df <- NULL

#loop to create decriptors chunk-by-chunk
while(high_index<=length(df$pubchem)){
  #take chunk from the entire data.frame
  small_df <- df[low_index:high_index,]
  
  # Get test compounds as SDF
  # This SDFset can get quite large, hence the breaking down to smaller chunks
  compounds <- pubchemCidToSDF(as.numeric(small_df$pubchem))
  
  # Calculate property matrix using ChemmineR
  propma <- data.frame(MF=MF(compounds, addH=FALSE), #Molecular formula 
                       MW=MW(compounds, addH=FALSE), #molecular weight
                       Ncharges=sapply(bonds(compounds, type="charge"), length),
                       atomcountMA(compounds, addH=FALSE), #count atoms 
                       groups(compounds, type="countMA"), #count groups
                       rings(compounds, upper=6, type="count", arom=TRUE)) #count rings with max 6 carbons
  
  # Write molecules in the chunk as a temporary SDF file
  write.SDF(compounds, file = "data/temp.sdf")
  # Read molecules for Rcpi-package
  mols <- readMolFromSDF(sdffile = "data/temp.sdf")
  
  # Calculate descriptors with Rcpi
  descriptors <- cbind(
    Rcpi::extractDrugALOGP(mols), # Atom additive logP and molar refractivity values descriptor
    Rcpi::extractDrugApol(mols), # Sum of atomix polarizabilities
    Rcpi::extractDrugAromaticAtomsCount(mols), # Number of aromatic atoms
    Rcpi::extractDrugAromaticBondsCount(mols), # Number of aromatic bonds
    Rcpi::extractDrugTPSA(mols), # Topological polarizable surface area
    Rcpi::extractDrugFragmentComplexity(mols), # Complexity of a system
    Rcpi::extractDrugHBondAcceptorCount(mols), # Number of hydrogen bond acceptors
    Rcpi::extractDrugHBondDonorCount(mols), # Number of hydrogen bond donors
    Rcpi::extractDrugRotatableBondsCount(mols), # Number of non-rotatable bonds on a molecule
    Rcpi::extractDrugVABC(mols), # Volume of a molecule
    Rcpi::extractDrugZagrebIndex(mols), # Sum of the squared atom degrees of all heavy atoms
    Rcpi::extractDrugECI(mols), # the Eccentric Connectivity Index Descriptor
    Rcpi::extractDrugWienerNumbers(mols) # Wiener path number and wiener polarity number
  )
  
  #temporary data.frame with all descriptors for the chunk
  temp <- cbind(small_df, propma, descriptors)
  
  #add temporary data to the new data.frame
  if(length(new_df) == 0){
    #For the first chunk
    new_df <- temp
  } else {
    # full_join instead of rbind as temp data.frame can have varying amount of columns depending on the molecules
    new_df <- dplyr::full_join(new_df, temp)
  }

  
  # Increase indices by chunk size
  low_index <- low_index + chunk_size
  high_index <- high_index + chunk_size
}

## Calculate descriptors for final smaller chunk if necessary --------------------

if(length(new_df$pubchem) != length(df$pubchem)){
  # take the final chunk, which is smaller than the chunk size
  high_index <- length(df$pubchem)
  
  #take chunk from the entire data.frame
  small_df <- df[low_index:high_index,]
  
  # Get test compounds as SDF
  compounds <- pubchemCidToSDF(as.numeric(small_df$pubchem))
  
  # Calculate property matrix using ChemmineR
  propma <- data.frame(MF=MF(compounds, addH=FALSE), #Molecular formula 
                       MW=MW(compounds, addH=FALSE), #molecular weight
                       Ncharges=sapply(bonds(compounds, type="charge"), length),
                       atomcountMA(compounds, addH=FALSE), #count atoms 
                       groups(compounds, type="countMA"), #count groups
                       rings(compounds, upper=6, type="count", arom=TRUE)) #count rings
  
  # Write molecules in the chunk as a temporary SDF file
  write.SDF(compounds, file = "data/temp.sdf")
  # Read molecules for Rcpi-package
  mols <- readMolFromSDF(sdffile = "data/temp.sdf")
  
  # Calculate descriptors with Rcpi
  descriptors <- cbind(
    Rcpi::extractDrugALOGP(mols), # Atom additive logP and molar refractivity values descriptor
    Rcpi::extractDrugApol(mols), # Sum of atomix polarizabilities
    Rcpi::extractDrugAromaticAtomsCount(mols), # Number of aromatic atoms
    Rcpi::extractDrugAromaticBondsCount(mols), # Numer of aromatic bonds
    Rcpi::extractDrugTPSA(mols), # Topological polarizable surface area
    Rcpi::extractDrugFragmentComplexity(mols), # Complexity of a system
    Rcpi::extractDrugHBondAcceptorCount(mols), # Number of hydrogen bond acceptors
    Rcpi::extractDrugHBondDonorCount(mols), # Number of hydrogen bond donors
    Rcpi::extractDrugRotatableBondsCount(mols), # Number of non-rotatable bonds on a molecule
    Rcpi::extractDrugVABC(mols), # Volume of a molecule
    Rcpi::extractDrugZagrebIndex(mols), # Sum of the squared atom degrees of all heavy atoms
    Rcpi::extractDrugECI(mols), # the Eccentric Connectivity Index Descriptor
    Rcpi::extractDrugWienerNumbers(mols) # Wiener path number and wiener polarity number
  )
  
  #temporary data.frame with all descriptors
  temp <- cbind(small_df, propma, descriptors)
  
  # new data frame with all molecules and new descriptors
  new_df <- dplyr::full_join(new_df, temp)
}

## Dealing with NAs caused by missing molecules ------------------------------------------
# If certain atoms are missing in a chunk the column contains NAs
# Let's convert these NAs to zeros

#Columns with NA values
cols_with_NA <- colnames(new_df)[colSums(is.na(new_df)) > 0]

## Turning NAs into zeros
require(magrittr)
new_df <- new_df %>% 
  dplyr::mutate(dplyr::across(.cols = cols_with_NA, .fns = function(x) dplyr::if_else(is.na(x), 0, x)))


## Saving the final dataset for building the ML model ---------------
# Removing unnecessary molecule representation
new_df <- dplyr::select(new_df, -inchi)

write.csv(new_df, file = "data/final_dataset.csv", row.names = FALSE)
