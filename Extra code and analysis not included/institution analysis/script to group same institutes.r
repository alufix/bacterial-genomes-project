#R script for grouping addresses for same institution but with different format 


_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_                   
_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B100_Walnut_Creek_CA_94598-1698_USA_     
_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B310_Walnut_Creek_CA_94598-1698_USA_         
_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_

#create new column of institute names, populate it with the existing institute names 
#then for the duplicates (2nd, 3rd + 4th in above list), make it the same as the 1st 

#have same issue with 
_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_
_The_Institute_for_Genomic_Research_9712_Medical_Center_Dr_Rockville_MD_20850_USA_
_The_J_Craig_Venter_Institute_9704_Medical_Center_Dr_Rockville_MD_20850_USA_
_The_J_Craig_Venter_Institute_Rockville_MD_USA_

#the 2nd is a 'precursor' to the 1st one (https://en.wikipedia.org/wiki/J._Craig_Venter_Institute) 
#and the 3rd + 4th are a different format but same address as the 1st 
#so replace these values with the 1st 


overviewData$newInst=overviewData$institution 
overviewData$newInst[overviewData$institution=="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B100_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
overviewData$newInst[overviewData$institution=="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B310_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
overviewData$newInst[overviewData$institution=="_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
overviewData$newInst[overviewData$institution=="_The_Institute_for_Genomic_Research_9712_Medical_Center_Dr_Rockville_MD_20850_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_" 
overviewData$newInst[overviewData$institution=="_The_J_Craig_Venter_Institute_9704_Medical_Center_Dr_Rockville_MD_20850_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_"
overviewData$newInst[overviewData$institution=="_The_J_Craig_Venter_Institute_Rockville_MD_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_"

institution <- overviewData$newInst 

	mostFreq <- subset(overviewData, reassessed = 1, select=c(newInst, reassessed)) 
	mostFreq6 <- summary(mostFreq) 
	write.csv(mostFreq6, file="mostFreqInst.csv") 
