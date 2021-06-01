#see present working directory, and set new one: 
	getwd() 
	setwd("C:/Users/Anna/OneDrive/Documents/University work/Project/workspace")

#remove any existing files: 
	if (file.exists("failvnonBarChart.pdf")) {file.remove("failvnonBarChart.pdf")} 
	if (file.exists("sequencingYears.pdf")) {file.remove("sequencingYears.pdf")} 	
	if (file.exists("sizeVyear.pdf")) {file.remove("sizeVyear.pdf")}
	if (file.exists("sizeVnumgene.pdf")) {file.remove("sizeVnumgene.pdf")}
	if (file.exists("numGenevsizeNT.pdf")) {file.remove("numGenevsizeNT")}
	if (file.exists("sizeVfail.xls")) {file.remove("sizeVfail.xls")}
	if (file.exists("sizeVpcfail.pdf")) {file.remove("sizeVpcfail.pdf")}
	if (file.exists("sizeVpcfail8F.pdf")) {file.remove("sizeVpcfail8F.pdf")}
	if (file.exists("sizeVpcfail1F.pdf")) {file.remove("sizeVpcfail1F.pdf")}
	if (file.exists("sizeVpcfail0.5F.pdf")) {file.remove("sizeVpcfail0.5F.pdf")}
	if (file.exists("numgeneVpcfail.pdf")) {file.remove("numgeneVpcfail.pdf")}
	if (file.exists("numgeneVpcfail0.5F.pdf")) {file.remove("numgeneVpcfail0.5F.pdf")}
	if (file.exists("instTable1.pdf")) {file.remove("instTable1.pdf")}
	if (file.exists("instTable2.pdf")) {file.remove("instTable2.pdf")}
	if (file.exists("mostFreqInst.csv")) {file.remove("mostFreqInst.csv")}
	if (file.exists("mostFailsInst.csv")) {file.remove("mostFailsInst.csv")}
	if (file.exists("allInstFreq.csv")) {file.remove("allInstFreq.csv")}
	if (file.exists("yearVpcfail.pdf")) {file.remove("yearVpcfail.pdf")}
	if (file.exists("yearVpcfail8F.pdf")) {file.remove("yearVpcfail8F.pdf")}
	if (file.exists("yearVfail.xls")) {file.remove("yearVfail.xls")}
	if (file.exists("genomeFailcats.pdf")) {file.remove("genomeFailcats.pdf")}


#set list of files in the pwd to 'filelist' variable; print list of files  
	filelist <- list.files(path = ".", pattern = ".csv$") 
	filelist 

#load in the file (genomes spreadsheet): 
	overviewData <- read.csv(file = "overview.csv", header = TRUE)

#save the data in each column to the column headings 
	accn <- unlist(overviewData$accn) 
	organism <- unlist(overviewData$organism) 
	sizeNT <- unlist(overviewData$sizeNT) 
	numGene <- unlist(overviewData$numGene)
	numFail <- unlist(overviewData$numFail) 
	pcFail <- unlist(overviewData$pcFail) 
	fail1 <- unlist(overviewData$fail1)
	fail2 <- unlist(overviewData$fail2)
	fail3 <- unlist(overviewData$fail3)
	fail4 <- unlist(overviewData$fail4)
	fail5 <- unlist(overviewData$fail5)
	fail <- unlist(overviewData$fail) 
	GCcontent <- unlist(overviewData$GCcontent)
	baseFreqA <- unlist(overviewData$baseFreqA)
	baseFreqT <- unlist(overviewData$baseFreqT)
	baseFreqC <- unlist(overviewData$baseFreqC)
	baseFreqG <- unlist(overviewData$baseFreqG)	
	institution <- unlist(overviewData$institution)
	seqYear <- unlist(overviewData$seqYear) 
	reassessed <- unlist(overviewData$reassessed) 



	
#################################
# failed and non-failed genomes #
#################################

#know we have 650 in total 
#to check, can count how many are 'reassessed' (because know that every genome has 1 for reassessed) 
table(reassessed) 

#count how many genomes are failed and how many are non-failed 
	countFails <- table(fail) 
#could save this as table: 
	# write.csv(countFails, file="numFailedGenomes.csv") 
	
#make bar chart and save as pdf: 
pdf("failvnonBarChart.pdf")
	colours <- c("green", "red")
	barplot(table(fail), main="Number of failed and non-failed genomes", xlab="Non-failed (0) or failed (1)", ylab="Number of genomes", ylim=c(0,400), col=colours) 
dev.off() 


#######################
# SEQUENCES OVER TIME #
#######################

#bar chart to represent the spread of sequences over years: 
yeardata <- table(seqYear)

pdf("sequencingYears.pdf", height=5, width=8) 
	barplot(yeardata, ylim=c(0,100), cex.axis=0.7, main="Distribution of sequences over time", xlab="Year of sequencing", ylab="Number of genomes", col="steelblue3")
dev.off()


################################
## genome size (nt) v seqYear ##
################################

#variables: 
	slope.yearVsizeNT <- signif(coef(lm(sizeNT~seqYear))[2], digits = 4) 
	cor.yearVsizeNT<- cor.test(seqYear, sizeNT) 
	r.yearVsizeNT <- signif(cor.yearVsizeNT$estimate, digits = 4) 
	Pval.yearVsizeNT <- signif(cor.yearVsizeNT$p.value, digits = 4) 

#plot: 	
pdf("sizeVyear.pdf") 
	plot(seqYear, sizeNT, xlab = "Year of sequencing", ylab = "Genome size (NT)", col = "orange", cex = 0.5, cex.axis = 0.8, pch = 16, main = "Sequence year v genome size (NT)")
	abline(lm(sizeNT~seqYear), col = "black") 
	mtext(paste0("gradient = ", slope.yearVsizeNT), side = 3, adj = 0.05, line = -1.2, cex = 0.8) 
	mtext(paste0("r = ", r.yearVsizeNT), side = 3, adj = 0.05, line = -2.2, cex = 0.8)
	mtext(paste0("p-value = ", Pval.yearVsizeNT), side = 3, adj = 0.05, line = -3.2, cex = 0.8)
dev.off()


############################
### SIZE (NT) V NUM GENE ###
############################	

#variables: 
	slope.sizeVnumgene <- signif(coef(lm(numGene~sizeNT))[2], digits = 4) 
	cor.sizeVnumgene <- cor.test(sizeNT, numGene) 
	r.sizeVnumgene <- signif(cor.sizeVnumgene$estimate, digits = 4) 
	Pval.sizeVnumgene <- signif(cor.sizeVnumgene$p.value, digits = 4) 
	
#plot: 
pdf("sizeVnumgene.pdf") 
plot(sizeNT, numGene, xlab = "Genome size (NT)", ylab = "No. genes in genome", col = "blue", cex = 0.5, pch = 18, main = "Genome size (NT) v no. genes in genome") 
	abline(lm(numGene~sizeNT), col = "black") 
	options(scipen=999)
		#scipen=999 converts x axis labels from scientific notation (e) to standard notation 
	mtext(paste0("gradient = ", slope.sizeVnumgene), side = 1, adj = 0.95, line = -3.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.sizeVnumgene), side = 1, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVnumgene), side = 1, adj = 0.95, line = -1.3, cex = 0.8, col = "black")
dev.off()
	

#cannot really say whether genome size (NT) or number of genes is the dependent or independent variable 
#(although likely sizeNT is independent (x) variable, and numGene is dependent (y) variable, if assume that dependent (y) is usually the more complex variable) 
#so likely that the above is the correct approach 
#however, do another plot with x and y variables switched: 

	slope.numgeneVsize <- signif(coef(lm(sizeNT~numGene))[2], digits = 4) 
	#correlation variables are same as above (order does not matter for cor.test, only for regression): 
	cor.sizeVnumgene <- cor.test(sizeNT, numGene) 
	r.sizeVnumgene <- signif(cor.sizeVnumgene$estimate, digits = 4) 
	Pval.sizeVnumgene <- signif(cor.sizeVnumgene$p.value, digits = 4) 
	
pdf("numGenevsizeNT.pdf") 
plot(numGene, sizeNT, xlab = "No. genes in genome", ylab = "Genome size (NT)", col = "blue", cex = 0.5, pch = 18, main = "No. genes in genome v genome size (NT)") 
	abline(lm(sizeNT~numGene), col = "black") 		
	mtext(paste0("gradient = ", slope.numgeneVsize), side = 1, adj = 0.95, line = -3.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.sizeVnumgene), side = 1, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVnumgene), side = 1, adj = 0.95, line = -1.3, cex = 0.8, col = "black")
dev.off()


########################
## size V fail or not ##
########################

#types of variables: continuous (size) v dichotomous (fail or not) 
#therefore stats test: point-biserial 

#correlation: 
	#syntax: ltm::biserial.cor(continuous, dichotomous, level=2) 
	ltm::biserial.cor(sizeNT, fail, level=2) 

#variables: 
	#point-biserial uses P value from Pearson's test: 
		cor.sizeVfail <- cor.test(sizeNT, fail) 
		Pval.sizeVfail <- signif(cor.sizeVfail$p.value, digits = 4) 
	#point-biserial uses cor coefficient (r pb, rather than r in Pearson's) 
		r.sizeVfail <- ltm::biserial.cor(sizeNT, fail, level=2) 
		r.sizeVfail <- format(r.sizeVfail, digits = 4) 
		#only one estimate given for biserial.cor test, therefore can assign the whole output to r value 

	
#save as table output: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "sizeVfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Genome size (NT) v fail", r.sizeVfail, Pval.sizeVfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "sizeVfail.xls", ncolumns=3, append=TRUE, sep="\t") 


######################
### SIZE V PC FAIL ###
###################### 

#variables: 
	slope.sizeVpcfail <- signif(coef(lm(pcFail~sizeNT))[2], digits = 4)  
	cor.sizeVpcfail <- cor.test(sizeNT, pcFail) 
	r.sizeVpcfail <- signif(cor.sizeVpcfail$estimate, digits = 4) 
	Pval.sizeVpcfail <- signif(cor.sizeVpcfail$p.value, digits = 4) 

#plot: 
pdf("sizeVpcfail.pdf")
plot(sizeNT, pcFail, xlab = "Genome size (NT)", ylab = "% of genes failed in genome", col = "violetred3", cex=0.5, pch=18, main = "Genome size v % genes failed")
	abline(lm(pcFail~sizeNT), col = "black")  
	mtext(paste0("gradient = ", slope.sizeVpcfail), side = 3, adj = 0.95, line = -1.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.sizeVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")
dev.off() 


##....LOOKING AT 8% FAIL AND UNDER....##

#focus on distribution of points under 8% fail (pdf output has '8F' at end for "8% focus") 
#variables: 
	slope.sizeVpcfail <- signif(coef(lm(pcFail~sizeNT))[2], digits = 4) 
	cor.sizeVpcfail <- cor.test(sizeNT, pcFail) 
	r.sizeVpcfail <- signif(cor.sizeVpcfail$estimate, digits = 4) 
	Pval.sizeVpcfail <- signif(cor.sizeVpcfail$p.value, digits = 4) 

#plot graph: 
pdf("sizeVpcfail8F.pdf")
plot(sizeNT, pcFail, xlab = "Genome size (NT)", ylab = "% of genes failed in genome", ylim = c(0,8), col = "violetred3", cex=0.5, pch=18, main = "Genome size v % genes failed (8% and under)")
	#abline(lm(pcFail~sizeNT), col = "black") 
		#not plotting linear regression for the 8% focus graph 
		#(because calculated using the parameters of the full-scale set of data, therefore not representative of smaller scale)   
	#mtext(paste0("gradient = ", slope.sizeVpcfail), side = 3, adj = 0.95, line = -1.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.sizeVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")
dev.off() 
		

##....LOOKING AT 1% FAIL AND UNDER....##

#to see if year of sequencing (older or more recent) is related to the distribution of the points

	overviewData$colour="black" 
	yearCol <- overviewData$colour 
	overviewData$colour[overviewData$seqYear==2003 | overviewData$seqYear==2004 | overviewData$seqYear==2005 | overviewData$seqYear==2006 | overviewData$seqYear==2007 | overviewData$seqYear==2008 | overviewData$seqYear==2009]=4 
	overviewData$colour[overviewData$seqYear==2010 | overviewData$seqYear==2011 | overviewData$seqYear==2012 | overviewData$seqYear==2013 | overviewData$seqYear==2014 | overviewData$seqYear==2015]=2 
	#assigns colours (4=blue, 2=red) to the years: 2003-2009 are coded blue(4), and 2010-2015 are coded red(2) 
	#this will allow us to visualise the spread of points by year sequenced (e.g. may expect most of the small genome & high % fail values to be older sequences) 

pdf("sizeVpcfail1F.pdf")
plot(sizeNT, pcFail, xlab = "Genome size (NT)", ylab = "% of genes failed in genome", ylim = c(0,1), col = overviewData$colour, cex=0.5, pch=16, main = "Genome size v % genes failed (1% and under)")
	mtext(paste0("r = ", r.sizeVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")		
dev.off() 


##....LOOKING AT 0.5% FAIL AND UNDER....##

pdf("sizeVpcfail0.5F.pdf")
plot(sizeNT, pcFail, xlab = "Genome size (NT)", ylab = "% of genes failed in genome", ylim = c(0,0.5), col = overviewData$colour, cex=0.5, pch=16, main = "Genome size v % genes failed (0.5% and under)")
	mtext(paste0("r = ", r.sizeVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.sizeVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")		
dev.off()


##########################
### NUM GENE V PC FAIL ###
##########################

#numGene v pcFail: should be very similar to 'size V pcFail' graph 
#variables: 
	slope.numgeneVpcfail <- signif(coef(lm(pcFail~numGene))[2], digits = 4) 
	cor.numgeneVpcfail <- cor.test(numGene, pcFail) 
	r.numgeneVpcfail <- signif(cor.numgeneVpcfail$estimate, digits = 4) 
	Pval.numgeneVpcfail <- signif(cor.numgeneVpcfail$p.value, digits = 4) 

#plot: 
pdf("numgeneVpcfail.pdf") 
	plot(numGene, pcFail, xlab = "No. genes in genome", ylab = "% of genes failed in genome", col = "purple4", cex=0.5, pch=18, main = "No. genes in genome v % genes failed")
	abline(lm(pcFail~numGene), col = "black") 
	mtext(paste0("gradient = ", slope.numgeneVpcfail), side = 3, adj = 0.95, line = -1.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.numgeneVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.numgeneVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")
dev.off()


##....LOOKING AT 0.5% FAIL AND UNDER....##

pdf("numgeneVpcfail0.5F.pdf") 
plot(numGene, pcFail, xlab = "No. genes in genome", ylab = "% of genes failed in genome", ylim = c(0,0.5), col = "purple4", cex=0.5, pch=18, main = "No. genes in genome v % genes failed (0.5% and under)")
	#don't want the lm line calculated from the full-scale set of data - not representative 
		#abline(lm(pcFail~numGene), col = "black") 
		#mtext(paste0("gradient = ", slope.numgeneVpcfail), side = 3, adj = 0.95, line = -1.3, cex = 0.8, col = "black") 
	mtext(paste0("r = ", r.numgeneVpcfail), side = 3, adj = 0.95, line = -2.3, cex = 0.8, col = "black")
	mtext(paste0("p-value = ", Pval.numgeneVpcfail), side = 3, adj = 0.95, line = -3.3, cex = 0.8, col = "black")
dev.off()


#############################
### INSTITUTION V PC FAIL ###
#############################

#NOTE: this section was not included in analysis due to too many problems encountered, however could be interesting route in the future 

#need to clean up institution name values: 
	#some of the top institutes with the most genomes sequenced are duplicate names of the same place - have minor differences in the address name  
	#want to put them all under the same name - have to do this manually, but important because they are significant institutions 
	overviewData$newInst=overviewData$institution 
	overviewData$newInst[overviewData$institution=="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B100_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
	overviewData$newInst[overviewData$institution=="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_B310_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
	overviewData$newInst[overviewData$institution=="_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_"]="_US_DOE_Joint_Genome_Institute_2800_Mitchell_Drive_Walnut_Creek_CA_94598-1698_USA_" 
	overviewData$newInst[overviewData$institution=="_The_Institute_for_Genomic_Research_9712_Medical_Center_Dr_Rockville_MD_20850_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_" 
	overviewData$newInst[overviewData$institution=="_The_J_Craig_Venter_Institute_9704_Medical_Center_Dr_Rockville_MD_20850_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_"
	overviewData$newInst[overviewData$institution=="_The_J_Craig_Venter_Institute_Rockville_MD_USA_"]="_J_Craig_Venter_Institute_9704_Medical_Center_Drive_Rockville_MD_20850_USA_"

	institution <- overviewData$newInst 

	#replace underscores with spaces:  
	institution <- gsub("_", " ", institution) 

summary(institution) 


#note that I only did the above address values because they showed up in the 'most frequently appearing institutions' top 6 
#therefore there is a chance that there are other spellings/differently-formatted versions of the same institution elsewhere with 1 or 2 
#and as these ones are not being picked up, they are not being added to the total quantity for that institution 

#even though not able to sort through the institutions properly because of problems with consistency of addresses used, 
	#this could perhaps be an indicator of 'quality' of submitted sequences in itself - institutions are not using consistent addressed, therefore not being grouped correctly 
	#could go down this direction as well! e.g. 'address inconsistency as indicator of quality' 
	#e.g. could do a table of institutions and the number of 'alternative/inconsistent' addresses they used 
	#could be used either by the instituions to let them know their submitted addresses are being inconsistent 
		#although this may be no fault of theirs - it is the contact address, so maybe the address has simply changed over time, or has merged with another institution 
		#or may be from 2 different departments within a same university, and important to keep distinction? 
	#in which case, could maybe be used by the public database (e.g. EBI) itself - to let them know there should be an extra field in the annotation that uses a streamlined identifier for the institution 
	#or a field that specifies the city? so they can be grouped together? etc 

#some stats we can pull out: 

	institution[1:5] #pulls out first 5 entries 
	institution[which(overviewData$pcFail > 0)] #pulls out institutions which have fail percentage greater than 0% (i.e. for all genomes which have failed at all) 
	institution[which(overviewData$pcFail > 20)] #"    " greater than 20%

	#subset function: pull out data for which pcFail is over 20%, and limit the columns to what is relevant (select=c; can also reorder): 
	worstInst <- subset(overviewData, pcFail > 20, select=c(newInst, seqYear, pcFail, fail1, fail2, fail3, fail4, fail5, numGene, numFail)) 
	worstInst$newInst <- gsub("_", " ", worstInst$newInst)
		#(cleans up the underscores again in the 'newInst' column) 

	#'sort' function: 
		sort(worstInst$seqYear) #puts out values in the seqYear column, in order, but independent of data frame 
	#'order' function to sort through data WITHIN its data frame (dependent on other variables) 
		date1.worstInst <- worstInst[order(worstInst$seqYear), ] #sorts all of data frame, from first (smallest, earliest) year at top and latest at bottom 
		date2.worstInst <- worstInst[order(-worstInst$seqYear), ] #from largest (most recent) year to smallest 
		bypc.worstInst <- worstInst[order(-worstInst$pcFail), ] #from highest pcFail to lowest (bypc = by percentage (fail)) 

#now actually save these data: 
		
#table of institutes responsible for the 6 'worst' genome sequences, in order year (2002 to 2015) 
	date1.worstInst <- worstInst[order(worstInst$seqYear), ] 
	write.csv(date1.worstInst, file="instTable1.csv") 
	
#table of institutes responsible for the 6 'worst' genome sequences, in order pcFail (high to low) 
	bypc.worstInst <- worstInst[order(-worstInst$pcFail), ] 
	write.csv(bypc.worstInst, file="instTable2.csv") 


#############################
# INSTITUTION V FAIL AT ALL #
#############################

#subset table of all (as all have been reassessed); then summary to sum up the totals 
#i.e. table of top 6 most frequent institutions (with number of how many times they appear, i.e. how many sequenced genomes uploaded) 
	mostFreq <- subset(overviewData, reassessed==1, select=c(newInst, reassessed)) 
	#says make a new table (vector) 'mostFreq', made of a subset of the overviewData table, composed of only values of reassessed==1, and only 2 columns 
	mostFreq6 <- summary(mostFreq) 
	write.csv(mostFreq6, file="mostFreqInst.csv") 
	#***NOTE: output is not in standard tabular format, therefore need to reformat manually for report 
	
	worstInst2 <- subset(overviewData, fail==1, select=c(newInst, fail)) 
	mostFails <- summary(worstInst2) 
	write.csv(mostFails, file="mostFailsInst.csv") 
		#this is just for the failed genomes 
		#will give the number of failed genomes from each institution - can manually check with table before? somehow do correlation test? 
	
#could do bar chart for the top 6 institutes (in total genomes sequenced) which has one bar for total genomes frequency (reassessed=1) and one bar for failed genomes frequency (fail=1) 
	#or choose top 6 institutes based on highest pc fails sum? (sum of all their pc fails) 

#or manually make a table of the 6 most frequent ones, with total number of genomes sequenced, and then number of failed genomes (so 3 columns: institute name, total genomes, failed genomes) 
#'grouped bar plot': 
#however this is difficult to do - would be lengthy manual process of summing up all the totals by the institutions, because of the inconsistencies in address 
#(discussed this above) 


allInstFreq <- table(mostFreq)
write.csv(allInstFreq, file="allInstFreq.csv") #(all institute frequencies) 


	numTotalGenomes <- c(1st inst, 2nd inst, 3rd inst etc) 
	numFailedGenomes <- c(1st inst, 2nd inst, 3rd inst etc) 
	counts <- table(numTotalGenomes, numFailedGenomes) 
	barplot(counts, main=" ", xlab=" ", col=c(" "," "), legend = rownames("Total number of genomes sequenced", "Number of failed genomes sequenced"), beside=TRUE)  
	



###############################
# institution + num fail cats #
###############################

#overviewGenes has column for num failcats (from 0-5 hypothetically, but only max 4 because no sequences failed all 5 cats) 
#make similar column for genome analysis? 
#would need to count if greater than 0 for each fail column, then add new column saying if fail category has been met, then sum up total number of fail categories 
#this is possible in R but probably better done in Tcl - probably not worth doing it for this project 
#include in Discussion section for future directions of this research 
	#would use to see which institutions are involved in genomes of 4 fail categories 
	#and count how many genomes are in each numfailcats category (how many in 0, how many in 1, etc.) 


#for Discussion section in report: do more stats on institutions, e.g. calculate percentage of failed genomes over total genomes sequenced for the 'big' institutions 

######################
### YEAR V PC FAIL ###
######################

#variables: 
	slope.yearVpcfail <- signif(coef(lm(pcFail~seqYear))[2], digits = 4) 
	cor.yearVpcfail <- cor.test(seqYear, pcFail) 
	r.yearVpcfail <- signif(cor.yearVpcfail$estimate, digits = 4) 
	Pval.yearVpcfail <- signif(cor.yearVpcfail$p.value, digits = 4) 
	
#plot graph: 
pdf("yearVpcfail.pdf")
plot(seqYear, pcFail, xlab = "Year of sequencing", ylab = "% of genes failed in genome", col = "darkmagenta", cex = 0.5, cex.axis = 0.9, pch = 18, main = "Genome sequence year v % of genes failed")
	abline(lm(pcFail~seqYear), col = "black") 
	mtext(paste0("gradient = ", slope.yearVpcfail), side = 3, adj = 0.95, line = -2.5, cex = 0.8) 
	mtext(paste0("r = ", r.yearVpcfail), side = 3, adj = 0.95, line = -3.5, cex = 0.8)
	mtext(paste0("p-value = ", Pval.yearVpcfail), side = 3, adj = 0.95, line = -4.5, cex = 0.8)
dev.off() 	


##....FOCUS ON UNDER 8% OF GENES FAILED....## 
#zoom in on under 8% of genes failed (where most points are concentrated) to see trend closer up: 
pdf("yearVpcfail8F.pdf") 
plot(seqYear, pcFail, xlab = "Year of sequencing", ylab = "% of genes failed in genome", col = "darkmagenta", cex = 0.5, cex.axis = 0.9, ylim = c(0,8), pch = 18, main = "Genome sequence year v % of genes failed (8% and under)")
	mtext(paste0("r = ", r.yearVpcfail), side = 3, adj = 0.95, line = -1.5, cex = 0.8)
	mtext(paste0("p-value = ", Pval.yearVpcfail), side = 3, adj = 0.95, line = -2.5, cex = 0.8)
dev.off() 


################################
# YEAR V FAIL OR NOT (GENOMES) #
################################

#variables: discrete numeric (year) + dichotomous (fail or not) 
#no point in doing scatter graph because y is dichotomous/binary rather than a continuous variable 
#therefore stat test: point-biserial 

#ltm::biserial.cor(continuous, dichotomous, level=2) 
ltm::biserial.cor(seqYear, fail, level=2) 

#variables: 
	cor.yearVfail <- cor.test(seqYear, fail) 
	Pval.yearVfail <- signif(cor.yearVfail$p.value, digits = 4) 
	r.yearVfail <- ltm::biserial.cor(seqYear, fail, level=2) 
		r.yearVfail <- format(r.yearVfail, digits = 4) 
		#only one estimate given for biserial.cor test, therefore can assign the whole output to r value 
		#note: r is 'r pb' coefficient in point-biserial 
	

#save output of correlation test variables: 
topline <- c("variables", "r", "p-value") 
out1 <- print(topline, quote=FALSE) 
write(out1, file = "yearVfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Year of sequencing genome v fail", r.yearVfail, Pval.yearVfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "yearVfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	
#######################################
# number of genomes with each failcat #
#######################################

#first, count number of failed genomes in each fail category: 
	#fail1 genomes: 
		#table(fail1) to return table of how many genomes have each number of fail1 genes (e.g. 0 = 604, 1 = 17, 352 = 1, etc.) 
		#(i.e. 604 genomes had 0 fail 1 genes, 17 genomes had just 1 fail 1 gene, 1 genome had 352 failed genes) 
		#but just want to know how many genomes in total had 1 or more fail 1 gene, don't care HOW many genes within each genome there were 
		
		fail1s <- subset(overviewData, fail1 > 0, select=c(reassessed, sizeNT, numGene, numFail, fail1)) 
			#subset table just containing select columns of interest 
		fail1sNUM <- table(fail1s$reassessed) 
			#gives us way of counting how many genomes have a fail1 category classification 
			#uses 'fail1s' table and 'reassessed' column to count how many (because all have same 'reassessed' value, i.e. 1, therefore give number of total rows) 
			#saves to variable 'fail1sNUM' which will use for constructing the bar chart 
			
	#fail2 genomes: 
		fail2s <- subset(overviewData, fail2 > 0, select=c(reassessed, sizeNT, numGene, numFail, fail2))
		fail2sNUM <- table(fail2s$reassessed) 
		
	#fail3 genomes: 
		fail3s <- subset(overviewData, fail3 > 0, select=c(reassessed, sizeNT, numGene, numFail, fail3))
		fail3sNUM <- table(fail3s$reassessed) 
	
	#fail4 genomes: 
		fail4s <- subset(overviewData, fail4 > 0, select=c(reassessed, sizeNT, numGene, numFail, fail4))
		fail4sNUM <- table(fail4s$reassessed) 
		
	#fail5 genomes: 
		fail5s <- subset(overviewData, fail5 > 0, select=c(reassessed, sizeNT, numGene, numFail, fail5))
		fail5sNUM <- table(fail5s$reassessed) 
		
#bar chart of these values to visualise how many genomes fell into each failcat: 
	#vector of values and colours: 
		failCounts <- c(fail1sNUM, fail2sNUM, fail3sNUM, fail4sNUM, fail5sNUM) 
		colours <- c("brown3", "chocolate1", "chartreuse3", "darkorchid", "aquamarine4") 
	#plot: 
	pdf("genomeFailcats.pdf")
		barplot(failCounts, main="Number of genomes classified with each fail category", xlab="Fail categories", ylab="Number of genomes", ylim=c(0,350), cex.axis=1, cex.name=1, col=colours, names.arg=c("fail 1","fail 2","fail 3","fail 4","fail 5")) 
	dev.off() 


#if enough time - sort out institution values, perhaps as further inidication of quality 
#also - could count number of failcats for each gene/genome and do analysis on that level (see if have same pattern of correlation as with instance/% of fail 
#(note: chose not to analyse the above) 