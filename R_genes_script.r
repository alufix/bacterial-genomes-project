#this script is for gene-level analysis 
#importing overviewGenes file, however when saving to variables, using letter 'G' as prefix 
#this is to make sure the gene-level variables are kept distinct to genome-level variables (in case decide to import and include this data) 


#set new working directory 
	getwd() 	
	setwd("C:/Users/Anna/OneDrive/Documents/University work/Project/workspace")

#remove any existing files: 
	if (file.exists("genetypesBarChart1.pdf")) {file.remove("genetypesBarChart1.pdf")} 
	if (file.exists("genetypesBarChart2.pdf")) {file.remove("genetypesBarChart2.pdf")} 
	if (file.exists("GcomplVGfail.xls")) {file.remove("GcomplVGfail.xls")} 
	if (file.exists("GcomplVGnumfail.xls")) {file.remove("GcomplVGnumfail.xls")} 
	if (file.exists("GpseudoVGfail.xls")) {file.remove("GpseudoVGfail.xls")} 
	if (file.exists("GpseudoVGnumfail.xls")) {file.remove("GpseudoVGnumfail.xls")} 
	if (file.exists("GnonspecVGcompl.xls")) {file.remove("GnonspecVGcompl.xls")} 
	if (file.exists("GnonspecVGpseudo.xls")) {file.remove("GnonspecVGpseudo.xls")} 
	if (file.exists("GnonspecVGfail.xls")) {file.remove("GnonspecVGfail.xls")} 
	if (file.exists("GnonspecVGnumfail.xls")) {file.remove("GnonspecVGnumfail.xls")} 
	
	
#install relevant R packages -  
#ltm package:  
	install.packages("ltm")
	#(package for point-biserial correlation test) 
	

#set list of the files in pwd to 'filelist' variable, and print the list of files 
	filelist <- list.files(path = ".", pattern = ".csv$") 
	filelist 

#load in overviewGenes 
	overviewGenes <- read.csv(file = "overviewGenes.csv", header = TRUE)

#save in the values in each column 
#(put 'G' for 'gene' before headings to distinguish them from the overview (genomes) file headings) 
		
		Gaccn <- unlist(overviewGenes$accn) 
		Gorganism <- unlist(overviewGenes$organism)
		Gloctag <- unlist(overviewGenes$loctag) 
		Gexact <- unlist(overviewGenes$exact)
		Gname <- unlist(overviewGenes$name)
		Gcompl <- unlist(overviewGenes$compl)
		Gpseudo <- unlist(overviewGenes$pseudo)
		Gfail <- unlist(overviewGenes$fail) #(fail or not) 
		GnumFailcats <- unlist(overviewGenes$numFailcats) #(number of fail categories) 
		Gfail1 <- unlist(overviewGenes$fail1)
		Gfail2 <- unlist(overviewGenes$fail2)
		Gfail3 <- unlist(overviewGenes$fail3)
		Gfail4 <- unlist(overviewGenes$fail4)
		Gfail5 <- unlist(overviewGenes$fail5)
		Ginstitution <- unlist(overviewGenes$institution)
		GseqYear <- unlist(overviewGenes$seqYear) 


#change values of Gexact to Gnonspec (for partial CDS variable) 
	#e.g. for correlating Gexact and Gcompl - cor.test(Gexact, Gcompl) doesn't work because 'Gexact' doesn't consist of numeric values 
	#need to convert 'exact' to '0' and 'non-specific' to '1' 
	
	Gexact = gsub("exact", 0, Gexact) 
	Gexact = gsub("non-specific", 1, Gexact) 
	
	#then convert this to numeric values (so in numeric mode rather than character mode) 
	#also rename the variable Gnonspec so that it works with being TRUE/1 
	
	Gnonspec <- as.numeric(Gexact)
	summary(Gnonspec)
		#note: this doesn't change the value of overviewGenes$Gexact - this is the csv column 



###########################
# bar chart of gene types #
###########################

#table() builds contingency table for each type: Gcompl, Gpseudo etc. 	
	numCompl <- table(Gcompl)[c(2)] #saves column 2 (the number of 1/positive for compl) 	
	numPseudo <- table(Gpseudo)[c(2)]	
	numNonspec <- table(Gnonspec)[c(2)]
	numFailed <- table(Gfail)[c(2)] 

#make graph (by first making vector consisting of above values) 
pdf("genetypesBarChart1.pdf") 
	Gtypes <- c(numCompl, numPseudo, numNonspec, numFailed) 
	colours <- c("firebrick3", "darkorange2", "dodgerblue3", "limegreen")
	barplot(Gtypes, main="Number of genes in each 'gene type' category", xlab="Types of gene", ylab="Number of genes", ylim=c(0,1200000), cex.axis=0.7, cex.name=0.7, col=colours, names.arg=c("complement","pseudogene","non-specific pos","failed")) 
	options(scipen=999)
dev.off() 

#have many more complement than other types - therefore do another graph ignoring compl 
#just to look at proportion of other types to one another 
pdf("genetypesBarChart2.pdf") 
	Gtypes <- c(numPseudo, numNonspec, numFailed) 
	colours <- c("darkorange2", "dodgerblue3", "limegreen")
	barplot(Gtypes, main="Number of genes in each 'gene type' category (excl. complement)", xlab="Types of gene", ylab="Number of genes", ylim=c(0,12000), cex.main=0.85, cex.axis=0.8, cex.name=0.75, col=colours, names.arg=c("pseudogene","non-specific pos","failed")) 
dev.off() 	


#########################
## compl v fail or not ##
#########################

#dichotomous (compl) + dichotomous (fail or not) 
#therefore stats test: Pearson's (because same cor coefficient as Phi, but also provides P value and is consistent)
cor.test(Gcompl, Gfail) 

	cor.complVfail <- cor.test(Gcompl, Gfail) 
	r.complVfail <- signif(cor.complVfail$estimate, digits = 4) 
	Pval.complVfail <- signif(cor.complVfail$p.value, digits = 4)  
		Pval.complVfail <- format(Pval.complVfail, scientific=FALSE) 

#to make a contingency table: 
	chi.complVfail <- chisq.test(Gcompl, Gfail) 
	cont.complVfail <- chi.complVfail$observed
	#table output doesn't retain formatting, therefore manually recreate contingency table using these values (along with totals of columns/rows) 
	#and manually record 

	
#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GcomplVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Complement v fail", r.complVfail, Pval.complVfail)  
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GcomplVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

#and save contingency table 
	cont.complVfail 

	
##########################
## compl v num failcats ##
##########################

#dichotomous (compl) + numeric (discrete; num failcats) 
#therefore stats test: point-biserial (but do Pearson's cor.test for getting P-value) 
	#syntax: ltm::biserial.cor(continuous, dichotomous, level=2) 
	ltm::biserial.cor(GnumFailcats, Gcompl, level=2) 

#summary of number of fails for each gene: 
	summary(GnumFailcats) 
	#tells us max is 4 - therefore no single gene has failed all 5 categories 

#variables: 
	#P value from Pearson's cor test: 
		cor.complVnumfail <- cor.test(Gcompl, GnumFailcats) 
		Pval.complVnumfail <- signif(cor.complVnumfail$p.value, digits = 4) 
		Pval.complVnumfail <- format(Pval.complVnumfail, scientific=FALSE) 
	#r pb value from point-biserial test: 
		r.complVnumfail <- ltm::biserial.cor(GnumFailcats, Gcompl, level=2) 
		r.complVnumfail <- format(r.complVnumfail, digits = 4) 
		#(test only gives one 'estimate', therefore can assign whole outcome to 'r') 


#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GcomplVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Complement v no. fail categories", r.complVnumfail, Pval.complVnumfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GcomplVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 


##########################
## pseudo v fail or not ##
##########################
	
#dichotomous (pseudo) + dichotomous (fail or not) 
#therefore use Pearson's (same coefficient as Phi which is for binary variables) 

#variables: 
	cor.pseudoVfail <- cor.test(Gpseudo, Gfail) 
	r.pseudoVfail <- signif(cor.pseudoVfail$estimate, digits = 4) 
	Pval.pseudoVfail <- signif(cor.pseudoVfail$p.value, digits = 4) 

#contingency table: 
	chi.pseudoVfail <- chisq.test(Gpseudo, Gfail) 
	cont.pseudoVfail <- chi.pseudoVfail$observed 

#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GpseudoVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Pseudogene v fail", r.pseudoVfail, Pval.pseudoVfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GpseudoVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

#and save contingency table: 
	cont.pseudoVfail 

	
###########################
## pseudo v num failcats ##
###########################

#dichotomous (pseudo) + numeric (num failcats) 
#therefore use point-biserial 

#variables: 
	cor.pseudoVnumfail <- cor.test(Gpseudo, GnumFailcats) 
	Pval.pseudoVnumfail <- signif(cor.pseudoVnumfail$p.value, digits = 4) 
	r.pseudoVnumfail <- ltm::biserial.cor(GnumFailcats, Gpseudo, level=2) 
		r.pseudoVnumfail <- format(r.pseudoVnumfail, digits = 4) 

	
#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GpseudoVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Pseudogene v no. fail categories", r.pseudoVnumfail, Pval.pseudoVnumfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GpseudoVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	

######################
### pseudo v compl ###
###################### 

#variables: dichotomous (pseudo) v dichotomous (compl) 
#therefore use Pearson's (same r value as Phi coefficient for binary variables) 

#variables: 
	cor.pseudoVcompl <- cor.test(Gpseudo, Gcompl) 
	r.pseudoVcompl <- signif(cor.pseudoVcompl$estimate, digits = 4) 
	Pval.pseudoVcompl <- signif(cor.pseudoVcompl$p.value, digits = 4) 
		Pval.pseudoVcompl <- format(Pval.pseudoVcompl, scientific=FALSE) 

#contingency table: 
	chi.pseudoVcompl <- chisq.test(Gpseudo, Gcompl) 
	cont.pseudoVcompl <- chi.pseudoVcompl$observed 

#save output: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GpseudoVGcompl.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Pseudogene v complement", r.pseudoVcompl, Pval.pseudoVcompl) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GpseudoVGcompl.xls", ncolumns=3, append=TRUE, sep="\t") 
	
#and manually save contingency table: 
	cont.pseudoVcompl 

		
##########################
## non-specific v compl ##
########################## 

#variables: dichotomous (non-specific) v dichotomous (compl) 
#therefore use Pearson's (same r value as Phi coefficient for binary variables) 

#variables: 
	cor.nonspecVcompl <- cor.test(Gnonspec, Gcompl) 
	r.nonspecVcompl <- signif(cor.nonspecVcompl$estimate, digits = 4) 
	Pval.nonspecVcompl <- signif(cor.nonspecVcompl$p.value, digits = 4) 
		Pval.nonspecVcompl <- format(Pval.nonspecVcompl, scientific=FALSE) 

#contingency table: 
	chi.nonspecVcompl <- chisq.test(Gnonspec, Gcompl) 
	cont.nonspecVcompl <- chi.nonspecVcompl$observed 

#save output: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GnonspecVGcompl.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Non-specific position v complement", r.nonspecVcompl, Pval.nonspecVcompl) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GnonspecVGcompl.xls", ncolumns=3, append=TRUE, sep="\t") 
	
#and manually save contingency table: 
	cont.nonspecVcompl 


###########################
## non-specific v pseudo ##
########################### 

#variables: dichotomous (non-spec) + dichotomous (pseudo) 
#therefore use Pearson's 

#variables: 
	cor.nonspecVpseudo <- cor.test(Gnonspec, Gpseudo) 
	r.nonspecVpseudo <- signif(cor.nonspecVpseudo$estimate, digits = 4) 
	Pval.nonspecVpseudo <- signif(cor.nonspecVpseudo$p.value, digits = 4) 
		Pval.nonspecVpseudo <- format(Pval.nonspecVpseudo, scientific=FALSE) 	
	
#contingency table: 
	chi.nonspecVpseudo <- chisq.test(Gnonspec, Gpseudo) 
	cont.nonspecVpseudo <- cor.nonspecVpseudo$observed 

#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GnonspecVGpseudo.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Non-specific position v pseudogene", r.nonspecVpseudo, Pval.nonspecVpseudo) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GnonspecVGpseudo.xls", ncolumns=3, append=TRUE, sep="\t") 

#and save contingency table manually: 
	cont.nonspecVpseudo 


#########################
## non-specific v fail ##
######################### 	

#variables: dichotomous (nonspec) + dichotomous (fail) 
#therefore use Pearson's 

#variables: 
	cor.nonspecVfail <- cor.test(Gnonspec, Gfail) 
	r.nonspecVfail <- signif(cor.nonspecVfail$estimate, digits = 4) 
	Pval.nonspecVfail <- signif(cor.nonspecVfail$p.value, digits = 4) 
		Pval.nonspecVfail <- format(Pval.nonspecVfail, scientific=FALSE) 
	
#contingency table: 
	chi.nonspecVfail <- chisq.test(Gnonspec, Gfail) 
	cont.nonspecVfail <- cor.nonspecVfail$observed 	
	
#save output: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GnonspecVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Non-specific position v fail", r.nonspecVfail, Pval.nonspecVfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GnonspecVGfail.xls", ncolumns=3, append=TRUE, sep="\t") 

#also manually save contingency table: 
	cont.nonspecVfail 

		
#################################
## non-specific v num failcats ##
################################# 

#variables: dichotomous (nonspec) + discrete numeric (num failcats) 
#therefore use point-biserial 

#variables: 
	cor.nonspecVnumfail <- cor.test(Gnonspec, GnumFailcats) 
	Pval.nonspecVnumfail <- signif(cor.nonspecVnumfail$p.value, digits = 4) 
	#and get cor value (r pb) from point-biserial test: 
	r.nonspecVnumfail <- ltm::biserial.cor(GnumFailcats, Gnonspec, level=2) 
		r.nonspecVnumfail <- format(r.nonspecVnumfail, digits = 4) 

#save output:  
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GnonspecVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Non-specific position v no. fail categories", r.nonspecVnumfail, Pval.nonspecVnumfail) 
	out2 <- print(line2, quote=FALSE) 
	write(out2, file = "GnonspecVGnumfail.xls", ncolumns=3, append=TRUE, sep="\t") 

