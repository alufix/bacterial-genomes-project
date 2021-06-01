#to run after 'R_genes_script.r', therefore gene-level overview document will already be loaded, thus contents and packages/functions will already be saved 


#remove any existing files: 
	if (file.exists("pcFailcatsBarChart.pdf")) {file.remove("pcFailcatsBarChart.pdf")}
	if (file.exists("numFailcatsBarChart.pdf")) {file.remove("numFailcatsBarChart.pdf")}
	if (file.exists("GnonspecVGallcats.xls")) {file.remove("GnonspecVGallcats.xls")}
	if (file.exists("GcomplVGallcats.xls")) {file.remove("GcomplVGallcats.xls")}
	if (file.exists("GpseudoVGallcats.xls")) {file.remove("GpseudoVGallcats.xls")}


     ###########################
*****  SPECIFIC FAIL CATEGORIES *****
     ###########################

###################################
# bar chart of all failcat totals # 
###################################

#sum total number of genes: 
	totalGenes <- table(Gfail1)[c(1)] + table(Gfail1)[c(2)] 
	#('totalGenes' is sum of fail1 = 0 and fail1 = 1, therefore way of counting all genes) 
	
#table() builds contingency table 
#so for table of fail1 values: table(Gfail1) ; etc 
	f1num <- table(Gfail1)[c(2)] #saves column 2 (the number of 1/positive for fail1 values) 
	f1pc <- f1num/totalGenes * 100 #percentage of f1 genes out of total num genes 

	f2num <- table(Gfail2)[c(2)] 
	f2pc <- f2num/totalGenes * 100 

	f3num <- table(Gfail3)[c(2)]
	f3pc <- f3num/totalGenes * 100 

	f4num <- table(Gfail4)[c(2)]
	f4pc <- f4num/totalGenes * 100 

	f5num <- table(Gfail5)[c(2)] 
	f5pc <- f5num/totalGenes * 100 
	
#make vector of these pc values to plot, and make graph 
pdf("pcFailcatsBarChart.pdf") 
	pcF <- c(f1pc, f2pc, f3pc, f4pc, f5pc) 
	colours <- c("brown3", "chocolate1", "chartreuse3", "darkorchid", "aquamarine4")
	barplot(pcF, main="Number of genes within each fail category, as percentage of total genes", xlab="Fail categories", ylab="% of genes failed", cex.main=0.9, ylim=c(0,0.25), col=colours, names.arg=c("fail 1","fail 2","fail 3","fail 4","fail 5")) 
dev.off() 


#if want to have bar plot showing the proportions of gene fail categories, but as numbers rather than percentages: 
#(will show same proportions) 
pdf("numFailcatsBarChart.pdf")
	fALLnum <- c(f1num, f2num, f3num, f4num, f5num) 
	barplot(fALLnum, col=colours, main="Number of genes within each fail category", xlab="Fail categories", ylab="Number of genes within fail category", ylim=c(0,6000), names.arg=c("fail 1","fail 2","fail 3","fail 4","fail 5")) 
dev.off() 
		#just for sake of visualising difference 
		#% plot gives sense of total number of genes (as % is calculate from all genes, not just failed ones) 
		#but this one includes numbers also - may be interesting 


###############################

#nonspec / compl / pseudo and failcats 

#################################
# nonspec v particular failcats #
#################################

#all are dichotomous (nonspec) + dichotomous (fail category) 
#therefore use Pearson's (same r as Phi coefficient for two binary variables) 
#save onto same table 

#fail 1 variables: 
	cor.nonspecVfail1 <- cor.test(Gnonspec, Gfail1) 
	r.nonspecVfail1 <- signif(cor.nonspecVfail1$estimate, digits = 4) 
	Pval.nonspecVfail1 <- signif(cor.nonspecVfail1$p.value, digits = 4) 
			Pval.nonspecVfail1 <- format(Pval.nonspecVfail1, scientific=FALSE) 

#fail 2 variables: 
	cor.nonspecVfail2 <- cor.test(Gnonspec, Gfail2) 
	r.nonspecVfail2 <- signif(cor.nonspecVfail2$estimate, digits = 4) 
	Pval.nonspecVfail2 <- signif(cor.nonspecVfail2$p.value, digits = 4) 
		Pval.nonspecVfail2 <- format(Pval.nonspecVfail2, scientific=FALSE) 

#fail 3 variables: 
	cor.nonspecVfail3 <- cor.test(Gnonspec, Gfail3) 
	r.nonspecVfail3 <- signif(cor.nonspecVfail3$estimate, digits = 4) 
	Pval.nonspecVfail3 <- signif(cor.nonspecVfail3$p.value, digits = 4) 
		Pval.nonspecVfail3 <- format(Pval.nonspecVfail3, scientific=FALSE)

#fail 4 variables: 
	cor.nonspecVfail4 <- cor.test(Gnonspec, Gfail4) 
	r.nonspecVfail4 <- signif(cor.nonspecVfail4$estimate, digits = 4) 
	Pval.nonspecVfail4 <- signif(cor.nonspecVfail4$p.value, digits = 4) 
		Pval.nonspecVfail4 <- format(Pval.nonspecVfail4, scientific=FALSE)

#fail 5 variables: 
	cor.nonspecVfail5 <- cor.test(Gnonspec, Gfail5) 
	r.nonspecVfail5 <- signif(cor.nonspecVfail5$estimate, digits = 4) 
	Pval.nonspecVfail5 <- signif(cor.nonspecVfail5$p.value, digits = 4) 
		Pval.nonspecVfail5 <- format(Pval.nonspecVfail5, scientific=FALSE)

#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Non-specific position v fail 1", r.nonspecVfail1, Pval.nonspecVfail1) 
		out2 <- print(line2, quote=FALSE) 
	line3 <- c("Non-specific position v fail 2", r.nonspecVfail2, Pval.nonspecVfail2) 
		out3 <- print(line3, quote=FALSE) 
	line4 <- c("Non-specific position v fail 3", r.nonspecVfail3, Pval.nonspecVfail3) 
		out4 <- print(line4, quote=FALSE) 
	line5 <- c("Non-specific position v fail 4", r.nonspecVfail4, Pval.nonspecVfail4) 
		out5 <- print(line5, quote=FALSE) 
	line6 <- c("Non-specific position v fail 5", r.nonspecVfail5, Pval.nonspecVfail5) 
		out6 <- print(line6, quote=FALSE) 
	
	write(out2, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out3, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out4, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out5, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out6, file = "GnonspecVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 


###############################
# compl v particular failcats #
###############################

#all are dichotomous (nonspec) + dichotomous (fail category) 
#therefore use Pearson's 
#save onto same table 

#fail 1 variables: 
	cor.complVfail1 <- cor.test(Gcompl, Gfail1) 
	r.complVfail1 <- signif(cor.complVfail1$estimate, digits = 4) 
	Pval.complVfail1 <- signif(cor.complVfail1$p.value, digits = 4) 
		Pval.complVfail1 <- format(Pval.complVfail1, scientific=FALSE)

#fail 2 variables: 
	cor.complVfail2 <- cor.test(Gcompl, Gfail2) 
	r.complVfail2 <- signif(cor.complVfail2$estimate, digits = 4) 
	Pval.complVfail2 <- signif(cor.complVfail2$p.value, digits = 4) 
		Pval.complVfail2 <- format(Pval.complVfail2, scientific=FALSE)

#fail 3 variables: 
	cor.complVfail3 <- cor.test(Gcompl, Gfail3) 
	r.complVfail3 <- signif(cor.complVfail3$estimate, digits = 4) 
	Pval.complVfail3 <- signif(cor.complVfail3$p.value, digits = 4) 
		Pval.complVfail3 <- format(Pval.complVfail3, scientific=FALSE)

#fail 4 variables: 
	cor.complVfail4 <- cor.test(Gcompl, Gfail4) 
	r.complVfail4 <- signif(cor.complVfail4$estimate, digits = 4) 
	Pval.complVfail4 <- signif(cor.complVfail4$p.value, digits = 4) 
		Pval.complVfail4 <- format(Pval.complVfail4, scientific=FALSE)

#fail 5 variables: 
	cor.complVfail5 <- cor.test(Gcompl, Gfail5) 
	r.complVfail5 <- signif(cor.complVfail5$estimate, digits = 4) 
	Pval.complVfail5 <- signif(cor.complVfail5$p.value, digits = 4) 
		Pval.complVfail5 <- format(Pval.complVfail5, scientific=FALSE)

#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Complement v fail 1", r.complVfail1, Pval.complVfail1) 
		out2 <- print(line2, quote=FALSE) 
	line3 <- c("Complement v fail 2", r.complVfail2, Pval.complVfail2) 
		out3 <- print(line3, quote=FALSE) 
	line4 <- c("Complement v fail 3", r.complVfail3, Pval.complVfail3) 
		out4 <- print(line4, quote=FALSE) 
	line5 <- c("Complement v fail 4", r.complVfail4, Pval.complVfail4) 
		out5 <- print(line5, quote=FALSE) 
	line6 <- c("Complement v fail 5", r.complVfail5, Pval.complVfail5) 
		out6 <- print(line6, quote=FALSE) 
	
	write(out2, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out3, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out4, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out5, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out6, file = "GcomplVGallcats.xls", ncolumns=3, append=TRUE, sep="\t")
	

################################
# pseudo v particular failcats #
################################

#all are dichotomous (nonspec) + dichotomous (fail category) 
#therefore use Pearson's 
#save onto same table 

#fail 1 variables: 
	cor.pseudoVfail1 <- cor.test(Gpseudo, Gfail1) 
	r.pseudoVfail1 <- signif(cor.pseudoVfail1$estimate, digits = 4) 
	Pval.pseudoVfail1 <- signif(cor.pseudoVfail1$p.value, digits = 4) 
		Pval.pseudoVfail1 <- format(Pval.pseudoVfail1, scientific=FALSE) 

#fail 2 variables: 
	cor.pseudoVfail2 <- cor.test(Gpseudo, Gfail2) 
	r.pseudoVfail2 <- signif(cor.pseudoVfail2$estimate, digits = 4) 
	Pval.pseudoVfail2 <- signif(cor.pseudoVfail2$p.value, digits = 4) 
			Pval.pseudoVfail2 <- format(Pval.pseudoVfail2, scientific=FALSE) 

#fail 3 variables: 
	cor.pseudoVfail3 <- cor.test(Gpseudo, Gfail3) 
	r.pseudoVfail3 <- signif(cor.pseudoVfail3$estimate, digits = 4) 
	Pval.pseudoVfail3 <- signif(cor.pseudoVfail3$p.value, digits = 4) 
			Pval.pseudoVfail3 <- format(Pval.pseudoVfail3, scientific=FALSE) 

#fail 4 variables: 
	cor.pseudoVfail4 <- cor.test(Gpseudo, Gfail4) 
	r.pseudoVfail4 <- signif(cor.pseudoVfail4$estimate, digits = 4) 
	Pval.pseudoVfail4 <- signif(cor.pseudoVfail4$p.value, digits = 4) 
			Pval.pseudoVfail4 <- format(Pval.pseudoVfail4, scientific=FALSE) 

#fail 5 variables: 
	cor.pseudoVfail5 <- cor.test(Gpseudo, Gfail5) 
	r.pseudoVfail5 <- signif(cor.pseudoVfail5$estimate, digits = 4) 
	Pval.pseudoVfail5 <- signif(cor.pseudoVfail5$p.value, digits = 4) 
			Pval.pseudoVfail5 <- format(Pval.pseudoVfail5, scientific=FALSE) 

#save values to table: 
	topline <- c("variables", "r", "p-value") 
	out1 <- print(topline, quote=FALSE) 
	write(out1, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 

	line2 <- c("Pseudogene v fail 1", r.pseudoVfail1, Pval.pseudoVfail1) 
		out2 <- print(line2, quote=FALSE) 
	line3 <- c("Pseudogene v fail 2", r.pseudoVfail2, Pval.pseudoVfail2) 
		out3 <- print(line3, quote=FALSE) 
	line4 <- c("Pseudogene v fail 3", r.pseudoVfail3, Pval.pseudoVfail3)
		out4 <- print(line4, quote=FALSE) 
	line5 <- c("Pseudogene v fail 4", r.pseudoVfail4, Pval.pseudoVfail4)
		out5 <- print(line5, quote=FALSE) 
	line6 <- c("Pseudogene v fail 5", r.pseudoVfail5, Pval.pseudoVfail5)
		out6 <- print(line6, quote=FALSE) 
	
	write(out2, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out3, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out4, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out5, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t") 
	write(out6, file = "GpseudoVGallcats.xls", ncolumns=3, append=TRUE, sep="\t")
