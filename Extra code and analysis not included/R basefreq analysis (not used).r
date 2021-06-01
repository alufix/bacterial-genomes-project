#extra analysis (base/GC frequency analysis, out of interest) 
#decided not to include in analysis or investigate any further (as need to be selective with tests) for interest of time 


######################
### SIZE (NT) V GC ###
######################

#sizeNT v GCcontent (out of interest, just because I can) 
#do not include in main report because not relevant to my investigation 

#variables: 
	slope.sizeVgc <- signif(coef(lm(GCcontent~sizeNT))[2], digits = 4) 
	cor.sizeVgc <- cor.test(sizeNT, GCcontent) 
	r.sizeVgc <- signif(cor.sizeVgc$estimate, digits = 4) 
	Pval.sizeVgc <- signif(cor.sizeVgc$p.value, digits = 4) 

#plot graph: 
pdf("sizeNTvGC.pdf")
plot(sizeNT, GCcontent, xlab = "genome size (NT)", ylab = "GC content", col = "red", cex = 0.5, pch = 18, main = "genome size (NT) v GC content")
	abline(lm(GCcontent~sizeNT), col = "black") 
	mtext(paste0("gradient = ", slope.sizeVgc), side = 1, adj = 0.95, line = -3.3, cex = 0.8) 
	mtext(paste0("r = ", r.sizeVgc), side = 1, adj = 0.95, line = -2.3, cex = 0.8)
	mtext(paste0("p-value = ", Pval.sizeVgc), side = 1, adj = 0.95, line = -1.3, cex = 0.8)
dev.off() 



###########################
### SIZE (NT) BASE FREQ ###
###########################

#sizeNT v baseFreqA/T/C/G (out of interest, just because I can) 
#do not include in main report because not relevant to my investigation 

#plot graph (with subplots): 
pdf("baseFreq_matrix.pdf", width = 21/2.54, height = 27.7/2.54) 

par(mfrow=c(2,2), oma = c(4,1,1,1) + 0.1, mar = c(3,2,3,2) + 0.1)
	#above is for customising 
	#par = parameters (of the plot) 
	#mfrow = c(2,2) - means want a plot that is 2 x 2 (so 4 altogether) 
	#other bits do control of margins etc. - can play around 
	
	#v baseFreqA 
		slope.sizeVfreqa <- signif(coef(lm(baseFreqA~sizeNT))[2], digits = 4) 
		cor.sizeVfreqa <- cor.test(sizeNT, baseFreqA) 
		r.sizeVfreqa <- signif(cor.sizeVfreqa$estimate, digits = 4) 
		Pval.sizeVfreqa <- signif(cor.sizeVfreqa$p.value, digits = 4) 
		
	plot(sizeNT, baseFreqA, xlab = "genome size (NT)", ylab = "base frequency A", col = "red", cex = 0.5, cex.lab = 0.7, cex.axis = 0.7, pch = 18, main = "genome size (NT) v frequency of base A", cex.main = 0.9)
	abline(lm(baseFreqA~sizeNT), col = "black") 
	mtext(paste0("gradient = ", slope.sizeVfreqa), side = 3, adj = 0.95, line = -1.3, cex = 0.7) 
	mtext(paste0("r = ", r.sizeVfreqa), side = 3, adj = 0.95, line = -2.3, cex = 0.7)
	mtext(paste0("p-value = ", Pval.sizeVfreqa), side = 3, adj = 0.95, line = -3.3, cex = 0.7)
	
	#v baseFreqT 
		slope.sizeVfreqt <- signif(coef(lm(baseFreqT~sizeNT))[2], digits = 4) 
		cor.sizeVfreqt <- cor.test(sizeNT, baseFreqT) 
		r.sizeVfreqt <- signif(cor.sizeVfreqt$estimate, digits = 4) 
		Pval.sizeVfreqt <- signif(cor.sizeVfreqt$p.value, digits = 4) 
		
	plot(sizeNT, baseFreqT, xlab = "genome size (NT)", ylab = "base frequency T", col = "green", cex = 0.5, cex.lab = 0.7, cex.axis = 0.7, pch = 18, main = "genome size (NT) v frequency of base T", cex.main = 0.9)
	abline(lm(baseFreqT~sizeNT), col = "black") 
	mtext(paste0("gradient = ", slope.sizeVfreqt), side = 3, adj = 0.95, line = -1.3, cex = 0.7) 
	mtext(paste0("r = ", r.sizeVfreqt), side = 3, adj = 0.95, line = -2.3, cex = 0.7)
	mtext(paste0("p-value = ", Pval.sizeVfreqt), side = 3, adj = 0.95, line = -3.3, cex = 0.7)
	
	#v baseFreqC 
		slope.sizeVfreqc <- signif(coef(lm(baseFreqC~sizeNT))[2], digits = 4) 
		cor.sizeVfreqc <- cor.test(sizeNT, baseFreqC) 
		r.sizeVfreqc <- signif(cor.sizeVfreqc$estimate, digits = 4) 
		Pval.sizeVfreqc <- signif(cor.sizeVfreqc$p.value, digits = 4) 
	
	plot(sizeNT, baseFreqC, xlab = "genome size (NT)", ylab = "base frequency C", col = "orange", cex = 0.5, cex.lab = 0.7, cex.axis = 0.7, pch = 18, main = "genome size (NT) v frequency of base C", cex.main = 0.9)
	abline(lm(baseFreqC~sizeNT), col = "black") 
	mtext(paste0("gradient = ", slope.sizeVfreqc), side = 1, adj = 0.95, line = -3.3, cex = 0.7) 
	mtext(paste0("r = ", r.sizeVfreqc), side = 1, adj = 0.95, line = -2.3, cex = 0.7)
	mtext(paste0("p-value = ", Pval.sizeVfreqc), side = 1, adj = 0.95, line = -1.3, cex = 0.7)
	
	#v baseFreqG
		slope.sizeVfreqg <- signif(coef(lm(baseFreqG~sizeNT))[2], digits = 4) 
		cor.sizeVfreqg <- cor.test(sizeNT, baseFreqG) 
		r.sizeVfreqg <- signif(cor.sizeVfreqg$estimate, digits = 4) 
		Pval.sizeVfreqg <- signif(cor.sizeVfreqg$p.value, digits = 4) 
	
	plot(sizeNT, baseFreqG, xlab = "genome size (NT)", ylab = "base frequency G", col = "blue", cex = 0.5, cex.lab = 0.7, cex.axis = 0.7, pch = 18, main = "genome size (NT) v frequency of base G", cex.main = 0.9)
	abline(lm(baseFreqG~sizeNT), col = "black") 
	mtext(paste0("gradient = ", slope.sizeVfreqg), side = 1, adj = 0.95, line = -3.3, cex = 0.7) 
	mtext(paste0("r = ", r.sizeVfreqg), side = 1, adj = 0.95, line = -2.3, cex = 0.7)
	mtext(paste0("p-value = ", Pval.sizeVfreqg), side = 1, adj = 0.95, line = -1.3, cex = 0.7)

	par(pty="s")
		#supposed to make the subplots square? 
	
dev.off() 

