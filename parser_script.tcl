

proc Complement s {

	global revs 
	 

	regsub -all {a} $s {1} s
	regsub -all {c} $s {2} s 
	regsub -all {t} $s {3} s
	regsub -all {g} $s {4} s

	regsub -all {1} $s {t} s
	regsub -all {2} $s {g} s
	regsub -all {3} $s {a} s
	regsub -all {4} $s {c} s

	set revs [string reverse $s]
	
return 
}



#make folders for failed sequences 
	if {![file exists failedsequences]} {file mkdir failedsequences}
		#this says if a file called failedsequences does not exist (!), make one 

	if {![file exists failcats]} {file mkdir failcats}
		#'failcats' folder to save 'fail category' files of each genome 	

#make file to keep log of fails: 
	set out_error [open errors.csv w] 
	puts $out_error "accn,locus_tag,fail" 


set flist [glob -nocomplain -directory filtered_genomes *.embl] 

#flist is filelist 
#'-directory genomes' means go into the directory 'genomes' and add to list of files 
#all input files end in .embl 


set num_genomes [llength $flist]
	#llength is 'list length' - returns the number of things in the list 

	
#make spreadsheet to keep record of genome details - 'overview' (of genomes) 
	set out_overview [open overview.csv w]
	puts $out_overview "accn,organism,sizeNT,numGene,numFail,pcFail,fail1,fail2,fail3,fail4,fail5,fail,GCcontent,baseFreqA,baseFreqT,baseFreqC,baseFreqG,institution,seqYear,reassessed,reasYear,seqMethod,QC"
		#makes header of excel sheet for listing overview of details for every genome 

#make spreadsheet to keep record of gene details - 'overviewGenes' 
	set out_overviewGenes [open overviewGenes.csv w]
	puts $out_overviewGenes "accn,organism,loctag,exact,name,compl,pseudo,fail,numFailcats,fail1,fail2,fail3,fail4,fail5,institution,seqYear" 

#make spreadsheet to keep record of failed gene details - 'overviewFails' 
	set out_overviewFails [open overviewFails.csv w]
	puts $out_overviewFails "accn,organism,loctag,exact,name,compl,pseudo,fail,numFailcats,fail1,fail2,fail3,fail4,fail5,institution,seqYear" 


#start counting genomes 
set cnt 1 


#handling a list 
#so to make a loop for list items: 'foreach' 

foreach genome $flist {

#have output telling us which number genome we are on (out of total), along with accession number of genome 
	regexp {genomes/(.*?)\.embl} $genome all accn 
	puts "Analysing $accn: $cnt of $num_genomes" 
	
	#(therefore if script breaks, can see which genome broke the script)

#make file for each genome 
	if {![file exists $accn]} {file mkdir $accn}


#make 'failcats' file (fail categories) 
	#if {![file exists $accn-failcats]} {
	set out_failcats [open $accn-failcats.xls w] 
	#} 

#pulling out sequences 
#first open file and set contents to 'read' 
	set in [open $genome r]
	set contents [read $in]
	close $in 


#create a list of 'failed genes' for each genome 
	set genomeFails [list ] 


#now pull out specific bits of code from the genome 

set firstXXOS [string first XX\nOS $contents]
	#find first XX newline 'OS' in the contents 
set OS_pos [string first OS $contents $firstXXOS]
	#search the string for the first 'OS' in contents, starting from the position of the firstXXOS 
set NL_pos [string first \n $contents $OS_pos] 

	set speciesline [string range $contents [expr $OS_pos + 2] $NL_pos]
	set speciesline [string trim $speciesline]
	set parts [split $speciesline]
	set genus [lindex $parts 0] 
	
	regsub -all { } $parts {_} parts 
	#(replacing space in organism name with underscore so .csv will read correctly)


#date of sequencing: 

set firstXXDT [string first XX\nDT $contents]
	#find first XX newline 'DT' in the contents 

	set nextXX [string first \nXX $contents $firstXXDT] 
	#find the next newline 'XX' in the contents, after the first position of XX\nDT which we just set 
	#this is essentially checkpointing the end of the date section 
	set dateInfo [string range $contents $firstXXDT $nextXX] 
	#setting this whole date section as $dateInfo (will use this string later to pull out reassessment by regexp) 
	
	set seqYear [string range $contents [expr $firstXXDT + 15] [expr $firstXXDT + 18]] 
	#$firstXXDT pos is set previously, and year starts at position 15 (counting from 0 position after firstXXDT), and ends at position 18 
	#therefore says: cut out the range, in $contents, starting at the position of the $firstXXDT + 15 and ending at the position + 18 
	#therefore have the whole year being saved 


#reassessment: 
#(note: all genomes ended up being positive for reassessment, so not used in analysis) 

	if {[regexp {Last updated} $dateInfo] ==1} {
	set reassessed 1 
		set reasPOS [string first {Last updated} $dateInfo] 
		set reasYear [string range $dateInfo [expr $reasPOS - 16] [expr $reasPOS - 13]] 
	
	} else { 
	set reassessed 0 
		set reasYear {NaN}  
	}


#institution:  

set firstRL [string first "\nRL   Submitted" $contents]
	#find first newline RL "submitted" in the contents 
	
	set nextRL [string first ".\nRL" $contents]
	#hopefully this will pull out the NEXT newline break RL AFTER the first one (so cut out the 'Submitted on' bit and just go to the address) 
	#depends on how string first works - run script and check 
	
	set endofRL [string first "\nXX" $contents $nextRL] 
	#this is finding the end of the RL institute bit: 
	#says set the checkpoint $endofRL at the first instance of \nXX in the $contents starting from the set position $firstRL (where the section starts) 
	#(before setting string range as the bit you want to pull out, need to define the start and end 'checkpoints', then position the text within) 
	
	set rough_inst [string range $contents $nextRL $endofRL] 
	#set rough_inst (rough institute) as the range of the $contents string, from the starting position '$nextRL' up until the end position '$endofRL'
	#this will come out as a 3-line set of text, including commas, which will separate the output into columns (as .csv) 
	#therefore need to clean up $rough_inst by removing linebreaks (replace with whitespace?), and removing any commas 
	#ideally, would want to clean up to the point where can just isolate the country and lab name 
	# however this is more complex - possible, but for this level of analysis, not worth it 
	# can also probably assume that the format of reporting the lab address is identical for each institute, so will be able to group the same institutes 
	
	regsub -all {,} $rough_inst {} rough_inst 
		#delete all commas (so that output is not separated at position of commas) 
	regsub -all {\.} $rough_inst {} rough_inst 
		#delete all fullstops (note: . has to have an escape char before)
	regsub -all {\n} $rough_inst { } rough_inst 
		#sub all line breaks with a space (so that all on one line) 
	regsub -all {RL   } $rough_inst {} rough_inst 
		#sub all RL indent bits with space 
	regsub -all { } $rough_inst {_} clean_inst 
		#sub all spaces with underscore (so have proper .csv formatting) 
	#don't worry TOO much about formatting here - but note that assuming that every institute's address is in exactly same format 
	#if all same format, then this code will do the same thing to each genome, then will be able to group identical ones 
		#(note: format of institutions ended up being messy and inconsistent, e.g. same address having different formats 
		#tried to reformat and group the same institutions together - see other files in USB - 
		#but ended up not including institution in analysis) 


#now get sequence: 
		
set XX_pos [string first "XX\nSQ" $contents] 
	#position of first instance of XX linebreak SQ 

set other_pos [string first "other;" $contents $XX_pos] 
	#now narrowing down to the first instance of "other;" after XX SQ 

set endline_pos [string first \n $contents $other_pos]
	#this is targeting it to the end of the line (\n) (the line that ends in 'other;')
	#i.e. targeting it to the line where the sequence begins 

set raw_seq [string range $contents $endline_pos end]
#pulling out the raw sequence 
#range is within the contents, starting at the 'endline_pos' position you just defined and ending at end of document 

regsub -all {[^A-Za-z]} $raw_seq {} clean_seq
#cleaning up the raw sequence 
#replaces all instances of anything other than letters (i.e. whitespace, tabs, numbers) with nothing (empty brackets)
#i.e. deleting everything other than letters 
#calling this new thing 'clean_seq' (clean sequence)




#count number of bases (each base and total): 
	set countA [regsub -all {a} $clean_seq "a" clean_seq] 
	set countT [regsub -all {t} $clean_seq "t" clean_seq]
	set countC [regsub -all {c} $clean_seq "c" clean_seq] 
	set countG [regsub -all {g} $clean_seq "g" clean_seq] 
	set sumNT [expr $countA + $countT + $countC + $countG]



#calculate base frequency (of each base): 
	set baseFreqA [expr $countA. / $sumNT]
	set baseFreqT [expr $countT. / $sumNT]
	set baseFreqC [expr $countC. / $sumNT]
	set baseFreqG [expr $countG. / $sumNT]


#calculate GC content: 
	set GCcontent [expr $baseFreqC + $baseFreqG] 

		#(note: base and GC content was considered but not included in final analysis) 
	

	#note that the above process that counts number of NT does not include any ambiguous (non-ACTG) bases 
	#therefore for the genes/genomes containing ambiguous bases, this is not entirely representative of the number of NT 
	#could change this to count any character 
	#however some 'messy' genomes may contain long stretches of NNN, and sometimes repeated stretches, therefore would not be 100% accurate in any case 
	#in addition, some sequences may also include other ambiguous base variations - e.g. ?A for representing 'either A or no nucleotide' - thus counting as '2 bases' would be inaccurate 
	#so in this project, number of genome in NT is only including ACTG bases 


#now want to extract the annotation (above the sequence) 
	#can use a split command to break up the important bits of the annotation 


set annot [string range $contents 0 $XX_pos]
#getting annotation: from 0 position to the end position where the sequence starts 
#line of code literally means 'setting annotation as the string range from the 'contents', 
#starts at the 0 position (very first character in the file) and ending at the XX_pos checkpoint we defined above 


unset raw_seq 
unset contents 
#unsetting variables no longer needed 



#now breaking up genes: first insert £, then split around £ (to turn from string into list) 

regsub -all {FT   CDS             } $annot "£FT   CDS             " annot 
#literally copy from the genome in BBedit from FT (indent bit) to CDS to the sequence 

set gene_list [split $annot £]
#create a new variable 'gene_list' 
	#define gene_list as the splitting of the 'annot' sequence at every £ symbol 
	#split around the £ signs you just put in (note: will remove them) 
	#puts these as elements in a list 

set gene_list [lrange $gene_list 1 end]
#setting list for starting at position 1 until the end of the section 



set gene_cnt 1 

	foreach gn $gene_list {
	#introducing 'gn' for this loop - so every time you mention gene_list you are substituting it into gn in this loop script 
	

	#get location:
	
		regexp {[^0-9]([0-9]+?)\.\.([0-9]+?)[^0-9]} $gn all strt nd 
	
		#pulling out the location of sequence (position) - e.g. 110616..111377 
		#(strt = 110616; nd = 111377) 

	
	
	#now pulling out locus tag: 

		if {[regexp {/locus_tag="(.*?)"} $gn]==1} {
			regexp {/locus_tag="(.*?)"} $gn all loctag 
			} else {set loctag $gene_cnt} 
		
		#pulls out the locus tag and sets as 'loctag'; if does not have loctag set then name loctag the number of the gene 
		#now clean up loctag (by removing commas and spaces) so that it doesn't mess up csv output file: 
			regsub -all {,} $loctag {} loctag
			regsub -all { } $loctag {_} loctag 


	#do the same for the gene name: 

		if {[regexp {/gene="(.*?)"} $gn]==1} {
			regexp {/gene="(.*?)"} $gn all genename  
			} else {set genename $gene_cnt}

		#now clean up genename (by removing commas and spaces) so that it doesn't mess up csv output file 
			regsub -all {,} $genename {} genename
			regsub -all { } $genename {_} genename  
		
	
		
	#now account for 'complement': 

		regexp {(.*?)\n} $gn all firstline 
		
		if {[regexp {complement} $firstline]} {
			set comp 1
		} else {
			set comp 0
		}
		
		#if it says 'complement' in the first line, then set comp as 1, otherwise set comp as 0 
		
		
	#now account for partial CDS (non-specific position, i.e. whether position contains > or <) 
	#(referred to in this case as 'posQuality' - i.e. referencing quality of positioning) 

		if {[regexp {>|<} $firstline] ==1} {
			set posQuality "non-specific" 
		} else { 
			set posQuality "exact" 
		} 
	

	#now account for 'pseudogene': 
		#look in whole of annotation ($gn) for pseudo or (pseudogene) 
		#(note: post-analysis, realised this may not be accounting for all identifiers of pseudogene in annotation - need to reconsider for future use) 

		if {[regexp {pseudo/|(pseudogene)} $gn] ==1} {
			set pseudo 1  
		} else { 
			set pseudo 0  
		} 

	
	#this covers the whole annotation section for each $gn	

	#now get sequence: 
		
		set seq [string range $clean_seq [expr $strt -1] [expr $nd -1]]
		set seq [string tolower $seq] 
		
		if {$comp ==1} { 
		Complement $seq 
		set seq $revs 
		} 
	
		#setting variable 'seq' (sequence of each $gn) to the $clean_seq variable we extracted above 
		# ranging from '$strt' (the first number in the gene location, i.e. '1426'..1498) -1 to account for the computer language 
		# (genome file has first codon at 1. but when searching for '1' location in Tcl, want it to output 0 position base) 
		# (therefore we do this calculation $strt -1 to put the number back 1 so that the Tcl search returns the right base) 
	
		#setting sequence to lowercase for sake of consistency 
	
		#then saying if comp is present (if $comp ==1), then run the Complement procedure for that sequence 
		# and the new 'seq' variable for this complement sequence is set to the value of $revs (at the end of the Complement proc)

	

#now quality control: telling us whether there is anything odd about the sequence 
	

	set fail 0 
	
		set fail1 0
		set fail2 0
		set fail3 0
		set fail4 0
		set fail5 0
	
	#create a list of 'fail categories' for each gene 
	set failcats [list ]
		#failcats == list of fail categories  
	
	
	while {$fail ==0} {
	
	#setting fail as 0 so that you can process each gene,
	# and then when you encounter a dodgy one, you can change the 'fail' variable to a value 1/2/3/4/5 
	
	#5 categories of fails 
		#1: is it divisible by 3 or not? 
		#2: are there letters other than atgc? 
		#3: does it end in a non-stop codon? 
		#4: does it contain an internal stop codon? 
		#5: does the first codon not end in TG? 
	#if any of these fails are present, we want to assign it with a fail value and move it somewhere else 
	#because don't want to be using bad gene files 
	#so we are building up a library for the type of fail 
	
	#FAIL 1: sequence should be a multiple of 3 - if not, want it to be a fail: 
		
		set seql [string length $seq] 
		if {[expr $seql % 3] >0} {set fail 1
		puts "$loctag: fail 1"
		
		set fail1 1
		
			#if the list failcats does not already contain a 1, append the list with " 1"
			if {![regexp "1" $failcats]==1} {
			append failcats " 1"
	
			}
		
		}
		
		#defining sequence length 'seql' 
		#second line says: if the length divided by 3 gives a remainder 
		# (i.e. if the length is not divisible by 3) (divisible by 3 should give remainder of 0)
		# then change the value of 'fail' variable to '1'  
		
		#last line: append failcats - add 1 to the list of fail categories 
		#says if the failcats list does not already contain a 1, then add one (as don't want duplicates) 

		
	#FAIL 2: if there are letters other than actg: 
	
		if {[regexp {[^a|c|t|g]} $seq]} {set fail 2
		puts "$loctag: fail 2"
		
		set fail2 1 
		
			if {![regexp "2" $failcats]==1} {
			append failcats " 2"
			}
		
		}

		#says if a pattern is found for anything other than a, c, t, or g; set fail 2 
	
	
	#FAIL 3: if it ends in a non-stop codon (something other than TGA, TAA, TAG): 
	
		set last_codon [string range $seq end-2 end]
		if {[regexp {[^tga|taa|tag]} $last_codon]==1} {set fail 3
		puts "$loctag: fail 3"
		
		set fail3 1
		
			if {![regexp "3" $failcats]==1} {
			append failcats " 3"
			}

		} 
		
		#starting by setting last_codon' so that we can isolate the last 3 bases and test that 
		#then says 'if the last codon is something other than the 3 stop codons (if this is true), set fail 3' 


	#FAIL 4: if it contains an internal stop codon: 
	
		for {set i 0} {$i < [expr $seql - 5]} {incr i +3} {
			set cdn [string range $seq $i [expr $i +2]]
			
			if {[regexp {tga|taa|tag} $cdn] ==1} {set fail 4
			puts "$loctag: fail 4"
			
			set fail4 1
		
				if {![regexp "4" $failcats]==1} {
				append failcats " 4"
				}
		
			
			} 
			}
		
		#setting i and incrementing so that script is able to read through codons (by multiple of 3 bases) 
		#then says 'if codon contains one of the stop codons, set fail as 4' 
		

	
	#FAIL 5: 
	
		set first_codon [string range $seq 1 2]
		if {[regexp {tg} $first_codon]==0} {set fail 5
		puts "$loctag: fail 5"
		
		set fail5 1
		
			if {![regexp "5" $failcats]==1} {
			append failcats " 5"
			}
		
		
		
		} 
	
		#setting first codon (but only looking at positions 1 and 2, not 0 - i.e. last 2 bases of first codon)
		#then says 'if 1st and 2nd position does not contain tg, set fail to 5' 

	
	
	if {$fail==0} {set fail -1} 
	
	#important for exiting quality control loop 
	#says 'anything that is still fail = 0 by the end of this loop, change fail to -1 and then exit' 
	


	} 
	
	
	
#now remaining code is creating the gene files that we have just processed in the above loop 

	puts $out_failcats "$loctag $failcats"
	#want it to record the gene locus_tag and list of fail categories independent of whether any fails were observed 
 

if {$fail >0}  {
	
	#i.e. if it is a fail (1,2,3,4,5)
	

#if fail is greater than 0, add the gene identifier (locus tag) to the $genomeFail list 
	append genomeFails " $loctag"

	set failornot 1 
		#if gene has failed, set $failornot 1 
		
	set numFailcats [llength $failcats] 
		#count the number of fails in the failcats list, and set it to $numFailcats 

	puts $out_error "$accn,$loctag,$failcats" 
	
		#says '(for each failed gene,) write down the genome accession number, locus tag, and type of fail' 
		#(in the csv document)
		# $failcats used (rather than $fail) so that it reports all of the fail categories encountered, not just most recent one 


	set out [open $accn-$loctag.tfa w] 
	puts $out "locus_tag: >$loctag\nposition: $strt-$nd ($posQuality)\nname: $genename\ncomplement?: $comp\npseudogene?: $pseudo\nfail category: $failcats\nsequence:\n$seq" 
	close $out 
	
		#says '(for each failed gene,) create a file named with the genome accession number and locus tag 
		#and in the document/file, write the relevant info (locus tag, etc.) 
		


	file rename -force $accn-$loctag.tfa failedsequences/$accn-$loctag.tfa 
		#this forces the failed sequence file to the 'failedsequences' folder which we created at the beginning 

		
#and save info to 'overview of fails' spreadsheet 
	puts $out_overviewFails "$accn,$parts,$loctag,$posQuality,$genename,$comp,$pseudo,$failornot,$numFailcats,$fail1,$fail2,$fail3,$fail4,$fail5,$clean_inst,$seqYear" 

	
	
} else { 

	#i.e. if it has not failed (if fail is not greater than 0) 

	set failornot 0 
		#if gene has not failed, set $failornot 0 
	
	set numFailcats 0 

	set out [open $loctag.tfa w] 
	puts $out "locus_tag:>$loctag; position:$strt-$nd; compl=$comp; name:$genename\n$seq"
	close $out 	
	
		#says 'make a file named with the locus tag' 
		#(note: ended up not using these files in analysis) 

	file rename -force $loctag.tfa $accn/$loctag.tfa 
	
		#moves the un-failed gene files to a folder named with its accession number 

	}
	

#save info to 'ovewview of genes' spreadsheet (regardless of fail) 	
#(so spreadsheet will contain both failed genes and non-failed genes) 
	puts $out_overviewGenes "$accn,$parts,$loctag,$posQuality,$genename,$comp,$pseudo,$failornot,$numFailcats,$fail1,$fail2,$fail3,$fail4,$fail5,$clean_inst,$seqYear" 

	#(note that this file is so massive that it does not open in Excel) 

	
incr gene_cnt 	
	}
	#now we increment the count of the genes processed by 1 (at the end of each run, in anticipation of the next one)


#NEW: as the overview output is below, i.e. after the loop for each gene within the genome, gene_cnt is currently +1 extra 
#in order to report the correct number of genes, need to decrease this count by 1 
#(will use 'true_gene_cnt' variable) 
	set true_gene_cnt [expr $gene_cnt -1] 


#want to close 'failcats.xls' file as you reach the end of the genome loop 
#note: a few lines down, we are going to open this file again to 'read' it 

	close $out_failcats 



#before we go to the next genome, we want to put a line out on the overview document for each genome (to save info for each genome) 
#this uses the variables set in the above loop for each genome: the accession number, number of genes, and length of genomeFails etc. 	
	set numFails [llength $genomeFails] 
		#llength counts the number of failed genes in the genomeFails list 
	set percent [expr (($numFails./$true_gene_cnt)*100)]
		#expresses the percentage of failed genes per total number of genes in each genome 
		#need the . after the term to make sure calculation is performed as a float, not an integer (otherwise will be 0) 
		#don't bother putting % in output as the header already has one 
		

#counting number of fail categories encountered for the single genome (i.e. how many genes encountered each type of fail) 
	set in_failcats [open $accn-failcats.xls r]
	set catcontents [read $in_failcats] 
	close $in_failcats 
	#above i am saving the contents of the failcats excel file so that i can read the contents below 
	#catcontents = contents of 'fail categories' file 
	#(note that i initially opened the failcats file with write permissions, then closed this, now opening with read again - perhaps not efficient, but still works) 


	set num_fail1 [regsub -all { 1} $catcontents " 1" catcontents] 
	set num_fail2 [regsub -all { 2} $catcontents " 2" catcontents]
	set num_fail3 [regsub -all { 3} $catcontents " 3" catcontents] 
	set num_fail4 [regsub -all { 4} $catcontents " 4" catcontents] 
	set num_fail5 [regsub -all { 5} $catcontents " 5" catcontents] 

		#need to count the number of fail categories in the genome-specific failcats document 
		#here we are saving, to the variable num_failX, the number of times each fail category is included in the genome doc 
		#(the regsub command returns the number of substitutions it has performed 
		#(so if we are not substituting with anything different, will just count and return the total to be saved to the num_failX variable)  



#now we can move the file to the 'failcats' folder 
	file rename -force $accn-failcats.xls failcats/$accn-failcats.xls 
	
	if {$percent==0} {set ifFail 0} else {set ifFail 1} 
	
	puts $out_overview "$accn,$parts,$sumNT,$true_gene_cnt,$numFails,$percent,$num_fail1,$num_fail2,$num_fail3,$num_fail4,$num_fail5,$ifFail,$GCcontent,$baseFreqA,$baseFreqT,$baseFreqC,$baseFreqG,$clean_inst,$seqYear,$reassessed,$reasYear,seqMethod,QC"
#"accn,organism,size(nt),numGene,numFail,pcFail,fail1,fail2,fail3,fail4,fail5,ifFail,GCcontent,baseFreqA,baseFreqT,baseFreqC,baseFreqG,institution,seqYear,reassessed,reasYear,seqMethod,QC"



incr cnt
	#we do the same for the count of genomes processed (e.g. if have 7 genomes, after 1st one, $cnt is increased to 2)


}



close $out_error 
close $out_overview 
close $out_overviewGenes
close $out_overviewFails 
