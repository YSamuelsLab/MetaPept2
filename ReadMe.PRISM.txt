
1. Index a genome (has to be done once for the human genome)

Run from your commandline:

java -jar IndexGenome.jar -organism homo_sapiens -version 90 -n h.ens90 -p -nostar -nobowtie -nokallisto

java -jar IndexGenome.jar -a Ecoli_MRE600.gtf -s Ecoli_MRE600.fasta -nostar -nobowtie -nokallisto

java -jar IndexGenome.jar -a lisa.gb.gtf -s lisa.gb.fasta -nostar -nobowtie -nokallisto



This will download the human genome and its annotation from ensembl, and create the indices needed by FindPep. You can do this anywhere you want (preferrably some folder where you would like to store the genome etc.)

Estimated time: 10min.

2. Search a Peaks-result file (.csv) in the genome

java -jar FindPep.jar -in XYZ.csv -seq h.ens90 -progress -D -plot

This will do the following steps:

1. Search the genome and transcriptome for the peptides in XYZ.csv.
2. Produce an output file name XZY.csv.annotated.csv. All peptides sequences without any hit are discarded. If the sequence is found in several locations, several lines are reported. The difference of the ALC of each hit to the max ALC is computed (Delta first). If Delta First>15 (can be configured by the command line paramter -df), the line is discarded. The output contains the additional columns: 
	- PSM rank: Original rank for this PSM
	- Location count: Number of locations, this peptide is found.
	- Delta first: The ALC difference of this PSM to the (potentially unmatched) top PSM of this spectrum
	- Decoy: Is it a match to a decoy (D) or target (T)
	- Location: The genomic mapping location
	- Annotation: The annotation of the mapping location (CDS/UTR3/UTR5/OffFrame/ncRNA/Intronic/Intergenic)
	- Sequence: The matched sequence (i.e. I/L resolved)
3. Go through this file, only consider lines with CDS, and determine the best spectrum for each matched sequence (or the first one, if there are several with the same ALC). Go through the file once more (again only considering CDS lines, and the best spectra per matched sequence as determined before), and compute the difference of the best matched sequence to the second best matched sequence (Delta next). If it is less than 16 (or parameter -dn), discard the whole spectrum.
4. Go through the result of 3. and collect statistics for FDR: For each peptide length and ALC value: compute the number of target PSMs, decoy PSMs and ambiguous PSMs (where the same sequence occurs in target and decoy DB). Based on that, compute for each peptide length and ALC score two q values: 
	- Q1: [Number of decoys with >ALC] / [Number of targets with >ALC]
	- Q2: [Number of decoys with >ALC] / ([Number of targets with >ALC] + [Number of ambiguous hits with >ALC]) * f
f=([Total number of identified targets] + [Total number of identified ambiguous hits]) / [Total number of identified targets]
Thus, to get all peptides with FDR>1%, use all lines where Q1>0.01 (or Q2>0.01).

Steps 3 and 4 are also done for lines with (CDS/OffFrame/ncRNA/UTR3/UTR5), i.e. everything in the transcriptome, and once more with all lines.

Finally, the produced tables are used to generate some pdfs (R has to be installed with the packages ggplot2, cowplot, plyr, uniReg, reshape2, snow ), try to run "Rscript" from command line; otherwise, omit the -plot parameter).

Instead of specifying a genome name (as produced by IndexGenome), you can directly specify a fasta file (containing DNA or protein). CDS is reported for all hits!

Estimated time: 15 min

