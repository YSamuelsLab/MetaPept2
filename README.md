# MetaPept V2.0

## The integrated pipeline for immunopeptidomics.

### Dependencies:

- Java Runtime Environment (JRE) >= 8.0
- Python >= 3.10.6 has to be installed with the packages: pandas, numpy, difflib, sqlite3, itertools, Levenshtein, shutil, pathlib, Bio (Biopython), rpy2
- R >= 4.1.2 has to be installed with the packages: protViz, ggplot2, cowplot, plyr, uniReg, reshape2, snow
- Peptide-PRISM has to be installed in 'tools/Peptide-PRISM' directory
- netMHCpan has to be installed in 'tools/netMHCpan' directory (gawk needs to be installed)

### Installation on local linux computer:

1) install R 
2) install R packages:
    >$ R  
    >install.packages(c("protViz", "ggplot2", "cowplot", "plyr", "uniReg", "reshape2", "snow"))  
    >q()
3) install Pyhon packages from the requirements:
    >$ pip3 install -r requirements.python3.txt
4) install Java Runtime Environment (JRE) >= 8.0
5) copy Peptide-PRISM into 'tools/Peptide-PRISM' directory
6) set genome version you need in 'initialization.PRISM.sh' script and run it
7) copy netMHCpan into 'tools/netMHCpan' directory
8) set 'NMHOME' variable in 'tools/netMHCpan/netMHCpan' file (see netMHCpan.readme)
9) make sure that gawk is on the system (if not, install it)

### Installation on LSF cluster:

1) load the next modules under your account:
   - jre/8.121
   - Python/3.10.8-GCCcore-12.2.0
   - Biopython/1.81-foss-2022b
   - R/4.3.1-foss-2022b

2) install R packages:
    >$ R  
    >install.packages(c("protViz", "ggplot2", "cowplot", "plyr", "uniReg", "reshape2", "snow"))  
    >q()  
3) and AFTER that install Pyhon packages from the requirements:
    >$ pip3 install --user -r requirements.python3.txt

4) copy Peptide-PRISM into 'tools/Peptide-PRISM' directory
5) set genome version you need in 'initialization.PRISM.sh' file and run 'initialization.PRISM.on_cluster.sh' script
6) copy netMHCpan into 'tools/netMHCpan' directory 
7) set the 'NMHOME' variable in 'tools/netMHCpan/netMHCpan' file (see netMHCpan.readme)
8) make sure that gawk is on the system (if not install it)

### Usage

To run the pipeline:
1) copy MaxQuant output files 'msms.txt' and 'peptides.txt' into 'input/maxquant' directory (FragPipe output must be transformed to MaxQuant output - see MetaPept.sh file)
2) copy PEAKS output files 'all de novo candidates.csv' and 'de novo peptides.csv' in 'input/peaks' directory (they can be in sub-directories)
3) create a tab-delimited sample description file with the columns:
   - Source_File (*.raw file name)
   - Sample_Name
   - Sample_Replica (1, 2, 3, etc.)
   - Sample_Type (ko - knock out, wt - wild type, etc.)
   - Group (samples with the same group name will be processed together and saved in the individual output file)
   - Experiment (experiment name from the MaxQuant experimental design file)

    The script 'make_sample_description.sh' helps to create a pre-built sample description file.
    Names of RAW files will be added automatically by the script.
    All you need is to add additional information to the file and move it to the input directory.
    The file 'EXAMPLE.sample_description.csv' in the main directory can be used as an example.  
    DO NOT USE any symbols in `Sample_Name` column except letters and numbers.
    
4) Uncomment pipeline steps you need in the main pipeline script 'MetaPept.sh'
5) Run 'MetaPept.sh' script (if you use the pipeline on LSF cluster, run 'MetaPept.on_cluster.sh' script instead of 'MetaPept.sh') 

#### The default aliases for categories which are used for running Peptide-PRISM:
- frameshift (or prio1): CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic   
- prio2: CDS,Extra,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic  
- prio3: CDS,UTR5,OffFrame,UTR3,ncRNA,Extra,Intronic,Intergenic  

At the step 4 of the pipeline you can specify you own categories using the notation: prio{4-9}=category_1,category_2,...,category_n  
For example: prio4=Extra,CDS,Extra,UTR5,OffFrame

#### The default filter for the pipeline report file:
- 'combined' hit has Best_ALC >= 80
- 'combined' hit has 70 <= Best_ALC < 80 and Delta >= 10 and Best_coverage >= 80
- 'PRISM_unique' hit has Best_Q <= 0.1 and Best_ALC >= 80 and HLA_Rank < 2.0

At the step 8 of the pipeline you can separately specify threshold values for combined and PRISM unique hits:   
-alc_suff (sufficient Best ALC threshold for combined hits - there are no other filters for this threshold)  
-alc_comb (Best ALC threshold for combined hits with the additional filters for Coverage and Delta)  
-cov_comb (Best coverage threshold - it is used together with alc_combined)  
-delta_comb (Delta threshold - it is used together with alc_combined)  
-alc_uni (Best ALC threshold for PRISM unique hits)  
-q_uni (Best Q threshold for PRISM unique hits)  
-rank_uni (HLA rank for PRISM unique hits - HLA rank value for MaxQuant hits are filtered several steps before)  

### Notes

For some peptides PRISM sets Intensity = 0

### TO DO

Add column 'location_count' in the PRISM combined file
Replace MaxQuant with MSFragger 
