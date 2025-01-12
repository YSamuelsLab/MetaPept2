#!/bin/bash

### MetaPept Version 2.0

### module section (UNCOMMENT IT ONLY IF YOU USE THE SCRIPT ON A CLUSTER!!!!) ###

#module load jre/8.121
#module load Python/3.10.8-GCCcore-12.2.0
#module load Biopython/1.81-foss-2022b
#module load R/4.3.1-foss-2022b


### main section (uncomment lines you need to use in the pipeline) ###

echo "MetaPept is running..."

### WARNING: the MetaPept pipeline uses `sample_description.csv` file to describe input datasets. DO NOT USE any symbols in `Sample_Name` column except letters and numbers.
### If you prepare `sample_description.csv` file in M$ Windows, please, remove Windows specific new-line symbols from the file.
### you can use these commands:
### mv sample_description.csv sample_description.win64.csv
### perl -pe 's/\r//g' sample_description.win64.csv > sample_description.new.csv

### Notes for MSFragger: there is no a dedicated branch for FragPipe in the MetaPept.
### Instead of that the MetaPept has a set of scripts to transform FragPipe output (combined_peptide.tsv and sample_name/psm.tsv) to MaxQuant-like output (msms.txt and peptides.txt)

### step 1.1 Optional (MSFragger): 
# If the run has only one replica, MSFragger does not create combined_peptide.tsv file. Use the script to create combined_peptide.tsv file based on psm.tsv file
#python3 scripts/MSF_combined_peptide_maker.py -s sample_description.csv -i path_to_FragPipe_output_dirSampleName/psm.tsv -o path_to_FragPipe_output_dir/combined_peptide.tsv

### step 1.2 (MSFragger): before the FragPipe output transformation, remove decoy sequences from the FragPipe output 
### (reference.fasta - the initial reference FASTA file you used with FragPipe; reference.fasta.fas - the FASTA file which was created by FragPipe after adding decoy sequences):
#python3 scripts/ContaminationSearch.py -fa path_to_reference.fasta -fas path_to_reference.fasta.fas -i path_to_FragPipe_output_dir/combined_peptide.tsv -o path_to_FragPipe_output_dir/combined_peptide.no_contaminant.tsv

### step 1.3 (MSFragger): the transformation of the MSFragger output to the MaxQuant output (msms.txt and peptides.txt)
#python3 scripts/MSF_combiner.py -d path_to_FragPipe_output_dir/ -c path_to_FragPipe_output_dir/combined_peptide.no_contaminant.tsv

### step 1.4 (MaxQuant/MSFragger): MaxQuant IMP (FASTA file names can be specified with '-c' and '-n' command options)
### to run the script without binding prediction use `--dummy` key; the maximum peptide length can be set with `--max-len` key
### for example:
### python3 scripts/IMP.py --dummy --max_len 100 -b tools/netMHCpan/Linux_x86_64 -d data/ORFs -c CDS.fasta -n nuORFdb.fasta -s input/sample_description.csv -a input/hla_alleles.csv -i input/maxquant/Canonical -o output/maxquant/Canonical
#python3 scripts/IMP.py -b tools/netMHCpan/Linux_x86_64 -d data/ORFs -c CDS.fasta -n nuORFdb.fasta -s input/sample_description.csv -a input/hla_alleles.csv -i input/maxquant/Aeffect -o output/maxquant/Aeffect
#python3 scripts/IMP.py -b tools/netMHCpan/Linux_x86_64 -d data/ORFs -c CDS.fasta -n nuORFdb.fasta -s input/sample_description.csv -a input/hla_alleles.csv -i input/maxquant/Canonical -o output/maxquant/Canonical
#python3 scripts/IMP.py -b tools/netMHCpan/Linux_x86_64 -d data/ORFs -c CDS.fasta -n nuORFdb.fasta -s input/sample_description.csv -a input/hla_alleles.csv -i input/maxquant/Peffect -o output/maxquant/Peffect

### step 2: MaxQuant/MSFragger scan validation
#python3 scripts/scan_validation.py -a output/maxquant/Aeffect/IMP_filtered.csv -b output/maxquant/Canonical/IMP_filtered.csv -o output/maxquant/Aeffect/
#python3 scripts/scan_validation.py -a output/maxquant/Peffect/IMP_filtered.csv -b output/maxquant/Canonical/IMP_filtered.csv -o output/maxquant/Peffect/

### step 3: combine MaxQuant/MSFragger validated files
#python3 scripts/mq_combiner.py -s input/sample_description.csv -i output/maxquant/Aeffect/IMP_scan_validation.csv output/maxquant/Peffect/IMP_scan_validation.csv -n A_effect P_effect -o output/maxquant/mq_combined.csv

### step 4: preparing a batch file for Peptide-PRISM
### default aliases for PRISM categories:
### frameshift=CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic
### prio1=CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic  (the same as 'frameshift')
### prio2=CDS,Extra,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic
### prio3=CDS,UTR5,OffFrame,UTR3,ncRNA,Extra,Intronic,Intergenic
###
### to specify non-default aliases use the notation: prio{4-9}=comma_delimited_list_of_categories
### for example -cat prio4=CDS,UTR5,OffFrame prio5=CDS,Extra,OffFrame prio6=CDS,UTR5,OffFrame,ncRNA,Extra
#python3 scripts/prism_batch_file_maker.py -g hs.90 -hla input/hla_alleles.csv -extra data/proxyPhe/A375_PA_all.fasta -r run_prism.sh -cat frameshift prio2 prio3 -i input/peaks -o output/prism

### step 5: run Peptide-PRISM
#echo "running PRISM..."
#bash run_prism.sh
#rm run_prism.sh
#echo "PRISM: done"

### step 6.1 (optional): modification of Peptide-PRISM output without binding prediction to make it compartable with the downstream pipeline
#python3 scripts/binding_prediction.py -i output/prism

### step 6.2: combine Peptide-PRISM output including decoy peptides
#python3 scripts/prism_combiner.py -s input/sample_description.csv -cat frameshift prio2 prio3 -i output/prism -r 2.0 -o output/prism/prism_combined.binders.csv
#python3 scripts/prism_combiner.py -D -s input/sample_description.csv -i output/prism -o output/prism/prism_combined.decoy.csv

### step 7: scan integration
#python3 scripts/pipeline_integrator.py -s input/sample_description.csv -mq output/maxquant/mq_combined.csv -prism output/prism/prism_combined.binders.csv -o output/scan_integration.csv

### step 8: data filtering and combining with PRISM unique peptides (Intensities are taken from MQ output for 'combined' hits and from PRISM output for 'PRISM unique' hits)
### applying default filters:
### 1. 'combined' hit has Best_ALC >= 80; 
### 2. 'combined' hit has 70 <= Best_ALC < 80 and Delta >= 10 and Best_coverage >= 80; 
### 3. 'PRISM_unique' hit has Best_Q <= 0.1 and Best_ALC >= 80 and HLA_Rank < 2.0
###
### to specify non-default thresholds use the command options: -alc_suff, -alc_comb, -cov_comb, -delta_comb, -alc_uni, -q_uni, -rank_uni (see README.md)
#python3 scripts/integration_filter.py -com output/combined_scan_integration.csv -uni output/prism_unique_scan_integration.csv -o output/scan_integration.filtered.csv

### step 9: searching the peptides in a FASTA file taking into account I2L modifications 
### 'strict' mode prevents I2L variant searching, 'peptide' and 'protein' modes are for FASTA files with peptides and proteins respectively
### 'peptide' mode (only for FASTA files with peptides) implements exact match between query and target peptides
### 'protein' mode implements full-text search across sequences
#python3 scripts/search_in_fasta.py -i output/scan_integration.filtered.csv -p Sequence -n HLA_atlas_WEB -f data/HLA_atlas/WEB/HLA_atlas.2020_12.HLA_I_and_HLA_IplusII.fa -t peptide -o output/scan_integration.filtered.HLA_atlas_WEB.csv
#python3 scripts/search_in_fasta.py -i output/scan_integration.filtered.HLA_atlas_WEB.csv -p Sequence -n HLA_atlas_PRISM -f data/HLA_atlas/PRISM/HLA_atlas.Q_less_0.1__ALC_more_70__Rank_less_2.0.fasta -t peptide -o output/scan_integration.filtered.HLA_atlas_WEB_PRISM.csv
#python3 scripts/search_in_fasta.py -strict -i output/scan_integration.filtered.HLA_atlas_WEB_PRISM.csv -p Sequence -n IEDB -f data/IEDB/epitope_table_export_1668681696.linear_peptides.fa -t peptide -o output/scan_integration.filtered.HLA_atlas_WEB_PRISM.IEDB.csv

### step 10 (optional): cryptic peptides filtering
#python3 scripts/search_in_fasta.py -i output/scan_integration.filtered.HLA_atlas_WEB_PRISM.IEDB.csv -p Sequence -n Extra -f data/proxyPhe/A375_PA_all.fasta -t protein -o output/scan_integration.filtered.HLA_atlas_WEB_PRISM.IEDB.Extra.csv
#python3 scripts/integration_filter.exp.py -i output/scan_integration.filtered.HLA_atlas_WEB_PRISM.IEDB.Extra.csv -o output/scan_integration.experimental.csv

echo -e "\n...all tasks are completed"
