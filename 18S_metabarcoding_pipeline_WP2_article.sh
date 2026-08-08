#!/bin/bash

'''this pipeline trims primer sequences, filter short reads and select the 18S rRNA gene
regions, if any, for individual taxonomic assignment with PR2. The selected sequences 
are concatenated in a single output file for further clustering with the script 
18S_metabarcoding_clustering_WP2_article.'''


# usage : <path_to_script_directory>/18S_metabarcoding_pipeline_WP2_article.sh <path_to_workdir> <path_to_primer_file> <path_to_database_dir>
# built : 2025.01.16

########################################## INPUT #################################################

barcode=$1
WORKDIR=$2
PRIMER_FILE=$3
DATABASE_FOLDER=$4

######################################## PIPELINE  ###############################################


#________________________0. CREATE OUTPUT SUBDIRECTORIES_______________________#

mkdir $WORKDIR/input_folder_fastq
mkdir $WORKDIR/input_folder_fasta
mkdir $WORKDIR/logs
mkdir $WORKDIR/blast_primers
mkdir $WORKDIR/adapter_cut
mkdir $WORKDIR/barrnap
mkdir $WORKDIR/usearch_18S_selection_pr2v51


#________________________1. RAW INPUT FILES PREPROCESSING______________________#

for barcode_dir in $WORKDIR/fastq_pass/barcode*/
do 

  # Extract barcode name from path
  barcode=$(basename "$barcode_dir")
    
  echo "Processing $barcode..."
  
  RAW_FASTQ=$WORKDIR/input_folder_fastq/${barcode}/${barcode}_pass_all.fastq
  RAW_FASTA=$WORKDIR/input_folder_fasta/${barcode}/${barcode}_pass_all.fasta


  # create a folder for merged fastq files 
  mkdir $WORKDIR/input_folder_fastq/${barcode}

  # TODO : cat the zipped files
  # unzip fastq files
  cd $WORKDIR/fastq_pass/${barcode}

  for filename in $(ls *.fastq.gz)
    do 
    gzip -d $filename
    done

  # group fastq cat und unzip raw fastq files
  for filename in $(ls *.fastq)
    do 
      filebase=${filename%'.fastq'}
      cat ${filebase}.fastq >> $RAW_FASTQ
      gzip $filename
    done

  # create a fasta file 
  mkdir $WORKDIR/input_folder_fasta/${barcode}

  sed -n '1~4s/^@/>/p;2~4p' $RAW_FASTQ > $RAW_FASTA

  # count the total number of raw sequences and store it in a log file
  seq_count=$(grep -c "^>" $RAW_FASTA)

  echo "passed_sequences, ${seq_count}" >> $WORKDIR/logs/${barcode}_recap_count_sequences.csv


  # check the presence of primers 
  '''checking the presence of the primer sequences in the raw sequences can be
  useful to detect potential anomalies(presence of primer tandem repeat, multiple
  amplicons in one sequence, etc.)'''
 
  # NCBI BLAST v2.16.0+
  conda activate blast_env
  mkdir $WORKDIR/blast_primers/${barcode}

  # search for matches with a minimum percentage id of 90% 
  # only the matches with at least 85% length are kept in the sumary table
  18S_metabarcoding_pipeline/blast_get_primer_table/main.py \
    -i $RAW_FASTA \
    -f $PRIMER_FILE \
    -p 90 \
    -l 85 \
    -o $WORKDIR/blast_primers/${barcode}	

  conda deactivate


#_______________________________2. PRIMER REMOVAL_____________________________#

''' Adapter trimming and filtering sequences shorter than 500bp, followed by two
rounds supplementary adapter trimming to remove remaining adapters due to primer 
tandem repeats or sequences containing multiple amplicon'''

mkdir $WORKDIR/adapter_cut/${barcode}
mkdir $WORKDIR/adapter_cut/${barcode}/input_cut_file

# CUTADAPT v5.0
conda activate cutadapt_env
#TODO : read the primers from the file

# first round of primer removal
cutadapt -e 0.2 \
  --overlap 16 \
  --no-indels \
  --rc \
  -a GGCAAGTCTGGTGCCAG...AAGGTAGCCAAATGCCTCGTC \
  -a GGCAAGTCTGGTGCCAG...AYTWGTGACGYGCATGAATGG \
  -m 500 \
  -o $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq \
  $RAW_FASTQ \
  > $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.stats 

# second round of primer removal if there are primer tandem repeats
cutadapt -e 0.2 \
  --overlap 16 \
  --no-indels \
  --rc \
  -g "GGCAAGTCTGGTGCCAG;rightmost" \
  -a "AAGGTAGCCAAATGCCTCGTC" \
  -a "AYTWGTGACGYGCATGAATGG" \
  -m 500 \
  -o $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt_round2.fastq \
  $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq \
  > $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt_round2.stats 

# last round removing the reads that still have F primers at this step 
cutadapt -e 0.2 \
  --overlap 16 \
  --no-indels \
  --rc \
  -g "GGCAAGTCTGGTGCCAG;rightmost" \
  --untrimmed-output $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fastq \
  -m 500 \
  -o $WORKDIR/adapter_cut/${barcode}/${barcode}_many_F_primers_sequences.fastq \
  $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt_round2.fastq \
  > $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt_final.stats 
 
conda deactivate


# convert the trimmed file into fasta 
sed -n '1~4s/^@/>/p;2~4p' $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fastq > $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta


# count the number of trimmed sequences and store this information in the log file
seq_count2=$(grep -c "^>" $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta)

echo "adapter_removal, ${seq_count2}" >> $WORKDIR/logs/${barcode}_recap_count_sequences.csv


#_________________________3. selection of 18S region___________________________#

'''we use barrnap to extract the rRNA sequences from covariance models.
First we search for all the rRNAs, then we select the sequences that have the
best match with the 18S model, the 18S region of the selected sequences is 
trimmed and put in a file '''

mkdir $WORKDIR/barrnap/${barcode}

# BARRNAP v0.9
conda activate barrnap_env 

# detect eukaryotic 18S and 23S rRNA using barrnap v0.9
barrnap -k euk --lencutoff 0.6 --threads 12 \
  $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta \
  > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_euk.csv \
  2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_euk.log

# detect eukaryotic 18S and 23S rRNA using barrnap v0.9
barrnap -k bac --lencutoff 0.6 --threads 12 \
  $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta \
  > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_bact.csv \
  2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_bact.log

# detect eukaryotic 18S and 23S rRNA using barrnap v0.9
barrnap -k arc --lencutoff 0.6 --threads 12 \
  $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta \
  > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_arc.csv \
  2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_arc.log
  
# detect eukaryotic 18S and 23S rRNA using barrnap v0.9  
barrnap -k mito --lencutoff 0.6 --threads 12 \
  $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta \
  > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_mito.csv \
  2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_mito.log
  
conda deactivate


mkdir $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID

# Biopython v1.85
# select and trim 18S rRNA sequences
18S_metabarcoding_pipeline/extract_18S_sequences_with_barcodeID/main_biopython.py \
    -f $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fastq \
    -b $WORKDIR/barrnap/${barcode}/ \
    -o $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fastq \
    -s ${barcode} \
    -t fastq

# add the selected sequence to the file grouping the sequences from all the barcode for the clustering
SELECTED_FASTQ=$WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fastq 
cat $SELECTED_FASTQ >> $WORKDIR/$CLUSTERING_DIR/input_sequences/grouped_18S_selection.fastq


#_____________________4. individual taxonomic assignment______________________#


mkdir $WORKDIR/usearch_18S_selection_pr2v51/${barcode}

# we select the sequences that are at least 1000bp 
./18S_metabarcoding_pipeline/select_seq_by_size.py \
  -f $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fasta \
  -t fasta -l 1000 \
  -o $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection_min1000.fasta
  

# we use vsearch --usearch_global to compare individual sequences with the database, both orientation are compared, only the best one will be kept

# VSEARCH v2.29.4
conda activate vsearch_env
vsearch \
  --usearch_global $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection_min1000.fasta \
  -db $DATABASE_FOLDER/pr2_version_5.1.0_SSU_mothur.fasta \
  --id 0.1 \
  -query_cov 0.9 \
  --strand both \
  --blast6out $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90_strand_both.txt \
  --output_no_hits
conda deactivate


# usearch_global compare the strands in two directions. Now we select the best (ie the correct) one.
./18S_metabarcoding_pipeline/select_best_strand_in_usearch_output.py \
  -i $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90_strand_both.txt \
  -o $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90_best_strand.txt
 
''' the files $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90_best_strand.txt 
can be used to build taxonomic tables.'''

#_____________________5. zipping files______________________#

gzip $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fasta

gzip $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fastq
gzip $WORKDIR/adapter_cut/${barcode}/input_cut_file/${barcode}_pass_all_cutadapt_final.fasta

gzip $WORKDIR/adapter_cut/${barcode}/${barcode}_many_F_primers_sequences.fastq
gzip $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt_round2.fastq
gzip $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq

done



