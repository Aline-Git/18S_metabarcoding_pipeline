#!/bin/bash

'''this pipeline trims primer sequences, filter short reads and select the 18S rRNA gene
regions, if any, for individual taxonomic assignment with PR2. The selected sequences 
are concatenated in a single output file for further clustering with the script 
18S_metabarcoding_clustering_WP2_article.'''


# usage : <path_to_script_directory>/18S_metabarcoding_pipeline_WP2_article.sh <barcode_name> <path_to_workdir> <path_to_primer_file> <min_final_seq_len>
# built : 2025.01.16
# INPUT : the fastq_pass output of a ONT sequencing run with multiple barcode folders

########################################## INPUT #################################################

barcode_list=$1 # name of the current barcode folder
WORKDIR=$2 # path to the directory where all the results will be written
PRIMER_FILE=$3 # a file with the primer sequences in fasta format
min_final_seq_len=$4 # input parameter, should be an integer between 100 and 4000


######################################## PIPELINE  ###############################################


#________________________0. CREATE OUTPUT SUBDIRECTORIES_______________________#

mkdir $WORKDIR/input_folder_fastq
mkdir $WORKDIR/input_folder_fasta
mkdir $WORKDIR/logs
mkdir $WORKDIR/blast_primers
mkdir $WORKDIR/read_length
mkdir $WORKDIR/adapter_cut
mkdir $WORKDIR/barrnap
mkdir $WORKDIR/usearch_18S_selection_pr2v51


#________________________1. RAW INPUT FILES PREPROCESSING______________________#

for barcode in $barcode_list
do

  RAW_FASTQ=$WORKDIR/input_folder_fastq/${barcode}/${barcode}_pass_all.fastq
  RAW_FASTA=$WORKDIR/input_folder_fasta/${barcode}/${barcode}_pass_all.fasta


  # create a folder for merged fastq files 
  mkdir $WORKDIR/input_folder_fastq/${barcode}

  # create a merged fastq file
  cd $WORKDIR/fastq_pass/${barcode}

  for filename in $(ls *.fastq.gz)
    do 
    cat $fileanme >> ${RAW_FASTQ}.gz
  
    done

  gzip -d ${RAW_FASTQ}

  # create a fasta file 
  mkdir $WORKDIR/input_folder_fasta/${barcode}

  sed -n '1~4s/^@/>/p;2~4p' $RAW_FASTQ > $RAW_FASTA

  # count the total number of raw sequences and store it in a log file
  seq_count=$(grep -c "^>" $RAW_FASTA)

  echo "passed_sequences, ${seq_count}" >> $WORKDIR/logs/${barcode}_recap_count_sequences.csv


  # check the presence of primers (I like to do this step as a quality control but it's not 100% mandatory (can be included at the end) 
  '''checking the presence of the primer sequences in the raw sequences can be
  useful to detect potential anomalies(presence of primer tandem repeat, multiple
  amplicons in one sequence, etc.)'''
 
  # NCBI BLAST v2.16.0+
  conda activate blast_env
  mkdir $WORKDIR/blast_primers/${barcode}

  # search for matches with a minimum percentage id of 90% 
  # only the matches with at least 85% length are kept in the sumary table
  ~/scripts/blast_get_primer_table/main.py \
    -i $RAW_FASTA \
    -f $PRIMER_FILE \
    -p 90 \
    -l 85 \
    -o $WORKDIR/blast_primers/${barcode}	

  conda deactivate


  # get sequence lengths
  mkdir $WORKDIR/read_length/${barcode}

  # produce a tsv table with sequenceID{tab}sequence_length
  # uses Biopython v1.85
  ~/scripts/utils/get_seqlen.py \
    -f $RAW_FASTA \
    -l 0 \
    -t fastq \
    -o $WORKDIR/read_length/${barcode}/read_length.csv


  #_______________________________1. PRIMER REMOVAL_____________________________#

  ''' Adapter trimming and filtering sequences shorter than 500bp, followed by two
  rounds supplementary adapter trimming to remove remaining adapters due to primer 
  tandem repeats or sequences containing multiple amplicon'''

  # This part will need generalisation because primer sequences will vary.
  # I think about reading the primer from a file with a python script, then reverse-complementing
  # the reverse primer(s) and calling cutadapt from there 
  # Do you think it is a good idea ? Or do you have a proposition for something simpler ?

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

  # convert the trimmed file into fasta 
  sed -n '1~4s/^@/>/p;2~4p' $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq > $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta


  # count the number of trimmed sequences and store this information in the log file
  seq_count2=$(grep -c "^>" $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta)

  echo "adapter_removal, ${seq_count2}" >> $WORKDIR/logs/${barcode}_recap_count_sequences.csv


  #_________________________2. selection of 18S region___________________________#

  '''we use barrnap to extract the rRNA sequences from covariance models.
  First we search for all the rRNAs, then we select the sequences that have the
  best match with the 18S model, the 18S region of the selected sequences is 
  trimmed and put in a file '''


  mkdir $WORKDIR/barrnap/${barcode}

  # BARRNAP v0.9
  conda activate barrnap_env 

  # detect eukaryotic 18S and 23S rRNA using barrnap v0.9
  barrnap -k euk --lencutoff 0.6 --threads 12 \
    $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta \
    > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_euk.csv \
    2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_euk.log

  # detect eukaryotic 18S and 23S rRNA using barrnap v0.9
  barrnap -k bac --lencutoff 0.6 --threads 12 \
    $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta \
    > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_bact.csv \
    2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_bact.log

  # detect eukaryotic 18S and 23S rRNA using barrnap v0.9
  barrnap -k arc --lencutoff 0.6 --threads 12 \
    $$WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta \
    > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_arc.csv \
    2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_arc.log
  
  # detect eukaryotic 18S and 23S rRNA using barrnap v0.9  
  barrnap -k mito --lencutoff 0.6 --threads 12 \
    $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta \
    > $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_mito.csv \
    2> $WORKDIR/barrnap/${barcode}/${barcode}_barrnap_mito.log
  
  conda deactivate


  mkdir $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID

  # Biopython v1.85
  # select and trim 18S rRNA sequences (can be adapted to fit bacterial sequences)
  ~/scripts/extract_18S_sequences_with_barcodeID/main_biopython.py \
    -f $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq \
    -b $WORKDIR/barrnap/${barcode}/ \
    -o $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fastq \
    -s ${barcode} \
    -t fastq

  # count the number of trimmed sequences and store this information in the log file
  seq_count3=$(grep -c "^>" $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta)

  echo "adapter_removal, ${seq_count2}" >> $WORKDIR/logs/${barcode}_recap_count_sequences.csv
  
  
  
  #_____________________3. individual taxonomic assignment______________________#


  mkdir $WORKDIR/usearch_18S_selection_pr2v51/${barcode}


  ~/scripts/quality_and_length_filtering.py \
    -i $WORKDIR/barrnap/${barcode}/input_16S_file_with_barcodeID/${barcode}_16S_selection.fastq \
    -q 10 \
    -l $min_final_seq_len
  
  # we use vsearch --usearch_global to compare individual sequences with the database, both orientation are compared, only the best one will be kept

  # VSEARCH v2.29.4
  conda activate vsearch_env
  vsearch \
    --usearch_global $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection_q10_l${min_final_seq_len}.fastq \
    -db ~/reference/database/PR2_db/pr2_version_5.1.0_SSU_mothur.fasta \
    --id 0.1 \
    -query_cov 0.9 \
    --blast6out $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90.txt \
    --output_no_hits
  conda deactivate

done 


''' the files $WORKDIR/usearch_18S_selection_pr2v51/${barcode}/${barcode}_usearch_output_query_cov_90.txt 
can be used to build taxonomic tables (in R for example).'''


#___3. concatenate the individual quality and length filtered fastq files with____#

for barcode in $barcode_list
do

  # add the selected sequence to the file grouping the sequences from all the barcode for the clustering
  SELECTED_FASTQ=$WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection.fastq 
  cat $SELECTED_FASTQ >> $WORKDIR/$CLUSTERING_DIR/input_sequences/grouped_18S_selection.fastq

done

# gzip some 
for barcode in $barcode_list
do
gzip $RAW_FASTQ
gzip $RAW_FASTA

gzip $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fastq
gzip $WORKDIR/adapter_cut/${barcode}/${barcode}_pass_all_cutadapt.fasta 

gzip $WORKDIR/barrnap/${barcode}/input_16S_file_with_barcodeID/${barcode}_16S_selection.fastq
gzip $WORKDIR/barrnap/${barcode}/input_18S_file_with_barcodeID/${barcode}_18S_selection_q10_l${min_final_seq_len}.fastq

done

