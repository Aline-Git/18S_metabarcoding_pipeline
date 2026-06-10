#!/bin/bash

'''this pipeline takes as input the file with the sequences corresponding to
18S rRNA gene generated with the script 18S_metabarcoding_pipeline_WP2_article.sh and clusters the 
sequences in decreasing order of quality. A consensus sequence is created for each cluster larger
than five and compared with PR2 for taxonomic assignment.'''


# usage : <path_to_script_directory>/18S_metabarcoding_clustering_WP2_article.sh <barcode_name> <path_to_workdir> <path_to_primer_file>
# built : 2025.01.16


########################################## INPUT #################################################

# path to working directory, needs to be the same as the one used with 18S_metabarcoding_pipeline_WP2_article.sh
WORKDIR=$1

# name of the clustering directory
CLUSTERING_DIR=$2 


######################################## PIPELINE  ###############################################


#________________________0. CREATE OUTPUT SUBDIRECTORIES_______________________#

mkdir $WORKDIR/$CLUSTERING_DIR
mkdir $WORKDIR/$CLUSTERING_DIR/centroids
mkdir $WORKDIR/$CLUSTERING_DIR/clusters
mkdir $WORKDIR/$CLUSTERING_DIR/polishing
mkdir $WORKDIR/$CLUSTERING_DIR/polishing/polished_centroids
mkdir $WORKDIR/$CLUSTERING_DIR/polishing/logs
mkdir $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5


#_____________________1. INPUT SEQUENCES QUALITY SORTING ______________________#

# sort all the sequences by quality, remove sequences with quality lower than 10
# Biopython v1.85
~/scripts/sort_reads_by_quality.py \
    -i $WORKDIR/clustering/input_sequences/grouped_18S_selection.fastq \
    -q 10

# remove or gzip the original file
gzip $WORKDIR/clustering/input_sequences/grouped_18S_selection.fastq


#_________________________________2. CLUSTERING ______________________________#

# trim and cluster the sequence with vsearch 2.29.4
conda activate vsearch_env
vsearch --fastx_filter \
  $WORKDIR/$CLUSTERING_DIR/input_sequences/grouped_18S_selection_qual_sorted.fastq \
  --fastq_minlen 1000 \
  --fastqout $WORKDIR/$CLUSTERING_DIR/grouped_qualsorted_filter_1000.fastq \
  --fastqout_discarded $WORKDIR/$CLUSTERING_DIR/discarded_grouped_qualsorted_filter_1000_discarded.fastq \
  --fastq_qmax 90
  

vsearch --cluster_smallmem \
  $WORKDIR/$CLUSTERING_DIR/grouped_qualsorted_filter_1000.fastq \
  --strand both \
  --usersort \
  --centroids $WORKDIR/$CLUSTERING_DIR/centroids/grouped_qualsorted_id98_filtered_1000_centroids.fasta \
  --uc $WORKDIR/$CLUSTERING_DIR/vsearch_clusters_id_grouped_98_filtered_qualsorted_min1000.uc \
  -id 0.98 \
  > $WORKDIR/$CLUSTERING_DIR/vsearch_clusters_id_grouped_98_filtered_qualsorted_min1000.log
  
 conda deactivate
 

# select clusters containing at least 5 sequences 
# Biopython v1.85
~/scripts/make_cluster_files.py \
  -c $WORKDIR/$CLUSTERING_DIR/vsearch_clusters_id_grouped_98_filtered_qualsorted_min1000.uc \
  -f $WORKDIR/$CLUSTERING_DIR/grouped_qualsorted_filter_1000.fastq \
  -o $WORKDIR/$CLUSTERING_DIR/clusters \
  -t fastq \
  -s 5
   
gzip $WORKDIR/$CLUSTERING_DIR/grouped_qualsorted_filter_1000.fastq


#__________________3. MAKING CONSENSUS SEQUENCES FOR CLUSTER _________________#

# polish the clusters with medaka
cd $WORKDIR/$CLUSTERING_DIR/clusters

conda activate medaka_env

for centroid in $(ls ASV*_centroid.fasta)
do 
#get the cluster_id

CLUSTER_ID="${centroid%*_centroid.fasta}"
# get the number of sequences
seq_count=$(grep -c "^>" "${CLUSTER_ID}.fasta")

echo ""
echo "Cluster ${CLUSTER_ID} has $seq_count sequences. Making consensus with medaka ..."
echo ""

# change the format to prevent compatibility with medaka
sed -i 's/:/_/g' ${CLUSTER_ID}.fasta
sed -i 's/:/_/g' $centroid

# polishing cluster with MEDAKA v2.0.1
medaka_consensus -t 2 \
   -i ${CLUSTER_ID}.fasta \
   -d $centroid \
   -o $WORKDIR/clustering/polishing/polished_centroids/${CLUSTER_ID} \
   > $WORKDIR/clustering/polishing/logs/medaka_${CLUSTER_ID}_polishing.log

 # put the polished sequence in an input file
cat $WORKDIR/clustering/polishing/polished_centroids/${CLUSTER_ID}/consensus.fasta >> $WORKDIR/clustering/polishing/centroids_consensus.fasta 

done

conda deactivate 


#_______________4. TAXONOMIC ASSIGNMENT OF CONSENSUS SEQUENCES ________________#


# taxonomic assignment with vsearch usearch_global minimum coverage of query sequence 90%, the two strand directions are compared with PR2
# VSEARCH v2.29.4
conda activate vsearch_env
vsearch \
  --usearch_global $WORKDIR/clustering_bioinfo3/polishing/centroids_consensus.fasta \
  -db ~/reference/database/PR2_db/pr2_version_5.1.0_SSU_mothur.fasta \
  --id 0.1 \
  -query_cov 0.9 \
  --strand both \
  --blast6out $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5/usearch_output_query_cov_90_both_strand.txt \
  --output_no_hits
 
conda deactivate


# select the best match in case of both direction producing a match
~/scripts/select_best_strand_in_usearch_output.py \
  -i $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5/usearch_output_query_cov_90_both_strand.txt \
  -o $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5/usearch_output_query_cov_90_best_strand.txt



# make a taxonomic table from the output clustering file (.uc) and vsearch usearch_global output of the polished centroids.
# Biopython v1.85
~/scripts/make_grouped_taxonomic_table_from_usearch_output.py \
  -c $WORKDIR/clustering_bioinfo3/vsearch_clusters_id_grouped_98_filtered_qualsorted_min1000.uc \
  -o $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5 \
  -s 5 \
  -t ~/reference/database/PR2_db/pr2_version_5.1.0_SSU_mothur.tax \
  -b $WORKDIR/clustering_bioinfo3/usearch_global_cov90_entroids_18S_clusters_size5/usearch_output_query_cov_90_best_strand.txt

echo "**** END OF THE PIPELINE ****"












