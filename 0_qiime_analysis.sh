#! /bin/sh

# Diversity
"""

conda activate qiime2-2022.2

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ../manifest\ file.tsv \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data paired-end-demux.qza \
--o-visualization demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 210 \
  --p-trunc-len-r 130 \
  --p-n-threads 4 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza &

qiime metadata tabulate \
 --m-input-file stats-dada2.qza \
 --o-visualization stats-dada2.qzv

 #renaming some intermediate files
cp rep-seqs-dada2.qza rep-seqs.qza
cp table-dada2.qza table.qza

# create the feature table!!!!
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime tools export \
  --input-path table.qza \
  --output-path exported-asv-table

biom convert\
  --to-tsv\
  -i exported-asv-table/feature-table.biom\
  -o exported-asv-table/feature-table.txt\
  --table-type "OTU table"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# calculate diversity metrics
# note rarefaction level used here
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 5000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

# alpha rarefaction
# Lowest sample is SRR7690081 (13900 features - is not an outlier!)
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 5000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

# beta div significance  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Host_disease \
  --o-visualization core-metrics-results/unweighted-unifrac-host-disease-significance.qzv \
  --p-pairwise

# taxonomy analysis
qiime feature-classifier classify-sklearn \
  --i-classifier /media/data/class2022/silva-138-99-nb-weighted-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# taxonomy visualization
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots.qzv

qiime feature-table heatmap \
  --i-table table.qza \
  --m-sample-metadata-file metadata.tsv \
  --m-sample-metadata-column "Host_disease"\
  --p-method "ward"\
  --o-visualization asv-heatmap.qzv

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table genus-table.qza

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table family-table.qza

qiime feature-table heatmap \
  --i-table genus-table.qza \
  --m-sample-metadata-file metadata.tsv \
  --m-sample-metadata-column "Host_disease" \
  --p-method "ward" \
  --p-color-scheme "nipy_spectral" \
  --o-visualization genus-heatmap.qzv

qiime feature-table heatmap \
  --i-table family-table.qza \
  --m-sample-metadata-file metadata.tsv \
  --m-sample-metadata-column "Host_disease"\
  --p-method "ward"\
  --p-color-scheme "nipy_spectral" \
  --o-visualization family-heatmap.qzv

# Identify differentially present bacteria between male BP and HC
# No differentially present bacteria were found. DESeq2 was then used
# to account for covariates effect.
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[Host_disease]='male_HC' OR [Host_disease]='male_BP'" \
  --o-filtered-table male-table.qza

# transform to a compositional table (porportions, add pseudo count to all cells)
qiime composition add-pseudocount \
  --i-table male-table.qza \
  --o-composition-table comp-male-table.qza

qiime composition ancom \
  --i-table comp-male-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Host_disease \
  --o-visualization ancom-male.qzv

# Identify differentially present bacteria between female BP and HC
# No differentially present bacteria were found. DESeq2 was then used
# to account for covariates effect.
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --p-where "[Host_disease]='female_HC' OR [Host_disease]='female_BP'" \
  --o-filtered-table female-table.qza

# transform to a compositional table (porportions, add pseudo count to all cells)
qiime composition add-pseudocount \
  --i-table female-table.qza \
  --o-composition-table comp-female-table.qza

qiime composition ancom \
  --i-table comp-female-table.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Host_disease \
  --o-visualization ancom-female.qzv
