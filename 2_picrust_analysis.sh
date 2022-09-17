#! /bin/sh

conda activate qiime2-2021.11


#Perform Picrust analysis for males
nohup \
qiime picrust2 full-pipeline \
--i-table male-table.qza \
--i-seq rep-seqs.qza \
--output-dir q2-picrust2_output_male \
--p-threads 4 &

qiime feature-table summarize \
  --i-table q2-picrust2_output_male/ko_metagenome.qza \
  --o-visualization q2-picrust2_output_male/ko_metagenome.qzv

qiime feature-table summarize \
  --i-table q2-picrust2_output_male/ec_metagenome.qza \
  --o-visualization q2-picrust2_output_male/ec_metagenome.qzv

qiime feature-table summarize \
  --i-table q2-picrust2_output_male/pathway_abundance.qza \
  --o-visualization q2-picrust2_output_male/pathway_abundance.qzv

qiime feature-table filter-features \
  --i-table q2-picrust2_output_male/ko_metagenome.qza \
  --p-min-frequency 50 \
  --p-min-samples 4 \
  --o-filtered-table q2-picrust2_output_male/ko_metagenome_filtered.qza

#Perform Picrust analysis for females

nohup \
qiime picrust2 full-pipeline \
--i-table female-table.qza \
--i-seq rep-seqs.qza \
--output-dir q2-picrust2_output_female \
--p-threads 5 &

qiime feature-table summarize \
  --i-table q2-picrust2_output_female/ko_metagenome.qza \
  --o-visualization q2-picrust2_output_female/ko_metagenome.qzv

qiime feature-table summarize \
  --i-table q2-picrust2_output_female/ec_metagenome.qza \
  --o-visualization q2-picrust2_output_female/ec_metagenome.qzv

qiime feature-table summarize \
  --i-table q2-picrust2_output_female/pathway_abundance.qza \
  --o-visualization q2-picrust2_output_female/pathway_abundance.qzv

qiime feature-table filter-features \
  --i-table q2-picrust2_output_female/ko_metagenome.qza \
  --p-min-frequency 50 \
  --p-min-samples 4 \
  --o-filtered-table q2-picrust2_output_female/ko_metagenome_filtered.qza

# machine learning classifier: find features that are most discriminat between states
# can be used with any feature table 
# classifier for males
qiime sample-classifier classify-samples \
  --i-table q2-picrust2_output_male/ko_metagenome_filtered.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Host_disease \
  --p-random-state 666 \
  --p-n-jobs 5 \
  --output-dir ./sample-classifier-results_male/

# classifier for females
qiime sample-classifier classify-samples \
  --i-table q2-picrust2_output_female/ko_metagenome_filtered.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Host_disease \
  --p-random-state 666 \
  --p-n-jobs 5 \
  --output-dir ./sample-classifier-results_female/

#plot a heatmap with the results from the male classification (male)
qiime sample-classifier heatmap \
  --i-table q2-picrust2_output_male/ko_metagenome_filtered.qza \
  --i-importance sample-classifier-results_male/feature_importance.qza \
  --p-feature-count 100 \
  --m-sample-metadata-file metadata.tsv \
  --m-sample-metadata-column Host_disease \
  --p-group-samples \
  --p-method "ward"\
  --p-color-scheme "nipy_spectral" \
  --o-filtered-table q2-picrust2_output_male/ko_metagenome_important-feature-table-top30.qza \
  --o-heatmap q2-picrust2_output_male/ko_metagenome_important-feature-table-heatmap.qzv


# export the feature importance list (male)
qiime metadata tabulate \
  --m-input-file sample-classifier-results_male/feature_importance.qza \
  --o-visualization sample-classifier-results_male/feature_importance.qzv

#plot a heatmap with the results from the male classification (female)
qiime sample-classifier heatmap \
  --i-table q2-picrust2_output_female/ko_metagenome_filtered.qza \
  --i-importance sample-classifier-results_female/feature_importance.qza \
  --p-feature-count 100 \
  --m-sample-metadata-file metadata.tsv \
  --m-sample-metadata-column Host_disease \
  --p-group-samples \
  --p-method "ward"\
  --p-color-scheme "nipy_spectral" \
  --o-filtered-table q2-picrust2_output_female/ko_metagenome_important-feature-table-top30.qza \
  --o-heatmap q2-picrust2_output_female/ko_metagenome_important-feature-table-heatmap.qzv


# export the feature importance list (female)
qiime metadata tabulate \
  --m-input-file sample-classifier-results_female/feature_importance.qza \
  --o-visualization sample-classifier-results_female/feature_importance.qzv
