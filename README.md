# BipolarMicrobiome

Microbiome paired-end sequencing data from the European Nucleotide Archive (ENA) project PRJNA544731 was downloaded and pre-processed with QIIME 2. 
To generate an input, the following command was ran on Linux terminal in the same directory as all the FASTQ files were downloaded:

```
conda activate qiime2-2022.2

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ../manifest\ file.tsv \
--output-path paired-end-demux.qza \
--input-format PairedEndFastqManifestPhred33V2
```

The manifest file was created according to the [QIIME2 page recomendations](https://docs.qiime2.org/2022.2/tutorials/importing/).
