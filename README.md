# Lreuteri_strain_analysis

## Introduction
This repository contains codes and example data to follow workflows used in "Metagenomic association analysis of gut symbiont _Lactobacillus reuteri_ without host-specific genome isolation" manuscript.

## Requirement
You can build a conda environment with necessary libraries and programs as following:

conda create -n reuteri_strain -c anaconda -c bioconda -c conda-forge -c etetoolkit ete3 pandas numpy scipy biopython raxml prodigal mummer prokka roary mmseqs2 bracken eggnog-mapper blast;

conda install -n reuteri_strain -c bioconda kraken --force-reinstall;

## Script

#Pick genome-based strain types (GSTs) and build Kraken database from reference genome sequences

python createGSTdb.py -i genome_info.txt -o pick_GST -d KrakenGST

#Build pan-genome database from reference genome sequences

python createPangenome.py -i genome_info.txt -o Pangenome

#Profile GST abundance of metagenomic samples

python profileGST.py -i sample_info.txt -d KrakenGST -o gst_profile

#Profile gene composition of metagenomic samples

python profileGene.py -i sample_info.txt -d Pangenome -o gene_profile
