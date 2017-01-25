#!/bin/bash
# This is a work in progress, don't run as a shell script!

# get the docker image, hop into it
docker pull khturner/tn-seq
docker run -it khturner/tn-seq bash

# Start installing some dependencies
cd
apt-get install python-setuptools python-dev build-essential
easy_install pip
pip install -U FastaIndex

# Install the LAST dependency
wget http://last.cbrc.jp/last-830.zip
unzip last-830.zip
cd last-830
make
make install
cd

# Move into a working directory
mkdir scaffolding
cd scaffolding
git clone https://github.com/lpryszcz/pyScaf.git

# Get the K56 genome
aws s3 configure # put in aws access keys and the region is us-east-1
aws s3 cp s3://cenocepacia-tnseq/k56.fasta ./

# Next up, download the J2315 genome and try scaffolding it - IN PROGRESS
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/485/GCF_000009485.1_ASM948v1/GCF_000009485.1_ASM948v1_genomic.fna.gz
gunzip GCF_000009485.1_ASM948v1_genomic.fna.gz
pyScaf/pyScaf.py -f k56.fasta -t 4 --dotplot png -g 1000000 -r GCF_000009485.1_ASM948v1_genomic.fna -o k56_scaffolded.fna
