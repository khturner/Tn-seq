#!/bin/bash
# This is meant to be run on a fresh EC2 instance (from the home directory!) to prepare for TnSeq analyses for April's project

# Edit sources to install R
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -

# Install key software and dependencies with apt-get
sudo apt-get update
sudo apt-get install r-base r-base-dev \ # R is king
gdebi-core \ # Dependency for Rstudio Server
bowtie2 \ # For mapping
libtre-dev libtre5 zlib1g zlib1g-dev \ # Dependencies for fqgrep
wget git # Other tools

# Install Rstudio Server for interactive analyses
wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
sudo gdebi rstudio-server-1.0.136-amd64.deb

# FQGrep for read searching/filtering
git clone https://github.com/indraniel/fqgrep.git
cd fqgrep
make
cd ..

# Get flexbar for adapter trimming - IN PROGRESS
wget https://github.com/seqan/flexbar/releases/download/v2.5.0/flexbar_v2.5_linux64.tgz

