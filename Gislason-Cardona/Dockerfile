# This is meant to be run on a fresh t2.xlarge EC2 instance (from the home directory!) to prepare for TnSeq analyses for April's project
# Log into your new AWS instance and install Docker:
### sudo apt-get update && sudo apt-get install docker.io

# Then clone this repo:
### git clone https://github.com/khturner/Tn-seq.git

# Add your user to the docker group
### sudo gpasswd -a ${USER} docker

# Logout and log back in so your group membership is updated
# Then start up docker and login to Dockerhub
### sudo service docker start
### sudo docker login

# And finally build the docker image
### docker build -t khturner/Tn-seq Tn-seq/Gislason-Cardona/

FROM ubuntu:latest
MAINTAINER Keith H. Turner "khturner@gmail.com"

# Edit sources to install R
RUN sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | sudo apt-key add -

# Install key software and dependencies with apt-get
RUN sudo apt-get update
RUN sudo apt-get install r-base r-base-dev \ # R is king
gdebi-core \ # Dependency for Rstudio Server
bowtie2 \ # For mapping
libtre-dev libtre5 zlib1g zlib1g-dev \ # Dependencies for fqgrep
wget git # Other tools

# Install Rstudio Server for interactive analyses
RUN wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
RUN sudo gdebi rstudio-server-1.0.136-amd64.deb

# FQGrep for read searching/filtering
RUN git clone https://github.com/indraniel/fqgrep.git && cd fqgrep && make && cd ..

# Get flexbar for adapter trimming - IN PROGRESS
RUN wget https://github.com/seqan/flexbar/releases/download/v2.5.0/flexbar_v2.5_linux64.tgz

