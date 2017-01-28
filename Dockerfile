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

# And finally build the docker image and push it
### docker build -t khturner/tn-seq Tn-seq/Gislason-Cardona/
### sudo docker push khturner/tn-seq

FROM ubuntu:latest
MAINTAINER Keith H. Turner "khturner@gmail.com"
WORKDIR /root

# Edit sources to install R
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -

# Install key software and dependencies with apt-get
RUN apt-get -y update
RUN apt-get install -y r-base r-base-dev gdebi-core bowtie2 libtre-dev libtre5 zlib1g zlib1g-dev wget git awscli gawk libcurl4-openssl-dev libxml2-dev libssl-dev

# Install Rstudio Server for interactive analyses
RUN wget https://download2.rstudio.org/rstudio-server-1.0.136-amd64.deb
RUN gdebi --non-interactive rstudio-server-1.0.136-amd64.deb

# FQGrep for read searching/filtering
RUN git clone https://github.com/indraniel/fqgrep.git && cd fqgrep && make && mv fqgrep /usr/local/bin/ && cd ..

# Get flexbar for adapter trimming - do I need to move that library somewhere?
# Commented cause I don't think I actually need it...
# RUN wget https://github.com/seqan/flexbar/releases/download/v2.5.0/flexbar_v2.5_linux64.tgz && tar -xvzf flexbar_v2.5_linux64.tgz && mv flexbar_v2.5_linux64/flexbar /usr/local/bin/

# Install R packages
RUN echo 'install.packages(c("tidyverse", "mclust", "seqinr", "devtools"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R
RUN echo 'devtools::install_github("dgrtwo/fuzzyjoin")' >> /tmp/packages.R
RUN echo 'install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.4/bioc"); source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2")' >> /tmp/packages.R
RUN Rscript /tmp/packages.R

# Clone Tn-seq repo for necessary scripts
RUN git clone https://github.com/khturner/Tn-seq.git
# IN DEVELOPMENT - switch to dev branch - remove this when done developing
RUN cd Tn-seq && git checkout dockerize && chmod a+x trimmer && mv trimmer /usr/local/bin/

# Move to script directory
WORKDIR /root/Tn-seq/py
ENV PATH /root/Tn-seq/py:$PATH
