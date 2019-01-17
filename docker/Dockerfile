# Pull base image
FROM debian:testing
# FROM ubuntu:18.04

# Prepare directory structure
RUN mkdir /tools

# Change working directory to tools
WORKDIR /tools

# update path
ENV PATH="/tools:${PATH}"

### install software ###

# general software
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3.6 \
    # python-pip \
    default-jre \
    curl \
    alien \
    git \
    wget \
    mercurial \
    unzip \
    bzip2 \
    make \
    gcc \
    g++ \
    cmake \
    libncurses5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libtbb-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libcurl4-openssl-dev \
    ed \
    less \
    locales \
    ca-certificates \
    fonts-texgyre \
    hdf5-tools \
	libhdf5-dev \
	hdf5-helpers \
	libhdf5-serial-dev \
    apt-utils \
    nano

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
        && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
## Also install stringr to make dococt install (from source) easier
RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends \
		littler \
        r-cran-littler \
        r-cran-stringr \
		r-base \
		r-base-dev \
		r-recommended \
        && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# ENV R_BASE_VERSION 3.5.2
#
# ## Now install R and littler, and create a link for littler in /usr/local/bin
# ## Also set a default CRAN repo, and make sure littler knows about it too
# ## Also install stringr to make dococt install (from source) easier
# RUN apt-get update \
# 	&& apt-get install -t unstable -y --no-install-recommends \
# 		littler \
#         r-cran-littler \
#         r-cran-stringr \
# 		r-base=${R_BASE_VERSION}-* \
# 		r-base-dev=${R_BASE_VERSION}-* \
# 		r-recommended=${R_BASE_VERSION}-* \
#         && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /etc/R/Rprofile.site \
#         && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
# 	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
# 	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
# 	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
# 	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
# 	&& install.r docopt \
# 	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
# 	&& rm -rf /var/lib/apt/lists/*

# get google cloud sdk
WORKDIR /tools
RUN wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz && \
    tar -xzf google-cloud-sdk.tar.gz && \
    rm google-cloud-sdk.tar.gz && \
    google-cloud-sdk/install.sh
ENV PATH="${PATH}:/tools/google-cloud-sdk/bin"

# pigz (parallel gzip)
WORKDIR /tools
RUN wget http://www.zlib.net/pigz/pigz-2.4.tar.gz && \
    tar -xzf pigz-2.4.tar.gz && \
    rm pigz-2.4.tar.gz && \
    cd pigz-2.4 && \
    make
ENV PATH="${PATH}:/tools/pigz-2.4"

# fastqc
WORKDIR /tools
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    rm fastqc_v0.11.5.zip && \
    chmod 755 ./FastQC/fastqc
ENV PATH="${PATH}:/tools/FastQC"

# # htslib
# RUN git clone https://github.com/samtools/htslib && \
#     cd htslib && \
#     autoheader && \
#     autoconf && \
#     ./configure && \
#     make && \
#     make install
# ENV PATH="${PATH}:/tools/htslib"

# kallisto
WORKDIR /tools
RUN git clone https://github.com/pachterlab/kallisto.git && \
    cd kallisto && \
    autoreconf ext/htslib && \
    mkdir build && \
    cd build && \
    cmake .. && \
	make && \
	make install
ENV PATH="${PATH}:/tools/kallisto"

# pip
RUN wget https://bootstrap.pypa.io/get-pip.py && \
    python3.6 get-pip.py

# multiqc
RUN python3.6 -m pip install setuptools && \
    python3.6 -m pip install setuptools multiqc

# bsl2fastq
WORKDIR /tools
RUN wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    unzip bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    rm bcl2fastq2-v2-20-0-linux-x86-64.zip && \
    alien -k bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
    rm bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && \
    dpkg -i bcl2fastq2_0v2.20.0.422-1_amd64.deb && \
    rm bcl2fastq2_0v2.20.0.422-1_amd64.deb

# doing this again because otherwise apt-get doesn't seem to work interactively
# NOTE: should be unnecessary when using the docker with cromwell
RUN apt-get update -y

# install R packages
RUN mkdir -p /scripts
ADD ./R_scripts/install_packages.R /scripts
RUN Rscript /scripts/install_packages.R

# add R scripts
ADD ./R_scripts/get_kallisto_counts.R /scripts
ADD ./R_scripts/run_edgeR.R /scripts

# give 755 access (all to user, read/execute to group and other)
RUN chmod 755 /tools/*
