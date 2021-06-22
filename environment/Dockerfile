# hash:sha256:809e485e54bbd6a584c926eca69b9d2abd740902cff525d618258cd63d35d5d2
FROM registry.codeocean.com/codeocean/r-base:4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
        rstudio-server=1.2.5033 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("DT", "0.17")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-2")' \
    && Rscript -e 'remotes::install_version("cowplot", "1.1.1")' \
    && Rscript -e 'remotes::install_version("cutpointr", "1.0.32")' \
    && Rscript -e 'remotes::install_version("data.table", "1.13.6")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.0.3")' \
    && Rscript -e 'remotes::install_version("egg", "0.4.5")' \
    && Rscript -e 'remotes::install_version("ggdendro", "0.1.22")' \
    && Rscript -e 'remotes::install_version("ggplot2", "3.3.3")' \
    && Rscript -e 'remotes::install_version("ggpubr", "0.4.0")' \
    && Rscript -e 'remotes::install_version("ggridges", "0.5.3")' \
    && Rscript -e 'remotes::install_version("ggsci", "2.9")' \
    && Rscript -e 'remotes::install_version("ggthemes", "4.2.0")' \
    && Rscript -e 'remotes::install_version("glmnet", "4.1")' \
    && Rscript -e 'remotes::install_version("pROC", "1.17.0.1")' \
    && Rscript -e 'remotes::install_version("purrr", "0.3.4")' \
    && Rscript -e 'remotes::install_version("reshape2", "1.4.4")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-7")' \
    && Rscript -e 'remotes::install_version("survminer", "0.4.8")' \
    && Rscript -e 'remotes::install_version("tidyverse", "1.3.0")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "AnnotationHub", \
        "BSgenome.Hsapiens.UCSC.hg19", \
        "Biobase", \
        "DESeq2", \
        "GenomicRanges", \
        "Homo.sapiens", \
        "Repitools", \
        "TCGAbiolinks", \
        "biovizBase", \
        "ggbio", \
        "limma", \
        "minfi" \
    ))' # Original versions: 2.22.0 1.4.3 2.50.0 1.30.0 1.42.0 1.3.1 1.36.0 2.18.0 1.38.0 1.38.0 3.46.0 1.36.0
