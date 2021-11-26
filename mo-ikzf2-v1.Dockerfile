# pull base image
FROM cancerbits/dockr:4.0.2-a

# who maintains this image
LABEL maintainer mohamed shoeb "mohamed.shoeb@ccri.at"
LABEL version mo-ikzf2-v1

#install stack
RUN apt-get update
RUN curl -sSL https://get.haskellstack.org/ | sh

#install pandoc
WORKDIR ~  
RUN apt-get -qq install pandoc pandoc-citeproc
#install GNU scientific library (GSL)
RUN apt-get -qq install libgsl-dev

#update path
#RUN export PATH=/root/.local/bin:$PATH

#packages
RUN R -e "BiocManager::install(c('MAST'))"
RUN R -e "BiocManager::install(c('org.Hs.eg.db'))"
RUN R -e "BiocManager::install(c('EnrichmentBrowser'))"
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38'))"
RUN R -e "BiocManager::install('multtest')"
RUN install2.r --error kableExtra here SoupX sctransform metap gghighlight umap
Run R -e "remotes::install_github('satijalab/seurat', ref = 'release/4.0.0')"
Run R -e "remotes::install_github('jlmelville/uwot')"
Run R -e "remotes::install_github('mojaveazure/seurat-disk')"
