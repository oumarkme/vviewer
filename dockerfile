
FROM continuumio/miniconda3:25.3.1-1


# Update conda
RUN conda update -n base -c defaults conda

# Install dependencies
RUN conda install -c conda-forge -c bioconda -c r -c defaults \
    conda-forge::r-essentials=4.4 \
    conda-forge::r-base=4.4.3 \
    conda-forge::r-shiny=1.11.0 \
    conda-forge::r-data.table=1.17.6 \
    conda-forge::r-tidyverse=2.0.0 \
    bioconda::bioconductor-rsamtools=2.22.0 \
    conda-forge::r-bslib=0.9.0 \
    bioconda::bioconductor-genomicranges=1.58.0 \
    bioconda::bioconductor-iranges=2.40.0 \
    conda-forge::r-httr=1.4.7 \
    conda-forge::r-jsonlite=2.0.0 \
    conda-forge::r-xml2=1.3.8 \
    && conda clean --all -y

# Set the working directory
WORKDIR /

# COPY the app files
COPY app.R /app.R

# Expose the port for the Shiny app
EXPOSE 9876

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/app.R', host='0.0.0.0', port=9876, launch.browser=FALSE)"]
