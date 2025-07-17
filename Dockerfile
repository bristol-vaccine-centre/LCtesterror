# Base image with R and build tools
FROM rocker/r-ver:4.3.2

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (e.g., for R packages that use C++ or XML)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    pandoc \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages listed in DESCRIPTION using remotes
#RUN R -e "install.packages('remotes')"

# Copy your package to the image
#WORKDIR /usr/local/src
#COPY . LCtesterror
COPY . /usr/local/src/LCtesterror
WORKDIR /usr/local/src/LCtesterror


# Install the package and its dependencies
#RUN R -e "install.packages('future'); library(future); print('future installed OK')" && \
#    R -e "install.packages('remotes')" && \
#   R -e "remotes::install_deps('LCtesterror', dependencies = TRUE)" && \
#    R -e "remotes::install_local('LCtesterror')"

# Install renv and restore locked dependencies
# Confirm files exist
RUN ls -lah /usr/local/src/LCtesterror/
# Install renv and restore environment
RUN R -e "install.packages('renv')"
COPY renv.lock renv.lock
ENV RENV_PATHS_LIBRARY=renv/library
RUN R -e "renv::restore()"
#RUN R -e "renv::install('.')"

#install LCtesterror
RUN R -e "install.packages('devtools')" && \
    R -e "devtools::install_local('.', dependencies = TRUE)"

#check that key packages work
RUN R -e "library(future); library(LCtesterror); print('All good!')"

# Set default command
CMD ["R"]

