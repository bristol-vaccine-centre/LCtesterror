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
RUN R -e "install.packages('remotes')"

# Copy your package to the image
WORKDIR /usr/local/src
COPY . LCtesterror

# Install the package and its dependencies
RUN R -e "install.packages('future'); library(future); print('future installed OK')" && \
    R -e "install.packages('remotes')" && \
    R -e "remotes::install_deps('LCtesterror', dependencies = TRUE)" && \
    R -e "remotes::install_local('LCtesterror')"

# Set default command
CMD ["R"]


