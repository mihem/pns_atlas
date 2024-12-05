# Use R 4.3.1 with ubuntu 22.02 jammy
# https://hub.docker.com/layers/rstudio/r-base/4.3.1-jammy/images/sha256-82c069a26ec53b50226f68c7d3e2271e60aff857207c5f9784355d660a10f3ef?context=explore
FROM rstudio/r-base@sha256:82c069a26ec53b50226f68c7d3e2271e60aff857207c5f9784355d660a10f3ef 

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update --yes \
 && apt-get upgrade --yes \
 && apt-get install --yes \
libglpk-dev libxml2-dev git libhdf5-dev \
libgsl-dev cmake libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
libtiff-dev libudunits2-dev libproj-dev libgdal-dev libgeos-dev


# install R packages required 
# Change the packages list to suit your needs
RUN R -e 'install.packages("https://packagemanager.posit.co/cran/latest/src/contrib/Archive/renv/renv_1.0.3.tar.gz")'

RUN echo "options(Ncpus = 6)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site

# Copy renv files 
WORKDIR /pns_atlas
COPY renv.lock renv.lock
COPY .gitignore .gitignore

# Persist the renv library and cache directories across builds
# VOLUME /root/.local/share/renv/cache /root/.local/share/renv/library

ENV RENV_PATHS_LIBRARY renv/library

# Restore the R environment
RUN R -e "renv::restore()"


