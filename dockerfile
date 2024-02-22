
FROM rocker/rstudio:4.0.0
#FROM rocker/rstudio:3.4.4

ENV RENV_VERSION=0.9.3-106

RUN apt update
RUN apt install -y libz-dev libxml2-dev

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# should be in the same directory as this file
COPY renv.lock ./
RUN R -e 'renv::consent(provided = TRUE)'
RUN R -e 'renv::restore()'

