# Dockerfile 07/07/2020
# Julia environment

ARG BASE_CONTAINER=jupyter/datascience-notebook
FROM $BASE_CONTAINER

MAINTAINER Nadanai Laohakunakorn <nadanai.laohakunakorn@ed.ac.uk>

RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'using Pkg; Pkg.add("DelimitedFiles")' && \
    julia -e 'using Pkg; Pkg.add("DataFrames")' && \
    julia -e 'using Pkg; Pkg.add("CSV")' && \
    julia -e 'using Pkg; Pkg.add("Interpolations")' && \
    julia -e 'using Pkg; Pkg.add("Plots")' && \
    julia -e 'using Pkg; Pkg.add("DifferentialEquations")' && \
    julia -e 'using Pkg; Pkg.add("ProgressBars")' && \
    julia -e 'using Pkg; Pkg.add("Statistics")' && \
    julia -e 'using Pkg; Pkg.add("PyCall")' && \
    julia -e 'using Pkg; Pkg.add("PyPlot")' && \
    julia -e 'using Pkg; Pkg.add("Sundials")' && \
    julia -e 'using Pkg; Pkg.precompile()'

WORKDIR $HOME
