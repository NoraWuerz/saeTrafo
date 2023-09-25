# The R Package saeTrafo for Estimating unit-level Small Area Models under Transformations

The package saeTrafo supports estimating regional means based on the nested error regression model (Battese et al., 1988). Therefore, point estimation and mean squared error estimation (Prasad and Rao, 1990) for the classical model is offered. In addition to the classical model, the logarithmic and the data-driven log-shift transformation are provided. The core function NER_Trafo allows several options to enter population data: Either individual population data or only aggregates can be entered. If full population data is given, the method of  Molina and Martín (2018) is applied. Compared to other small area packages, these transformations are accessible in the absence of population micro-data. Only population aggregates (mean values, population sizes and preferably also covariances) need to be supplied. The methodology for point and mean squared error estimates is described in Wuerz et al. (2022) and is made available in a user-friendly way within saeTrafo.

# Contact
nora.wuerz@uni-bamberg.de
