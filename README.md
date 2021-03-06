# RIA-network-analysis
Data and R code for relative importance analysis (RIA) and network analysis: Carignan, C. (2019), "A network-modeling approach to investigating individual differences in articulatory-to-acoustic relationship strategies." Speech Communication, 108, 1-14, DOI: 10.1016/j.specom.2019.01.007.


This repository provides the data and R functions necessary to replicate the results from Carignan (2019, Speech Communication).


Save the following data set into your R working directory: imitation_data.csv


Import the following functions into R: oral.nasal.matching (from oral_nasal_matching.R), RIA (from RIA.R), lm.cross.validation (from cross_validation.R), network.clusters (from network_clusters.R), ED.measurements (from ED_measurements.R)


Use the code provided in network_analysis.R to perform the relative importance analysis (RIA) and subsequent network clustering, and to plot the results.

Additionally, the code provided in network_analysis_allV.R will perform the entire analysis for all three vowel pairs, combined. This analysis was not part of the original study, and is therefore provided as a separate code base.
