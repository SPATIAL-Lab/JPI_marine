# JPI_marine

Code and data files associated with Bowen et al. (submitted) Joint Proxy Inversion methods paper.

Top-level functions are contained in jpi_mgca18O.R (runs main analyses for paper), CalibFit.R (runs independent regressions used as exploratory and in specifying priors), and MgCa_sw_model.R (runs independent inversion of seawater Mg/Ca data). All figures in the manuscript were created with plots.R.

1123_HOL.csv: Coretop foraminiferal d18O and Mg/Ca values from Elderfield et al. (2010) used in the downcore Mg/Ca calibration test

1123_LGM.csv: Last Glacial Maximum foraminiferal d18O and Mg/Ca values from Elderfield et al. (2010) used in the downcore Mg/Ca calibration test

birner_2016.csv: Raw foraminiferal d18O and Mg/Ca values from Birner et al. (2016) used in JPI analyses

birner_2016_interp.csv: Processed and interpolated foraminiferal d18O and Mg/Ca data and interpreted BWT and seawater d18O from Birner et al. (2016), used only for plotting and comparison with JPI results

C_d18O_calib.csv: Compiled coretop d18O calibration data for C. spp

CalibFit.R: Top-level code used to conduct Bayesian regression modeling of independent proxy calibration relationships, used in developing prior estimates of proxy model parameter values

elderfield_2012.csv: Raw foraminiferal d18O and Mg/Ca values from Elderfield et al. (2012) used in JPI analyses

elderfield_2012_interp.csv: Processed and interpolated foraminiferal d18O and Mg/Ca data and interpreted BWT and seawater d18O from Elderfield et al. (2012), used only for plotting and comparison with JPI results

helpers.R: Functions used in preparing data and plotting results

jpi_mgca18O.R: Top-level code used to conduct JPI analyses reported in the manuscript

Lear_combined: Raw foraminiferal d18O and Mg/Ca values from Lear et al. (2002, 2015) used in JPI analyses

Lear_combined_interp.csv: Processed and interpolated foraminiferal d18O and Mg/Ca data and interpreted BWT and seawater d18O from Lear et al. (2002, 2015), used only for plotting and comparison with JPI results

mg_model.R: BUGS model used to run JPI on seawater Mg/Ca data independently

mgca_calib_model.R: BUGS model used to run Bayesian regression of coretop foraminiferal Mg/Ca data

mgca_dc_calib_model.R: BUGS model used to run Bayesian regression of downcore foraminiferal Mg/Ca data

mgca_sw.csv: Compiled proxy estimates of seawater Mg/Ca

MgCa_sw_model.R: Top-level code used to conduct independent proxy inversion of seawater Mg/Ca data

O_mgca_calib.csv: Compiled coretop Mg/Ca calibration data for O. spp

plots.R: Code used to create plots presented in the manuscript

split_temporal_birn.R: BUGS model used to run JPI on site U1385 proxy data

split_temporal_elder.R: BUGS model used to run JPI on site 1123 proxy data

split_temporal_elder_dc.R: BUGS model used to run JPI on site U1385 proxy data using downcore Mg/Ca calibration method

split_temporal_lear.R: BUGS model used to run JPI on site 806 proxy data

split_temporal_multi.R: BUGS model used to run simultaneous JPI on site U1385 and 1123 proxy data

U_d18O_calib.csv: Compiled coretop d18O calibration data for U. spp

U_mgca_calib.csv: Compiled coretop Mg/Ca calibration data for U. spp
