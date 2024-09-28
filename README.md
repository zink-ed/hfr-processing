# hfr-processing
Python toolbox to read and process High Frequency Radar (HFR) files in CODAR radial format. It can also combine multiple radial files into a processed Total file.

## Goals:
-	to observe surface currents
-	delayed-mode quality control of radial data utilizing QARTOD manual
-	create totals

### Quality Control Process:
1.	review daily and weekly radial distributions
2.	review radial diagnostics to identify any time periods that may require special attention
3.	plot all the spectra and check each timestep to ensure that the first order portion of doppler spectra are properly defined as identified by the first order line settings
a.	identify any sources of outside interference that may affect data quality
4.	recalculate the radial currents from the Doppler spectra with the best available measured antenna patterns and first order line settings to produce the best radial current vectors
5.	apply full suite of QUARTOD radial tests to reprocessed radials + flag vectors
6.	plot post-processed radial currents with original radial files from remote site

### Bilinear Interpolation Algorithm:
-	interpolating two dimensions by doing repeated linear interpolation (average of neighbors)



## Research Papers:

An Initial Evaluation of High-Frequency Radar Radial Currents in the Straits of Florida in Comparison with Altimetry and Model Products: https://digital-library.theiet.org/content/books/10.1049/sbra537e_ch5#:~:text=An%20analysis%20of%20the%20initial%20seven%20months%20of%20High-Frequency%20Radar

A Unified Approach to HF Radar Radial Quality Control for Understanding Gulf Ocean Systems: https://ieeexplore.ieee.org/document/9705910



## Resources: 

Information about SeaSondes: http://support.codar.com/Technicians_Information_Page_for_SeaSondes/Docs/GuidesToFileFormats/File_LonLatUV_RDL_TOT_ELP.pdf#:~:text=LLUV%20are%20allowed%20to%20have%20more%20than%20one,%E2%80%9CRD%E2%80%9D%20will%20be%20read%20in%20as%20radial%20data.

QARTOD: https://cdn.ioos.noaa.gov/media/2022/07/HFR_QARTOD_Manual_Update_Final-1b.pdf

