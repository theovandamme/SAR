# SAR
the scripts form my bachelorpaper used for the local analysis using SAR

Helper, Border_noise_correction, speckle_filter and terrain_flattening are part of a preprocessing framework for SAR amplitude imagery designed by Mullisa et al. (2021)
Version: v1.2
Date: 2021-02-11
Authors: Mullissa A., Vollrath A., Braun, C., Slagter B., Balling J., Gou Y., Gorelick N.,  Reiche J.

Mullissa, A., Vollrath, A., Odongo‐Braun, C., Slagter, B., Balling, J., Gou, Y., Gorelick, N., & Reiche, J. (2021). Sentinel‐1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine [Number: 10 Publisher: Multidisciplinary Digital Publishing Institute]. Remote Sensing, 13(10), 1954. https://doi.org/10.3390/rs13101954

Wrapper and CD_thresholding are used for the Change detection with thresholding methodology. Wrapper is an adapted version of the script by Mullisa et al. (2021) while CD_thresholding is an adapted version of a script designed by Michielsen (2022)

Michielsen, M. (2022). Flooding hazards in the plain of the ruzizi river studied through radar and optical satellite imagery [Master’s thesis, Vrije Universiteit Brussel].

Wrapper_SAR and SAR_indices are used for the calculation of a pre- and post-quality mosaic of a chosen index. It also returns the relative difference between the two mosaics.
