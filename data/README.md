# Input data

The pipeline reads these files. They are archived in the data deposit
(**Dryad DOI: to be added**); download them into this folder (or point
`CCE1_DATA_DIR` at the download) with this layout:

```
data/
  combined_diel_daily.mat              daily CCE1 record: CTD temperature/salinity (8 depths), dissolved
                                       oxygen, nitrate, chlorophyll a, FSLE, backscatter (NASC, day/night),
                                       cetacean detections, hydrographic periods
  CCE1_SB_01_all.mat                   hourly WBAT backscatter (NASC.All) used for DVM amplitude
  OMP/x_error_depth75_id1004.csv       OMP water-mass fractions and closure error at 75, 150 and 301 m
  OMP/x_error_depth150_id1003.csv
  OMP/x_error_depth301_id1002.csv
  MSL_final_per_night_seedWalk.csv     nightly mesopelagic scattering layer top edge, core depth and
                                       integrated backscatter (output of the MATLAB seed-and-walk detector)
```

Variable names and units are described in the deposit's data dictionary.
