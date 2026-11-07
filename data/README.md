# Input data

The pipeline reads these files. They are archived in the data deposit
(Zenodo, [10.5281/zenodo.24044208](https://doi.org/10.5281/zenodo.24044208)); download them into this folder (or point
`CCE1_DATA_DIR` at the download) with this layout:

```
data/
  CCE1_WBAT_HARP_daily.mat         daily CCE1 record: CTD temperature/salinity (8 depths), dissolved
                                   oxygen, nitrate, chlorophyll a, FSLE, backscatter (NASC, day/night),
                                   cetacean detections, hydrographic periods
  CCE1_WBAT_hourly.mat             hourly WBAT backscatter (NASC.All, 12 x 50 m bins) and day/night windows,
                                   used for DVM amplitude and, through `night_all`, for the layer
                                   tracing. Column 1 = shallowest bin, as in CCE1_WBAT_HARP_daily.mat;
                                   bin bounds (m below surface) in `depth`
  OMP/x_error_depth75_id1004.csv   OMP water-mass fractions and closure error at 75, 150 and 301 m
  OMP/x_error_depth150_id1003.csv
  OMP/x_error_depth301_id1002.csv
  MSL_night_seedWalk.csv           nightly mesopelagic scattering layer top edge, core depth and
                                   integrated backscatter (output of the MATLAB seed-and-walk detector,
                                   matlab/run_MSL_seedwalk.m)
  WBAT_export_5m/CCE1_WBAT_01_YYMM.csv   hourly, 5 m-layer WBAT exports (Echoview), one file per
                                   month, 1610-1711; layers ordered shallowest first, depths in m below
                                   the surface (transducer depth 504 m). Not read by the R scripts;
                                   input to matlab/run_MSL_seedwalk.m
```

Variable names and units for every file are in the `README.md` of the data deposit.
