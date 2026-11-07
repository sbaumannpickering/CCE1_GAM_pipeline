# Mesopelagic scattering layer (MSL) tracing

MATLAB code that produces `MSL_night_seedWalk.csv`, the nightly top edge, core
depth and integrated backscatter of the non-migratory mesopelagic scattering
layer, plus an interactive tool to hand-label the layer edge and a script that
scores the detector against those labels. `R/08_watermass_regression.R` reads
the CSV. The CSV is already in the data deposit, so the tracing is only needed
to regenerate it.

## Run

From MATLAB (R2025b; needs the Statistics and Machine Learning Toolbox):

```matlab
run('matlab/run_MSL_seedwalk.m')
```

Inputs are read from `data/` and the output is written to `results/MSL/`. Set
the environment variables `CCE1_DATA_DIR` and `CCE1_RESULTS_DIR` to use other
folders, as for the R scripts.

| Input (`data/`) | Used for |
|---|---|
| `WBAT_export_5m/CCE1_WBAT_01_YYMM.csv` | hourly NASC in 5 m layers, depth in m below the surface |
| `CCE1_WBAT_hourly.mat`, variable `night_all` | sunset and following sunrise of each of the 386 nights |

Output: `results/MSL/MSL_night_seedWalk.csv` (columns in the `README.md` of the data deposit).
Run on the deposited inputs it reproduces the deposited `MSL_night_seedWalk.csv` to within 1e-12.

## Label the layer edge and score the detector

```matlab
run('matlab/run_MSL_labeling.m')   % interactive; needs a MATLAB window
run('matlab/run_MSL_scoring.m')
```

`run_MSL_labeling.m` opens `labelMSLTopEdgeGUI.m` on 16 windows of 50 stitched
hours (about 5 nights each), spread evenly over the record. In TRACE mode,
left-click adds a point on the top edge and right-click removes the nearest
point. In MARK ABSENT mode, left-click-drag shades a span with no clear
boundary. Progress is saved whenever you move between windows and when you
close the window, and running the script again resumes at the first unlabeled
window. The window shows real calendar dates and grey dotted lines at the real
night boundaries. Labels are saved to `results/MSL/MSL_TopEdge_HandLabels.mat`
(set `CCE1_LABEL_FILE` to use another path).

`run_MSL_scoring.m` runs the detector with the settings in
`mslDetectorConfig.m` and compares its edge with your labels: mean absolute
error, bias (positive = detector too deep), RMS error, the miss rate inside
traced windows and the false-positive rate inside spans marked absent. It writes
`MSL_score_perWindow.csv`, `MSL_score_absentIntervals.csv` and `MSL_score.mat` to
`results/MSL/`. Change `mslDetectorConfig.m` and score again to test a parameter
change. The hand labels used in the paper are not part of this release.

## Files

| File | Purpose |
|---|---|
| `run_MSL_seedwalk.m` | Driver: grid, night-time selection, detector, per-night aggregation, gap filling, filtering, CSV |
| `run_MSL_labeling.m` | Opens the hand-labeling tool |
| `run_MSL_scoring.m` | Scores the detector against the hand labels |
| `labelMSLTopEdgeGUI.m` | Interactive labeling tool |
| `scoreMSLEdgeAgainstLabels.m` | Error, miss rate and false-positive rate of a detector against labels |
| `buildLabelWindows.m`, `initLabels.m` | Evenly spaced label windows and the empty label structure |
| `renderEchogramLabelStyle.m`, `pcolorLikeImagesc.m` | Echogram drawing for the labeling tool |
| `ensureDeploymentID.m` | Treats the record as one continuous segment when no deployment ID is given |
| `loadMSLNightGrid.m` | Builds the night-time grid used by all three run scripts |
| `mslDetectorConfig.m` | Detector settings |
| `buildNASC5mGrid.m` | Reads the monthly 5 m exports into one hourly NASC x depth grid; drops Layer −100 (the 5 m bin at the transducer face, 499–504 m depth, ringdown) and Layers 1 and higher (past the sea surface) |
| `restrictToNightTrimmed.m` | Keeps night-time rows, trimmed by 1 h at both ends of each night, stitched into one continuous series |
| `computeMSLEdgeAndCore.m` | Seed-and-walk detector for the layer top edge and core |
| `interpolateShortGaps.m` | Fills gaps of up to 2 nights |

## Method

1. **Grid.** 5 m layers from 1.5 to 496.5 m depth, hourly, 2016-10-13 to
   2017-11-03.
2. **Night-time rows.** Rows from 1 h after sunset to 1 h before sunrise, joined
   into one series so the detector can follow the layer across nights.
3. **Detector.** For each hour the layer core is the maximum of the smoothed
   NASC profile between 100 and 495 m. It must exceed z = 4 against the
   reference water (50–100 m) and the 50th percentile of positive NASC in the
   search window. The top edge is found by walking toward the surface from the
   core while z ≥ 2.5 and NASC ≥ 0.85 × the floor. A continuity gate, a timeout
   and a re-acquisition rule govern tracking between hours.
4. **Per night.** Median edge depth, core depth, core NASC and integrated NASC
   over the hours with a detection (`median*` columns), with the number of hours
   present.
5. **Gap filling.** Gaps of up to 2 nights that do not touch the ends of the
   record are interpolated (`*Filled` columns).
6. **Filtering.** Five-night running median of the gap-filled series
   (`filtered*` columns). The R analysis uses the raw `median*` columns.
