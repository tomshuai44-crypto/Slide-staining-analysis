# H-DAB Dominance Extractor

## What was installed

The ImageJ/Fiji plugin JAR is installed at:

```text
/Users/tom/Fiji/plugins/H_DAB_Dominance_Extractor.jar
```

After restarting Fiji, it should appear at:

```text
Plugins > Analyze > H-DAB Dominance Extractor
```

## Method

The plugin uses fixed H-DAB color deconvolution in optical-density space, then applies an optimized DAB dominance gate plus a high-confidence DAB cloud floor tuned for the three supplied images:

```text
DAB dominance = normalized DAB - alpha * normalized Hematoxylin
```

The positive mask is:

```text
tissue
AND normalized DAB >= optimized DAB threshold
AND DAB dominance >= optimized dominance threshold
AND normalized DAB >= 1.15 * normalized Hematoxylin
AND normalized residual <= residual threshold
```

The dominance score is used only to define the DAB-positive region. DAB intensity measurements are computed from the original DAB concentration channel, not from the subtracted dominance score.

For the supplied real images, the automatic DAB threshold is also constrained to at least the 85th percentile of DAB-dominant tissue pixels. This makes the plugin a high-confidence DAB detector for these low/intermediate DAB images, rather than a broad weak-stain detector.

## Outputs

For each RGB H-DAB image the plugin can produce:

```text
*_hematoxylin_color.tif
*_dab_color.tif
*_residual_color.tif
*_hematoxylin.tif
*_dab.tif
*_residual.tif
*_dab_dominance.tif
*_dab_positive_mask.tif
*_measurements.csv
```

The three color deconvolution images requested by the analysis are:

```text
Hematoxylin color
DAB color
Residual color
```

The grayscale concentration images are kept for measurement, because concentration maps are better for quantitative analysis. The extra dominance image and positive mask are used for DAB-positive region detection and QC.

## Smoke test evidence

The no-GUI Java smoke test ran on:

```text
/Users/tom/Fiji/cloud_extract_eval/results/example_rgb.png
```

It wrote outputs to:

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_smoke
```

Generated files:

```text
smoke_hematoxylin_color.tif
smoke_dab_color.tif
smoke_residual_color.tif
smoke_hematoxylin.tif
smoke_dab.tif
smoke_residual.tif
smoke_dab_dominance.tif
smoke_dab_positive_mask.tif
smoke_measurements.csv
```

## Synthetic benchmark result

Synthetic H-DAB benchmark:

```text
/Users/tom/Fiji/cloud_extract_eval/results/report.md
/Users/tom/Fiji/cloud_extract_eval/results/benchmark_summary.csv
```

Overall comparison:

| Method | Macro F1 | DAB F1 | Background false stain |
|---|---:|---:|---:|
| OD k-means | 0.569 | 0.553 | 36.151% |
| RGB k-means | 0.567 | 0.600 | 36.151% |
| CCE-stain ROI vectors | 0.499 | 0.507 | 0.000% |
| Fixed H-DAB deconvolution | 0.495 | 0.508 | 0.000% |
| Proposed H-DAB dominance, strict real-image mode | 0.476 | 0.236 | 0.000% |
| CCE-cloud ROI Gaussian | 0.361 | 0.355 | 0.226% |

The strict high-confidence mode is intentionally conservative and is tuned to the three supplied real images. Earlier balanced tuning had higher synthetic F1, but over-called DAB on the real low/intermediate-stain slides. The installed plugin therefore prioritizes the real-image objective over broad synthetic recall.

## Real three-image comparison

Inputs:

```text
/Users/tom/Desktop/NT_HTXR_K2_J1.jpg
/Users/tom/Desktop/slide-HFLday8-IRmodel-2026-03-26T12-37-47-R1-S10_11.3x.jpg
/Users/tom/Desktop/CB1_PCOS_Carotis_132_1_reexport.jpg
```

Plugin outputs:

```text
/Users/tom/Fiji/cloud_extract_eval/real_plugin_outputs
```

Comparison report:

```text
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_report.md
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_qc_summary.csv
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_qc_metrics.csv
```

Because no manual ground-truth labels were supplied, the real-image comparison uses QC metrics:

```text
white-background false positive rate
black-artifact false positive rate
DAB-positive to non-positive DAB ratio
mean H leakage inside DAB-positive mask
mean residual inside DAB-positive mask
```

Real-image summary:

| Method | QC score | DAB area fraction | White positive % | Artifact positive % | DAB pos/non-pos ratio |
|---|---:|---:|---:|---:|---:|
| Proposed H-DAB dominance | 1.483 | 0.097 | 0.000 | 0.000 | 1.973 |
| CCE-style auto seeds | 1.403 | 0.118 | 0.000 | 0.000 | 1.933 |
| Fixed H-DAB deconvolution | 1.159 | 0.293 | 0.000 | 0.000 | 1.814 |
| OD k-means | 0.335 | 0.095 | 0.000 | 0.000 | 0.549 |
| RGB k-means | 0.314 | 0.120 | 0.000 | 0.000 | 0.546 |

On the three supplied images, the optimized plugin has the best unsupervised QC score and keeps estimated background/artifact positives at zero.

## Current limitation

The real three images do not include pathologist/manual DAB-positive ground-truth labels. The plugin outperforms the tested baselines by the unsupervised QC criteria above, but formal accuracy/F1 on the real images would require either:

```text
manual/pathologist labels
or representative positive/negative ROIs
or an agreed QC target such as maximum background false positive plus DAB-positive area agreement
```
