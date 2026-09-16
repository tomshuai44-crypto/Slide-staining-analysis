# H-DAB Dominance Extractor: local tissue model installed

The default tissue model is now **Local H+DAB transitions**, version `local-hdab-1.0`, with the reviewed threshold **k=8**. It is implemented entirely in Java inside the Fiji plugin; Python is not needed to analyze an image in Fiji. The previous `h-spatial-1.2` method remains available as **Previous H-only bounded growth** for reproducible comparisons.

Restart an already running Fiji instance before using the new version. The installed JAR is `/Users/tom/Fiji/plugins/H_DAB_Dominance_Extractor.jar`. The menu remains **Plugins → Analyze → H-DAB Dominance Extractor**.

## How to use it

1. Open an RGB image and optionally select an area ROI.
2. Choose **Local H+DAB transitions**, leave **Local tissue threshold k = 8** initially, and choose an output folder.
3. Inspect the saved tissue mask, green tissue overlay, and H seed evidence before interpreting the CSV. Green means selected tissue; the cyan/magenta experimental overlays were comparisons against the previous mask, not the plugin's standard display convention.
4. For a comparison with the former algorithm, choose **Previous H-only bounded growth**. Its reach setting applies only to that model. Local transitions have no maximum growth radius.

The model choice and local thresholds are recorded in every measurement CSV. To reproduce a tissue selection within an ROI, retain the same source image, ROI geometry and parameters. Only the current plane is analyzed; the slice, ROI pixel count and calibration are recorded. Spatial parameters remain native pixels, not automatically calibrated micrometres.

Example macro:

```java
run("H-DAB Dominance Extractor",
    "tissuemodel=local tissuelocalk=8 tissuelocalminod=0.02 " +
    "tissueautowhite=true tissuesupportradius=8 tissuecloseradius=2 " +
    "show=true output=[/absolute/output/folder] prefix=sample");
```

Use `tissuemodel=h_only` to request the former model explicitly. Existing macros that omit `tissuemodel` now use the new default. ImageJ macro keys must be lowercase.

## What the model does

The new model retains the previous algorithm's independently estimated background reference and supported H seeds. The H seed threshold, 3×3 cluster check and sliding support window are unchanged. It then replaces the old distance-limited growth with the reviewed local-transition rule:

1. **Consistent diagnostic H and DAB:** correct OD with the estimated or manual RGB reference; solve an H/DAB/neutral model; clip diagnostic H and DAB coefficients nonnegative. Neutral OD belongs to the neutral component. This diagnostic DAB differs from the extractor's legacy H/DAB/residual measurement channel.
2. **Independent background samples:** identify the brightest 10% of local brightness means and the lowest quartile of brightness standard deviations within that set. Sample up to 100,000 eligible pixels without replacement, with the configured fixed seed. This does not condition on the tissue mask or H distance.
3. **Local distributions:** compute median, Q90−Q10 spread, and fraction exceeding a background Q99 pixel threshold in sliding 17×17 windows by default. Median/spread windows are sampled every four native pixels and linearly interpolated; fractions are computed on the full grid.
4. **Feature thresholds:** for each local feature F, use `T(F)=median_background(F)+k*max(1.4826*MAD_background(F), floor)`. Floors are 0.01 for medians, 0.015 for spreads and 0.02 for fractions. k is an engineering multiplier, not calibrated confidence. A stain supports growth if its median exceeds its threshold, or both its spread and excess fraction exceed theirs.
5. **Spatial support:** retain local H or DAB evidence, or visibility evidence accompanied by H/DAB texture. Nearby reliable H sustains DAB-negative H tissue. High OD alone does not sustain growth. Require 60% support in a 5×5 neighborhood and visible material with corrected OD norm at least 0.02. A radius-2 closing bridges small gaps; it does not generally fill lumens.
6. **No distance cap:** keep every 8-connected allowed region touching a supported H seed. This is an unrestricted connectivity search through locally supported material. It is not a 96-pixel dilation. Isolated DAB material cannot create its own seed.

If H evidence is absent, or there are fewer than 64 usable local-background samples, selection fails explicitly. There is no automatic fallback to darkness, DAB positivity, or the previous model. Empty selections export an empty mask and NaN area fraction/intensity summaries. The former reference/seed failures retain their existing status labels.

## Controls and provenance

| Control / macro key | Default | Effect |
|---|---:|---|
| `tissuemodel` | `local` | New model; `h_only` selects previous method |
| `tissuelocalk` | 8 | Higher values require stronger departure from background |
| `tissuelocalminod` | 0.02 | Minimum visible material for local growth |
| `tissuesupportradius` | 8 px | 17×17 seed/local feature neighborhoods; local maximum radius is 32 |
| `tissuecloseradius` | 2 px | Small-gap closing radius |
| `tissueminod` or `tissue` | 0.08 | Corrected visibility threshold for H seed detection |
| `tissuehmin` | 0.08 | H seed threshold floor |
| `tissueseed` | 20260916 | Reproducible background/seed sampling |
| `tissuereachpixels` | 96 px | Applies only to `h_only`; ignored by `local` |

Other H seed/reference options remain available as documented in [the H-only method description](H_DAB_Tissue_Mask_Correction.md). The local CSV includes the algorithm/model, diagnostic basis, sampler identity, background pool/sample counts, feature thresholds, robust centers/spreads, pixel excess thresholds, local k, material threshold and sampling stride. `tissue_reach_px` is NaN and `tissue_distance_metric` is `UNBOUNDED_8_CONNECTED` in local mode. Former-model settings retained in the CSV do not imply those settings constrain local growth.

The output filenames remain the same: H/DAB/residual color and concentration maps, dominance score, DAB-positive mask, tissue mask, tissue overlay, H evidence mask, and measurements CSV. Use a new prefix/output folder to preserve older measurements when comparing models.

## Measurement behavior

The strict DAB-positive optimizer, its p85 floor, stain matrices, raw H/DAB/residual concentration calculations, and positive-pixel intensity definitions are unchanged. The final local tissue mask now supplies the p99 normalization domain, threshold optimization domain, positive-mask restriction and area denominator.

Consequently DAB-positive counts and fractions can change even though the strict-positive policy is unchanged. This is an intended consequence of selecting a different analysis area; replacing a denominator in a historical CSV is not equivalent to rerunning the plugin.

The local model uses DAB as evidence for tissue and can therefore make selected area depend on marker staining. This is the model the user authorized adopting. The experiment's DAB-disabled ablation showed that dependence clearly. H-rich/DAB-negative tissue is retained, but sparse-H/DAB-negative or acellular tissue remains a limitation. Colored artifact seeds can also remain selected. The k=8 choice was reviewed on NT and is not a universally validated threshold for other stains, resolutions or specimens.

## Native Java validation

The Java implementation matches the reviewed NT k=8 experiment with mask IoU **0.999974917**: 112 of 4,971,616 pixels differ. H seed masks are identical. Java selects **4,465,082** pixels versus the experiment's 4,465,194. The Java port uses a fixed-seed partial shuffle instead of NumPy PCG64 sampling, and floating-point library details differ. These differences are disclosed rather than claiming bit-for-bit equality with Python.

All **81 Java checks** pass against the installed JAR: 36 former-model regressions and 45 local-model checks. They include analytic tissue geometry at k=2/3/4/5/8, DAB-supported extension beyond 250 pixels from H, DAB-negative H tissue, lumens, blanks/noise, disconnected unseeded DAB, the known colored-artifact limitation, repeatability, ROI exclusion and outside-ROI invariance, empty/failure handling, CSV provenance and measurement consistency. The synthetic geometry has IoU 1.0 under all five k settings. These controlled tests are not real-image biological accuracy measurements.

On the same saved JPEG copies:

| Image | Previous tissue pixels | Local k=8 tissue pixels | Previous DAB fraction | Local DAB fraction |
|---|---:|---:|---:|---:|
| NT | 4,500,255 | 4,465,082 | 7.9858% | 8.0538% |
| CB1 | 2,885,412 | 2,632,655 | 9.2017% | 9.5853% |
| HFL | 1,731,101 | 1,769,431 | 9.1340% | 9.0903% |

The previous-model option reproduces the earlier NT tissue mask exactly. Raw H/DAB/residual array hashes match across models on all three images. The resulting masks and overlays are available in [the installation validation folder](local_transition_plugin_validation/REPORT.md). The images are saved JPEG copies with no physical calibration; no real annotated accuracy is claimed. CB1's reduced area particularly merits review of pale stroma when using k=8.

## Build, installation and rollback

```bash
python3 /Users/tom/Fiji/cloud_extract_eval/build_h_dab_tissue.py --install
python3 /Users/tom/Fiji/cloud_extract_eval/validate_h_dab_local.py
```

The build uses the bundled JDK, Java 8-compatible bytecode and ImageJ 1.54p. It packages only H-DAB classes and its menu entry, tests the build, backs up the previous JAR, atomically replaces the installed JAR, and tests the installed file. `local_transition_plugin_validation/build_manifest.json` records source/JAR hashes and the precise timestamped backup path. The original H-only JAR is also retained as `local_transition_plugin_validation/backup/H_DAB_Dominance_Extractor.h_spatial_1_2.jar`.

Selecting `h_only` is the quickest comparison path. To restore the prior installation completely, copy that backup JAR over `plugins/H_DAB_Dominance_Extractor.jar` and restart Fiji. Existing open Fiji instances are not hot-reloaded by copying a JAR.
