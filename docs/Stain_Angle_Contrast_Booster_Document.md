# Stain Angle Contrast Booster

## Purpose

This plugin is a visualization tool, not a quantitative stain measurement tool.

The goal is to make weakly separated brightfield stains easier for a human observer to see. It is useful when two stain colors are present but their optical density is low, their apparent RGB colors are close, or the observer has to zoom in heavily to notice the difference.

The plugin does three things:

1. Estimates or accepts two stain color vectors in optical density space.
2. Separates the image into two stain-density maps.
3. Re-renders those maps with an artificially larger angle between the two display colors.

The output should be interpreted as an enhanced viewing image. It should not be used as a direct measurement of stain positivity.

## Literature Basis

### Color Deconvolution

Ruifrok and Johnston introduced color deconvolution for histochemical stain separation. The key idea is to transform RGB intensities into optical density space, where stain mixtures are closer to linear combinations of stain vectors.

Reference:

- Ruifrok AC, Johnston DA. Quantification of histochemical staining by color deconvolution. Analytical and Quantitative Cytology and Histology. 2001. PubMed: https://pubmed.ncbi.nlm.nih.gov/11531144/

### OD-Space Stain Vector Estimation

Macenko et al. described a practical way to estimate stain vectors automatically from a slide. Their method converts RGB to optical density, removes background pixels, uses SVD/PCA to find the main OD plane, and takes angular extremes in that plane as stain vectors.

Reference:

- Macenko M et al. A method for normalizing histology slides for quantitative analysis. ISBI 2009. PDF: https://wwwx.cs.unc.edu/~mn/sites/default/files/macenko2009.pdf

### Perceptual Re-Coloring for Histology

Kather et al. showed that standard blue-brown histology color maps are not always optimal for human observers. They extracted a bivariate color map from histology images and optimized it to increase perceptual contrast. This directly supports the idea that re-rendering stain channels in a more separable color space can improve visual interpretation.

Reference:

- Kather JN et al. New Colors for Histology: Optimized Bivariate Color Maps Increase Perceptual Contrast in Histological Images. PLoS One. 2015. https://pmc.ncbi.nlm.nih.gov/articles/PMC4696851/

### Decorrelation Stretch

Decorrelation stretch is a general image-enhancement method that applies principal component analysis to separate correlated color information, stretches the components, and maps them back to RGB. It is common in remote sensing and is conceptually related to this plugin, although this plugin works in stain optical-density space rather than raw RGB.

Reference:

- MathWorks documentation, `decorrstretch`: https://www.mathworks.com/help/images/ref/decorrstretch.html

## Algorithm

### 1. Estimate White Reference

For each RGB channel, the plugin estimates a white reference from a high percentile of the image intensity distribution.

This avoids assuming that pure background is exactly `(255, 255, 255)`.

Default:

```text
white percentile = 99.8
```

### 2. Convert RGB to Optical Density

For each pixel:

```text
OD_R = -log((R + 1) / (white_R + 1))
OD_G = -log((G + 1) / (white_G + 1))
OD_B = -log((B + 1) / (white_B + 1))
```

Negative OD values are clipped to zero.

The plugin marks tissue pixels using:

```text
sqrt(OD_R^2 + OD_G^2 + OD_B^2) >= tissue_min_OD
```

Default:

```text
tissue_min_OD = 0.06
```

### 3. Get Two Stain Vectors

The plugin supports four modes:

```text
Auto OD-PCA
H-DAB
H&E
Custom OD vectors
```

In auto mode:

1. Build a 3 by 3 OD covariance matrix from tissue pixels.
2. Compute the first two eigenvectors.
3. Project tissue pixels onto this 2D OD plane.
4. Convert projected pixels to angular coordinates.
5. Use low and high angular percentiles as the two stain vectors.

Default:

```text
extreme percentile = 1.0
```

This is Macenko-style vector estimation, simplified for a visualization plugin.

### 4. Compute Stain Densities

The plugin completes the two stain vectors with a residual vector:

```text
v3 = normalize(cross(v1, v2))
```

Then it solves:

```text
OD = C1 * v1 + C2 * v2 + C3 * v3
```

Only `C1` and `C2` are used for visualization.

The outputs `*_stain1_density.tif` and `*_stain2_density.tif` are the raw concentration-like maps.

### 5. Local Contrast Enhancement

For each stain density map, the plugin combines global percentile stretching and local contrast.

For stain 1:

```text
global = C1 / percentile(C1, 95)
local = max(0, (C1 - local_mean(C1)) / (local_std(C1) + floor))
boosted = global + local_weight * 0.35 * local
visual = asinh(chroma_gain * boosted) / asinh(chroma_gain * 1.5)
visual = visual ^ density_gamma
```

Defaults:

```text
scale percentile = 95
local radius = 24 px
local weight = 0.70
chroma gain = 2.00
density gamma = 0.72
```

This is why the plugin is not just a global contrast adjustment. It preserves global weak stain information but also boosts local deviations.

### 6. Increase the Display Angle

The original stain vectors may be close:

```text
angle(v1, v2) = small
```

The plugin computes the bisector:

```text
axis = normalize(v1 + v2)
```

Then it constructs two new display vectors around that axis:

```text
u1 = cos(theta') * axis + sin(theta') * delta
u2 = cos(theta') * axis - sin(theta') * delta
```

where:

```text
theta' = target_display_angle / 2
target_display_angle = max(original_angle * angle_gain, minimum_display_angle)
```

Defaults:

```text
angle gain = 1.8
minimum display angle = 55 degrees
maximum display angle = 118 degrees
```

The vectors are clamped to valid non-negative OD directions so that the output remains displayable as RGB.

### 7. Re-Render Visualization

The angle-boosted overlay is rendered as:

```text
display_OD = display_strength * visual_C1 * u1
           + display_strength * visual_C2 * u2

display_RGB = 255 * exp(-display_OD)
```

The plugin also renders stain-specific images:

```text
*_stain1_enhanced_color.tif
*_stain2_enhanced_color.tif
```

## Output Files

For each input image:

```text
*_angle_boosted_overlay.tif
*_stain1_enhanced_color.tif
*_stain2_enhanced_color.tif
*_stain1_density.tif
*_stain2_density.tif
*_angle_boost_report.csv
```

Additional preview panels were generated during testing:

```text
*_preview_panel.jpg
```

## Installation

Installed plugin:

```text
/Users/tom/Fiji/plugins/Stain_Angle_Contrast_Booster.jar
```

Fiji menu:

```text
Plugins > Analyze > Stain Angle Contrast Booster
```

Source:

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_src/Stain_Angle_Contrast_Booster.java
```

Build command:

```bash
cd /Users/tom/Fiji
JDK="$PWD/java/macos-arm64/zulu21.42.19-ca-jdk21.0.7-macosx_aarch64/zulu-21.jdk/Contents/Home"
"$JDK/bin/javac" --release 8 -cp jars/ij-1.54p.jar \
  -d cloud_extract_eval/angle_plugin_build/classes \
  cloud_extract_eval/plugin_src/Stain_Angle_Contrast_Booster.java
cp cloud_extract_eval/plugin_src/angle_plugins.config \
  cloud_extract_eval/angle_plugin_build/classes/plugins.config
"$JDK/bin/jar" cf plugins/Stain_Angle_Contrast_Booster.jar \
  -C cloud_extract_eval/angle_plugin_build/classes .
```

## Test Results on the Three Provided Images

Output directory:

```text
/Users/tom/Fiji/cloud_extract_eval/angle_plugin_outputs_auto
```

### NT_HTXR_K2_J1

Auto mode estimated two stain vectors with an original angle of about `29.15 degrees`.

The display angle was increased to `55.00 degrees`.

Preview:

```text
/Users/tom/Fiji/cloud_extract_eval/angle_plugin_outputs_auto/NT_HTXR_K2_J1_preview_panel.jpg
```

### slide-HFLday8-IRmodel

Auto mode estimated two stain vectors with an original angle of about `18.99 degrees`.

The display angle was increased to `55.00 degrees`.

Preview:

```text
/Users/tom/Fiji/cloud_extract_eval/angle_plugin_outputs_auto/slide_HFLday8_IRmodel_R1_S10_11_3x_preview_panel.jpg
```

### CB1_PCOS_Carotis

Auto mode estimated two stain vectors with an original angle of about `20.11 degrees`.

The display angle was increased to `55.00 degrees`.

Preview:

```text
/Users/tom/Fiji/cloud_extract_eval/angle_plugin_outputs_auto/CB1_PCOS_Carotis_132_1_reexport_preview_panel.jpg
```

## Why This Is Different From the Previous Plugin

The previous H-DAB plugin is for stain separation and DAB-positive extraction.

This plugin is different:

```text
Previous plugin: separate and quantify DAB positivity.
This plugin: visually exaggerate weak color differences.
```

The previous plugin tries to stay biologically conservative.

This plugin intentionally changes the display color geometry, so it is useful for visual inspection but not for direct measurement.

## Advantages

### 1. It Works in OD Space

Weak stains mix more linearly in OD space than in raw RGB space.

This makes the separation more meaningful than a simple saturation, contrast, or hue adjustment.

### 2. It Can Estimate Vectors Automatically

The plugin does not require H-DAB. It can estimate two dominant stain directions from the image itself using an OD-PCA/Macenko-style method.

This makes it useful for exploratory inspection of non-standard stains.

### 3. It Explicitly Increases Vector Angle

If the two stain vectors are close, the plugin does not just increase saturation. It constructs new display vectors with a larger angular separation.

This directly targets the problem:

```text
two weak stains occupy nearby directions in color space
```

### 4. It Uses Local and Global Enhancement Together

Global percentile stretching helps broad weak stain.

Local contrast helps subtle regional differences.

The combination is more useful than either pure global contrast or pure local edge enhancement.

### 5. It Keeps Raw Density Outputs

Even though the visualization is artificial, the plugin also saves raw stain-density maps. This allows checking whether the enhanced color image is exaggerating a real separated channel or just visual noise.

## Limitations

### Not for Quantification

The angle-boosted color image changes color geometry and local contrast. It is not a valid quantitative readout.

For quantification, use the density maps or the H-DAB dominance plugin.

### Auto Vector Estimation Can Pick Artifacts

If tissue is sparse or artifacts dominate the OD cloud, auto mode may estimate stain directions from dust, folds, or background debris.

In that case, use preset or custom OD vectors.

### Close Vectors Are Still Ill-Conditioned

The plugin can make close colors easier to see, but it cannot create true information that was not captured by RGB.

If two stains are genuinely indistinguishable in RGB, a visualization plugin cannot fully solve that. Multispectral imaging or supervised annotation would be required.

### Display Colors Are Artificial

The output colors are chosen for separability, not for matching the original slide.

This is intentional.

## Recommended Use

Use this plugin as a viewing aid:

1. Start with Auto OD-PCA.
2. If the result is visually plausible, use the angle-boosted overlay to inspect weak staining.
3. Check `stain1_density` and `stain2_density` if something looks suspicious.
4. For known stains such as H-DAB, compare Auto with the H-DAB preset.
5. Do not use the angle-boosted overlay for area fraction or intensity measurement.

