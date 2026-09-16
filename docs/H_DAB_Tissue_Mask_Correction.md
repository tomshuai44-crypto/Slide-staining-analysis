> **Current installation:** The default is now local H+DAB transitions, `local-hdab-1.0`, k=8. The H-only method below remains available as `tissuemodel=h_only`. See [current plugin guide and Java validation](H_DAB_Local_Transition_Plugin.md). Earlier installation/benchmark descriptions below are historical.

# H-guided tissue mask correction — 16 September 2026

The installed H-DAB Dominance Extractor now selects a **candidate tissue area from spatially supported hematoxylin evidence**, then uses that area for its normalization and DAB measurements. The strict DAB-positive policy has not been redesigned. The selected area is an inspectable estimate; it is not a validated biological ground truth.

The canonical implementation is `plugin_src/H_DAB_Dominance_Extractor.java` with `plugin_src/H_DAB_Tissue_Selector.java`. Version `h-spatial-1.2` replaces the previous `OD norm >= 0.08` denominator. Build and installation are reproducible through `build_h_dab_tissue.py`; original installed JARs and the original main source are retained in `tissue_mask_validation/backup/`.

## What changed and why

The previous method asked only whether a pixel was darker than an almost-white threshold. It had no stain-specific tissue evidence, no spatial support, no ROI restriction, and no saved tissue mask. On the available NT JPEG it counted the entire 2,584 × 1,924 frame, including visible gray background, as tissue.

Your proposal contributes three useful principles: hematoxylin as tissue evidence, spatial exploration across the image, and H-density distributions rather than a single brightness decision. Four adjustments make it usable:

1. **The image center is only the first sample.** A center-out random walk could miss peripheral or disconnected tissue and could start in a lumen. Stratified random samples cover the whole frame; the final spatial scan visits every pixel inside the analysis domain.
2. **H-rich nuclei are evidence, not the whole denominator.** Supported H regions seed growth through nearby visible material. Pale inter-nuclear tissue can be included even when its H signal is below the seed threshold.
3. **A pixel is the smallest spatial unit, but a pixel alone is insufficient evidence.** A local H cluster and a sliding multi-pixel support window suppress isolated noise. There is no subpixel search or resampling of the input for segmentation.
4. **No reliable H means an explicit failed selection.** The algorithm does not substitute DAB positivity, darkness, or a plausible-looking silhouette when evidence is absent. A manual white reference helps when background cannot be estimated, but does not override the H requirement.

## Step-by-step mathematics

### 1. Analysis domain and white reference

Let Ω be the current image plane, restricted to the selected area ROI when one exists. Arbitrary area ROIs are respected pixel by pixel. Line and point ROIs are rejected. A stack is not processed in full; the current slice and image calibration are recorded in CSV. Spatial parameters remain **pixels**, regardless of the calibration unit.

Use 16 × 16 blocks only for sampling. Visit the central block first; divide the block lattice into at most 64 × 64 strata and choose one block uniformly from each stratum using `java.util.Random(20260916)`. Take 32 pixel samples with replacement from each selected block, discarding samples outside Ω. There are at most 4,097 sampled blocks and 131,104 samples. Sampling estimates distributions; it does not exclude unsampled tissue from the later full scan. Very small or sparse ROIs can have inadequate samples and should use an independently measured reference or be enlarged for review.

For each sampled block, estimate mean RGB, mean brightness, and brightness standard deviation. Among the brightest 10% of sampled blocks, use those in the lowest quartile of brightness standard deviation. The componentwise median of their mean RGB is the estimated white/background reference W. It need not be RGB 255: NT has a tinted, dim background.

Automatic reference estimation is marked `AUTO_UNRELIABLE` and fails if the selected flatness limit exceeds 8 RGB levels or any reference channel is below 128. Passing these heuristic checks means `AUTO_ESTIMATED_REVIEW`, not proof that a true blank background was observed. A uniform stained image can supply no valid blank reference; manual reference mode is available.

For tissue selection only, define corrected optical density:

\[
o_c(p)=-\ln\frac{I_c(p)+1}{W_c+1},\qquad
v(p)=\sqrt{\sum_c\max(0,o_c(p))^2}.
\]

Retain signed OD values for the chromatic calculation. Clip negative values only for visibility magnitude v. Do not clip them before computing chromatic differences: doing so would break neutral-offset cancellation.

### 2. Isolate chromatic H evidence from neutral darkness

Use normalized fixed vectors

\[
h=\operatorname{unit}(0.650,0.704,0.286),\quad
d=\operatorname{unit}(0.268,0.570,0.776).
\]

For the tissue selector, use the model

\[
o\approx a_Hh+a_Dd+q(1,1,1).
\]

Subtract blue OD from red and green OD, eliminating neutral attenuation q:

\[
\begin{bmatrix}o_R-o_B\\o_G-o_B\end{bmatrix}
=\begin{bmatrix}h_R-h_B&d_R-d_B\\h_G-h_B&d_G-d_B\end{bmatrix}
\begin{bmatrix}a_H\\a_D\end{bmatrix}.
\]

Only the first row of the inverse is needed. Numerically,

\[
a_H\approx-1.50(o_R-o_B)+3.70(o_G-o_B).
\]

The exact coefficients are computed from the vectors in Java; these displayed coefficients are rounded. A purely neutral OD increase gives zero H evidence. Under this model pure DAB also gives zero H evidence, and adding DAB to an H mixture does not itself create an H seed. The D vector is used to separate H, **never as a positive tissue criterion**. Colored debris or stain-vector mismatch can still mimic H; an algebraic separation cannot determine biological identity.

### 3. Estimate an H-density threshold

Exclude pixels whose minimum RGB channel is at or below 10, a configurable dark/saturated-pixel safeguard. It is not a complete artifact detector.

From the Monte Carlo samples, collect A = H values at or above the user minimum (default 0.08). From low-visibility samples with v < 0.08, estimate a background H median m and robust spread

\[
\sigma_H=1.4826\operatorname{median}|a_H-m|.
\]

Choose

\[
T_H=\max\bigl(0.08,\;0.5Q_{0.90}(A),\;m+4\sigma_H\bigr).
\]

The 0.08 floor, half-upper-decile term, and four-MAD noise guard are explicit engineering defaults, not externally validated universal biological constants. An empty sample set contributes zero to its estimate. A full deterministic scan still checks every pixel against the resulting threshold. On these saved JPEGs, the final thresholds are about 0.15885 (NT), 0.11573 (CB1), and 0.14991 (HFL).

### 4. Require clustered, locally supported H

Define visible material C = {p ∈ Ω: v(p) ≥ 0.08 and not dark/saturated}. A strong-H pixel belongs to C and satisfies a_H ≥ T_H.

A supported-H pixel must have at least **three strong-H pixels, including itself, in its 3 × 3 neighborhood**. Then evaluate supported-H density in a sliding 17 × 17 neighborhood centered at every pixel (radius 8). A supported-H pixel becomes a seed only if that window contains at least

\[
\max\bigl(4,\lceil0.02|\text{window}\cap\Omega|\rceil\bigr)
\]

supported-H pixels. A full window therefore requires six pixels. Edge and ROI windows use their actual domain area. A single colored pixel or separated single-pixel speckle cannot seed tissue. A cluster of colored artifact pixels can, which is why the seed image is exported for review.

### 5. Include pale tissue by bounded spatial growth

Apply a square morphological closing of radius 2 pixels to C to bridge very small gaps. The implementation preserves existing candidate pixels and truncates neighborhoods at image edges; it does not erode edge tissue. The closed mask is clipped to Ω and excludes dark/saturated pixels again. Clear spaces wider than this small gap scale remain barriers; there is no general hole filling or convex hull.

Growth may pass through a closed candidate pixel if at least 20% of its sliding 17 × 17 domain window was visible in the original candidate mask. Seed pixels are retained even if this density condition is not met.

Run a multi-source shortest-path expansion from **all** seed pixels on this allowed mask. Cardinal steps cost 3 and diagonal steps cost 4. Select pixels at cumulative cost ≤ 3R, with default R = **96 pixels**. This chamfer metric approximates Euclidean distance with octagonal rather than square reach contours. It is not an exact physical-distance metric, particularly with anisotropic pixels. Bucketed distance propagation keeps this bounded computation efficient.

Thus the final tissue mask is

\[
T=\{p\in\text{allowed material}:d_{3,4}(p,\text{H seeds})\le 3R\}.
\]

Pale material need not pass the H threshold. It must be locally visible and reachable from H evidence. Background and lumens block growth. Completely unstained, optically indistinguishable tissue, very sparse nuclei, long acellular stroma, folds, and colored debris remain intrinsic limitations.

The initial uninstalled prototype used hard tile gates and a square reach metric. Close inspection showed both split evidence at tile boundaries and artificial square cutoff contours. The final version uses overlapping neighborhoods and pixel paths. The default reach was increased from the prototype's nominal 64-pixel scale to 96 pixels after inspecting sparse-H regions in NT; this is a review-informed engineering choice, not a ground-truth optimization. The validation folder includes 64- and 128-pixel sensitivity runs and enlarged prototype/final comparisons. Remaining small omissions are shown rather than smoothed away.

### 6. Apply the existing DAB measurement policy within T

The original fixed three-channel deconvolution remains unchanged, including its residual vector and original white reference:

\[
OD_c^{legacy}=-\ln((I_c+1)/256).
\]

Raw H, DAB, and residual concentration images remain full-frame diagnostic outputs. Tissue selection does not change their pixel values. For downstream statistics, channel p99 scales, automatic threshold selection, DAB-positive classification, and the area denominator all use T.

\[
P=T\cap\{D_n\ge D_{min},\ D_n-\alpha H_n\ge\delta,
\ D_n\ge1.15H_n,\ R_n\le R_{max}\}.
\]

The automatic p85 DAB floor and optimizer are unchanged. Positive-pixel mean and integrated DAB still use the original raw DAB concentrations over P. Changing T can change the p99 scales, optimized thresholds, and P, so this change is **more than replacing a denominator in an old CSV**. Re-run the entire analysis when comparing corrected results. `concentration_reference=LEGACY_255_UNCORRECTED` makes the intentional separation between tissue-reference correction and legacy intensity measurement explicit.

## Parameters and outputs

The dialog exposes corrected-OD minimum, H minimum, support radius, maximum reach, closing radius, and automatic/manual reference settings. Advanced macro keys follow ImageJ's **lowercase** convention:

| Macro key | Default | Meaning |
|---|---:|---|
| `tissue` or `tissueminod` | 0.08 | Minimum corrected visibility OD; must be positive |
| `tissuehmin` | 0.08 | H seed threshold floor |
| `tissuesupportradius` | 8 | Sliding window radius: side = 2r+1 |
| `tissuereachpixels` | 96 | Maximum chamfer-geodesic reach in pixel units |
| `tissuecloseradius` | 2 | Small-gap closing radius |
| `tissueminseedpixels` | 4 | Absolute minimum supported H pixels in window |
| `tissueseedfraction` | 0.02 | Supported H fraction required for a seed |
| `tissuevisiblefraction` | 0.20 | Visible fraction required for growth |
| `tissuedarkmax` | 10 | Reject any pixel with a channel ≤ this value |
| `tissueblocksize` | 16 | Sampling block size only; no final tile decisions |
| `tissueseed` | 20260916 | Signed 64-bit PRNG seed |
| `tissueautowhite` | true | Estimate tissue reference from this domain |
| `tissuewhiter`, `tissuewhiteg`, `tissuewhiteb` | 255 | Manual reference RGB, used when automatic is false |

Example:

```java
run("H-DAB Dominance Extractor", "tissue=0.08 tissuehmin=0.08 tissuesupportradius=8 tissuereachpixels=96 tissuecloseradius=2 tissueautowhite=true tissueseed=20260916 show=true output=[/absolute/output/folder] prefix=sample");
```

New files are `*_tissue_mask.tif` (255 selected, 0 excluded), `*_tissue_overlay.tif` (green = selected), and `*_tissue_h_evidence.tif` (actual seed pixels). The mask and overlay appear with normal output display. All three are saved alongside the previous outputs. CSV includes algorithm version, selection/reference status, current plane/ROI scope, domain pixel count, calibration, reference RGB, estimated H threshold/noise, all effective spatial settings, sampling settings, seed counts, and downstream normalization/optimization provenance.

`REVIEW_REQUIRED` means selection was produced but is not biologically validated. Failures distinguish empty ROI, unreliable automatic reference, no reliable H, and fewer than ten selected pixels. Failed selections export an empty tissue mask, zero positive pixels, and NaN area fraction, concentration summary values, and normalization scales. There is no forced fallback. A valid tissue selection with zero DAB positives retains the existing zero-positive summary convention.

## Validation and practical use

Read `tissue_mask_validation/validation_report.md` for final counts, test evidence, sensitivity, and links to the image panels. These tests exercise the actual built and installed Java implementation. Historical Python benchmark metrics used a different tissue rule and **are not validation of this correction or the old installed Java tissue mask**.

1. Restart Fiji to load the replaced JAR, open an RGB H-DAB image, and optionally draw an area ROI. Run Plugins → Analyze → H-DAB Dominance Extractor with output saving enabled.
2. Inspect the original, seed map, tissue mask, and overlay at native scale. Check outer background, lumens, pale stroma, folds, debris, and disconnected fragments before interpreting the CSV.
3. Keep pixel scales consistent across images at the same acquisition resolution. When resolution changes, use the recorded calibration to choose a comparable physical support/reach; automatic scale conversion is not implemented.
4. Use a measured blank-field reference when automatic reference estimation is suspect. A manual area ROI narrows the domain but still requires H evidence. Very weak or absent H requires an explicitly different validated tissue-selection approach.
5. Establish manual reference masks across representative tissues, stain batches, and magnifications before calling this a biologically accurate area measurement. Compare mask overlap and area bias, and tune parameters on a development set rather than on the final evaluation set.

Alternatives include a manually traced tissue mask, a calibrated blank-field plus supervised tissue classifier, or a multiscale texture/stain segmentation model trained on annotated examples. A nuclear convex hull is simpler but can fill vascular lumens or bridge separate fragments. Unbounded region growth may recover long acellular stroma but risks leaking into dim background. A center-out Monte Carlo walk is useful for exploration but is a poor sole mechanism for complete area coverage. None of these alternatives is implemented silently in this release.

## Scientific basis and limits of attribution

Optical-density stain unmixing is supported by [Ruifrok and Johnston, 2001](https://pubmed.ncbi.nlm.nih.gov/11531144/). The background/reference correction, fixed vector assumptions, and limitations of quantitative interpretation are discussed in [Landini, Martinelli and Piccinini, 2021](https://academic.oup.com/bioinformatics/article/37/10/1485/5913390). Hematoxylin-specific evidence and object-size filtering have also been used in [a primary nuclear-segmentation study](https://www.nature.com/articles/srep00503).

Those papers support the underlying concepts, not this exact tissue algorithm, its numerical defaults, or its accuracy. The neutral-invariant H projection, sampling scheme, spatial rules, failure handling, and parameter choices above are the implementation developed here. DAB intensity remains an image-derived stain measurement; this change does not establish antigen concentration or resolve the separate strict-positive detection policy.
