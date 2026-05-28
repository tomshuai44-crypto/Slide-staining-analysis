# Fiji Stain Analysis Plugins

This folder contains two experimental Fiji/ImageJ plugins for brightfield stain analysis and visualization.

## Plugins

### H-DAB Dominance Extractor

`H_DAB_Dominance_Extractor` is intended for H-DAB color deconvolution and conservative DAB-positive extraction.

It outputs:

- hematoxylin, DAB, and residual color deconvolution images
- hematoxylin, DAB, and residual concentration maps
- DAB dominance map
- DAB-positive mask
- measurement CSV

Compiled plugin:

```text
fiji-plugins/H_DAB_Dominance_Extractor.jar
```

Source:

```text
src/main/java/H_DAB_Dominance_Extractor.java
```

Detailed documentation:

```text
docs/H_DAB_Dominance_Engineering_Document.md
docs/H_DAB_Dominance_Plugin_README.md
```

### Stain Angle Contrast Booster

`Stain_Angle_Contrast_Booster` is a visualization-only tool. It estimates or accepts two stain vectors in optical-density space, separates the two stain-density maps, and re-renders them with a larger display angle so weak color differences become easier to see.

It is not intended for direct quantification.

Compiled plugin:

```text
fiji-plugins/Stain_Angle_Contrast_Booster.jar
```

Source:

```text
src/main/java/Stain_Angle_Contrast_Booster.java
```

Detailed documentation:

```text
docs/Stain_Angle_Contrast_Booster_Document.md
```

## Fiji Installation

Copy both JAR files into Fiji's `plugins` directory:

```text
Fiji.app/plugins/
```

After restarting Fiji, the plugins appear under:

```text
Plugins > Analyze > H-DAB Dominance Extractor
Plugins > Analyze > Stain Angle Contrast Booster
```

The Fiji menu config used during local builds is included at:

```text
fiji-config/plugins.config
```

## Build

From a Fiji installation root:

```bash
JDK="$PWD/java/macos-arm64/zulu21.42.19-ca-jdk21.0.7-macosx_aarch64/zulu-21.jdk/Contents/Home"

"$JDK/bin/javac" --release 8 -cp jars/ij-1.54p.jar \
  -d build/classes \
  src/main/java/H_DAB_Dominance_Extractor.java \
  src/main/java/Stain_Angle_Contrast_Booster.java

cp fiji-config/plugins.config build/classes/plugins.config

"$JDK/bin/jar" cf fiji-plugins/Stain_Analysis_Plugins.jar \
  -C build/classes .
```

The checked-in `fiji-plugins/*.jar` files are the already built plugin artifacts from the local Fiji environment.

## Status

These plugins are experimental research utilities.

- Use `H-DAB Dominance Extractor` for H-DAB deconvolution and conservative DAB-positive extraction.
- Use `Stain Angle Contrast Booster` only as a visual aid for weakly separated stain colors.

