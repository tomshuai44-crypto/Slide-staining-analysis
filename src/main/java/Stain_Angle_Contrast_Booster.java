import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class Stain_Angle_Contrast_Booster implements PlugIn {
    private static final double[] HEMA_REF = normalize(new double[] {0.650, 0.704, 0.286});
    private static final double[] DAB_REF = normalize(new double[] {0.268, 0.570, 0.776});
    private static final double[] EOSIN_REF = normalize(new double[] {0.072, 0.990, 0.105});

    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("Stain Angle Contrast Booster", "Open an RGB brightfield image first.");
            return;
        }
        if (imp.getType() != ImagePlus.COLOR_RGB) {
            IJ.error("Stain Angle Contrast Booster", "This plugin expects an RGB image.");
            return;
        }

        Params params = Params.fromOptions(arg);
        if (!params.headless && !params.fromMacro && !showDialog(params)) {
            return;
        }

        Result result = analyze(imp, params);
        if (result == null) {
            IJ.error("Stain Angle Contrast Booster", "Could not analyze this image.");
            return;
        }

        if (params.showOutputs) {
            result.angleBoostedImage.show();
            result.stain1EnhancedImage.show();
            result.stain2EnhancedImage.show();
            result.stain1DensityImage.show();
            result.stain2DensityImage.show();
            result.measurements.show("Stain Angle Contrast Booster");
        }

        if (params.outputDir.length() > 0) {
            try {
                result.save(params.outputDir, params.prefix.length() > 0 ? params.prefix : sanitizeTitle(imp.getTitle()));
            } catch (IOException e) {
                IJ.error("Stain Angle Contrast Booster", "Could not save outputs:\n" + e.getMessage());
            }
        }
    }

    public static Result analyze(ImagePlus imp, Params params) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        int n = width * height;
        int[] pixels = (int[]) imp.getProcessor().convertToRGB().getPixels();

        double[] white = estimateWhite(pixels, params.whitePercentile);
        boolean[] tissue = new boolean[n];
        int tissueCount = markTissue(pixels, white, tissue, params.tissueMinOD);
        if (tissueCount < 50) return null;

        StainPair stains = chooseStains(pixels, white, tissue, tissueCount, params);
        double[] residual = residualVector(stains.v1, stains.v2);
        double[][] inv = inverse(new double[][] {
                {stains.v1[0], stains.v2[0], residual[0]},
                {stains.v1[1], stains.v2[1], residual[1]},
                {stains.v1[2], stains.v2[2], residual[2]},
        });

        double[] c1 = new double[n];
        double[] c2 = new double[n];
        for (int i = 0; i < n; i++) {
            int rgb = pixels[i];
            double[] od = rgbToOD(rgb, white);
            c1[i] = Math.max(0.0, inv[0][0] * od[0] + inv[0][1] * od[1] + inv[0][2] * od[2]);
            c2[i] = Math.max(0.0, inv[1][0] * od[0] + inv[1][1] * od[1] + inv[1][2] * od[2]);
        }

        double c1Scale = Math.max(maskedPercentile(c1, tissue, params.scalePercentile), 1e-9);
        double c2Scale = Math.max(maskedPercentile(c2, tissue, params.scalePercentile), 1e-9);
        double[] displayVectors = expandedAngleVectors(stains.v1, stains.v2, params);
        double[] u1 = new double[] {displayVectors[0], displayVectors[1], displayVectors[2]};
        double[] u2 = new double[] {displayVectors[3], displayVectors[4], displayVectors[5]};
        double displayAngle = vectorAngleDegrees(u1, u2);

        float[] vis1 = visualDensity(width, height, c1, c1Scale, params);
        float[] vis2 = visualDensity(width, height, c2, c2Scale, params);

        ColorProcessor overlay = angleBoostedProcessor(width, height, vis1, vis2, tissue, u1, u2, params);
        ColorProcessor stain1 = singleStainProcessor(width, height, vis1, tissue, u1, params);
        ColorProcessor stain2 = singleStainProcessor(width, height, vis2, tissue, u2, params);
        FloatProcessor density1 = floatProcessor(width, height, c1);
        FloatProcessor density2 = floatProcessor(width, height, c2);

        String base = sanitizeTitle(imp.getTitle());
        ResultsTable rt = new ResultsTable();
        rt.incrementCounter();
        rt.addValue("image", imp.getTitle());
        rt.addValue("mode", stains.mode);
        rt.addValue("visualization_only", 1);
        rt.addValue("width", width);
        rt.addValue("height", height);
        rt.addValue("tissue_pixels", tissueCount);
        rt.addValue("white_R", white[0]);
        rt.addValue("white_G", white[1]);
        rt.addValue("white_B", white[2]);
        rt.addValue("stain1_R_OD", stains.v1[0]);
        rt.addValue("stain1_G_OD", stains.v1[1]);
        rt.addValue("stain1_B_OD", stains.v1[2]);
        rt.addValue("stain2_R_OD", stains.v2[0]);
        rt.addValue("stain2_G_OD", stains.v2[1]);
        rt.addValue("stain2_B_OD", stains.v2[2]);
        rt.addValue("original_angle_degrees", vectorAngleDegrees(stains.v1, stains.v2));
        rt.addValue("display_stain1_R_OD", u1[0]);
        rt.addValue("display_stain1_G_OD", u1[1]);
        rt.addValue("display_stain1_B_OD", u1[2]);
        rt.addValue("display_stain2_R_OD", u2[0]);
        rt.addValue("display_stain2_G_OD", u2[1]);
        rt.addValue("display_stain2_B_OD", u2[2]);
        rt.addValue("display_angle_degrees", displayAngle);
        rt.addValue("angle_gain", params.angleGain);
        rt.addValue("minimum_display_angle_degrees", params.minDisplayAngleDegrees);
        rt.addValue("maximum_display_angle_degrees", params.maxDisplayAngleDegrees);
        rt.addValue("stain1_density_p" + params.scalePercentile, c1Scale);
        rt.addValue("stain2_density_p" + params.scalePercentile, c2Scale);
        rt.addValue("local_radius_px", params.localRadius);
        rt.addValue("local_weight", params.localWeight);
        rt.addValue("chroma_gain", params.chromaGain);
        rt.addValue("density_gamma", params.densityGamma);

        return new Result(
                new ImagePlus(base + " - angle boosted visual overlay", overlay),
                new ImagePlus(base + " - stain 1 enhanced visual", stain1),
                new ImagePlus(base + " - stain 2 enhanced visual", stain2),
                new ImagePlus(base + " - stain 1 density", density1),
                new ImagePlus(base + " - stain 2 density", density2),
                rt);
    }

    private static boolean showDialog(Params params) {
        GenericDialog gd = new GenericDialog("Stain Angle Contrast Booster");
        String[] modes = new String[] {"Auto OD-PCA", "H-DAB", "H&E", "Custom OD vectors"};
        gd.addChoice("Stain vector source", modes, params.mode);
        gd.addNumericField("Tissue minimum OD", params.tissueMinOD, 3);
        gd.addNumericField("White reference percentile", params.whitePercentile, 2);
        gd.addNumericField("Auto vector extreme percentile", params.extremePercentile, 2);
        gd.addNumericField("Angle gain", params.angleGain, 2);
        gd.addNumericField("Minimum display angle (degrees)", params.minDisplayAngleDegrees, 1);
        gd.addNumericField("Maximum display angle (degrees)", params.maxDisplayAngleDegrees, 1);
        gd.addNumericField("Local contrast radius (pixels)", params.localRadius, 1);
        gd.addNumericField("Local contrast weight", params.localWeight, 2);
        gd.addNumericField("Color/chroma gain", params.chromaGain, 2);
        gd.addNumericField("Density gamma", params.densityGamma, 2);
        gd.addNumericField("Density scale percentile", params.scalePercentile, 1);
        gd.addNumericField("Display OD strength", params.displayStrength, 2);
        gd.addStringField("Custom stain 1 OD vector", params.customVector1, 24);
        gd.addStringField("Custom stain 2 OD vector", params.customVector2, 24);
        gd.addCheckbox("Show output images", params.showOutputs);
        gd.addStringField("Output directory (optional)", params.outputDir, 32);
        gd.addStringField("Output prefix (optional)", params.prefix, 24);
        gd.showDialog();
        if (gd.wasCanceled()) return false;

        params.mode = gd.getNextChoice();
        params.tissueMinOD = gd.getNextNumber();
        params.whitePercentile = gd.getNextNumber();
        params.extremePercentile = gd.getNextNumber();
        params.angleGain = gd.getNextNumber();
        params.minDisplayAngleDegrees = gd.getNextNumber();
        params.maxDisplayAngleDegrees = gd.getNextNumber();
        params.localRadius = gd.getNextNumber();
        params.localWeight = gd.getNextNumber();
        params.chromaGain = gd.getNextNumber();
        params.densityGamma = gd.getNextNumber();
        params.scalePercentile = gd.getNextNumber();
        params.displayStrength = gd.getNextNumber();
        params.customVector1 = gd.getNextString();
        params.customVector2 = gd.getNextString();
        params.showOutputs = gd.getNextBoolean();
        params.outputDir = gd.getNextString();
        params.prefix = gd.getNextString();
        return true;
    }

    private static StainPair chooseStains(int[] pixels, double[] white, boolean[] tissue, int tissueCount, Params params) {
        String key = normalizeMode(params.mode);
        if (key.equals("hdab")) {
            return new StainPair(HEMA_REF, DAB_REF, "H-DAB preset");
        }
        if (key.equals("he")) {
            return new StainPair(HEMA_REF, EOSIN_REF, "H&E preset");
        }
        if (key.equals("custom")) {
            return new StainPair(parseVector(params.customVector1, HEMA_REF), parseVector(params.customVector2, DAB_REF), "custom OD vectors");
        }
        return estimateMacenkoPair(pixels, white, tissue, tissueCount, params);
    }

    private static StainPair estimateMacenkoPair(int[] pixels, double[] white, boolean[] tissue, int tissueCount, Params params) {
        double[][] cov = new double[3][3];
        for (int i = 0; i < pixels.length; i++) {
            if (!tissue[i]) continue;
            double[] od = rgbToOD(pixels[i], white);
            for (int r = 0; r < 3; r++) {
                for (int c = r; c < 3; c++) cov[r][c] += od[r] * od[c];
            }
        }
        for (int r = 0; r < 3; r++) {
            for (int c = r; c < 3; c++) {
                cov[r][c] /= Math.max(1, tissueCount);
                cov[c][r] = cov[r][c];
            }
        }
        Eigen eig = eigenSymmetric3x3(cov);
        double[] e1 = eig.vector(0);
        double[] e2 = eig.vector(1);
        if (sum(e1) < 0) e1 = multiply(e1, -1.0);
        if (sum(e2) < 0) e2 = multiply(e2, -1.0);

        int maxSamples = 300000;
        int stride = Math.max(1, tissueCount / maxSamples);
        double[] angles = new double[Math.min(tissueCount, maxSamples + 1)];
        int k = 0;
        int seen = 0;
        for (int i = 0; i < pixels.length; i++) {
            if (!tissue[i]) continue;
            if ((seen++ % stride) != 0) continue;
            double[] od = rgbToOD(pixels[i], white);
            double x = dot(od, e1);
            double y = dot(od, e2);
            if (Math.sqrt(x * x + y * y) <= 1e-9) continue;
            if (k >= angles.length) break;
            angles[k++] = Math.atan2(y, x);
        }
        if (k < 50) return new StainPair(HEMA_REF, DAB_REF, "auto fallback H-DAB");
        angles = Arrays.copyOf(angles, k);
        Arrays.sort(angles);
        double low = percentileSorted(angles, params.extremePercentile);
        double high = percentileSorted(angles, 100.0 - params.extremePercentile);
        double[] v1 = positiveNormalize(add(multiply(e1, Math.cos(low)), multiply(e2, Math.sin(low))));
        double[] v2 = positiveNormalize(add(multiply(e1, Math.cos(high)), multiply(e2, Math.sin(high))));
        if (vectorAngleDegrees(v1, v2) < 3.0 || !allFinite(v1) || !allFinite(v2)) {
            return new StainPair(HEMA_REF, DAB_REF, "auto fallback H-DAB");
        }
        return new StainPair(v1, v2, "auto OD-PCA/Macenko-style");
    }

    private static double[] estimateWhite(int[] pixels, double pct) {
        int[][] hist = new int[3][256];
        for (int rgb : pixels) {
            hist[0][(rgb >> 16) & 0xff]++;
            hist[1][(rgb >> 8) & 0xff]++;
            hist[2][rgb & 0xff]++;
        }
        double[] white = new double[3];
        int target = (int) Math.ceil((pct / 100.0) * pixels.length);
        for (int c = 0; c < 3; c++) {
            int cumulative = 0;
            int value = 255;
            for (int i = 0; i < 256; i++) {
                cumulative += hist[c][i];
                if (cumulative >= target) {
                    value = i;
                    break;
                }
            }
            white[c] = Math.max(180.0, value);
        }
        return white;
    }

    private static int markTissue(int[] pixels, double[] white, boolean[] tissue, double tissueMinOD) {
        int count = 0;
        for (int i = 0; i < pixels.length; i++) {
            double[] od = rgbToOD(pixels[i], white);
            double norm = Math.sqrt(od[0] * od[0] + od[1] * od[1] + od[2] * od[2]);
            tissue[i] = norm >= tissueMinOD;
            if (tissue[i]) count++;
        }
        return count;
    }

    private static double[] rgbToOD(int rgb, double[] white) {
        double r = (rgb >> 16) & 0xff;
        double g = (rgb >> 8) & 0xff;
        double b = rgb & 0xff;
        return new double[] {
                Math.max(0.0, -Math.log((r + 1.0) / (white[0] + 1.0))),
                Math.max(0.0, -Math.log((g + 1.0) / (white[1] + 1.0))),
                Math.max(0.0, -Math.log((b + 1.0) / (white[2] + 1.0))),
        };
    }

    private static float[] visualDensity(int width, int height, double[] density, double scale, Params params) {
        float[] base = new float[density.length];
        float[] sq = new float[density.length];
        for (int i = 0; i < density.length; i++) {
            base[i] = (float) density[i];
            sq[i] = (float) (density[i] * density[i]);
        }
        FloatProcessor mean = new FloatProcessor(width, height, base.clone());
        mean.blurGaussian(params.localRadius);
        FloatProcessor meanSq = new FloatProcessor(width, height, sq);
        meanSq.blurGaussian(params.localRadius);
        float[] meanPix = (float[]) mean.getPixels();
        float[] meanSqPix = (float[]) meanSq.getPixels();
        float[] out = new float[density.length];
        double denomFloor = scale * 0.04 + 1e-9;
        for (int i = 0; i < density.length; i++) {
            double localMean = meanPix[i];
            double variance = Math.max(0.0, meanSqPix[i] - localMean * localMean);
            double localStd = Math.sqrt(variance);
            double global = density[i] / scale;
            double local = Math.max(0.0, (density[i] - localMean) / (localStd + denomFloor));
            double boosted = Math.max(0.0, global + params.localWeight * 0.35 * local);
            double stretched = asinh(params.chromaGain * boosted) / asinh(params.chromaGain * 1.5);
            out[i] = (float) Math.pow(clamp(stretched, 0.0, 1.0), params.densityGamma);
        }
        return out;
    }

    private static ColorProcessor angleBoostedProcessor(int width, int height, float[] vis1, float[] vis2,
                                                        boolean[] tissue, double[] v1, double[] v2, Params params) {
        int[] out = new int[vis1.length];
        for (int i = 0; i < out.length; i++) {
            if (!tissue[i]) {
                out[i] = 0xffffff;
                continue;
            }
            double a = params.displayStrength * vis1[i];
            double b = params.displayStrength * vis2[i];
            int r = clamp255(255.0 * Math.exp(-(a * v1[0] + b * v2[0])));
            int g = clamp255(255.0 * Math.exp(-(a * v1[1] + b * v2[1])));
            int bl = clamp255(255.0 * Math.exp(-(a * v1[2] + b * v2[2])));
            out[i] = (r << 16) | (g << 8) | bl;
        }
        return new ColorProcessor(width, height, out);
    }

    private static ColorProcessor singleStainProcessor(int width, int height, float[] vis, boolean[] tissue,
                                                       double[] vector, Params params) {
        int[] out = new int[vis.length];
        for (int i = 0; i < out.length; i++) {
            if (!tissue[i]) {
                out[i] = 0xffffff;
                continue;
            }
            double a = params.displayStrength * vis[i];
            int r = clamp255(255.0 * Math.exp(-a * vector[0]));
            int g = clamp255(255.0 * Math.exp(-a * vector[1]));
            int b = clamp255(255.0 * Math.exp(-a * vector[2]));
            out[i] = (r << 16) | (g << 8) | b;
        }
        return new ColorProcessor(width, height, out);
    }

    private static double[] expandedAngleVectors(double[] v1, double[] v2, Params params) {
        double original = vectorAngleDegrees(v1, v2);
        double target = Math.max(original * params.angleGain, params.minDisplayAngleDegrees);
        target = clamp(target, original, params.maxDisplayAngleDegrees);
        double half = Math.toRadians(target / 2.0);
        double[] axis = normalize(add(v1, v2));
        double[] delta = subtract(v1, v2);
        delta = subtract(delta, multiply(axis, dot(delta, axis)));
        if (norm(delta) < 1e-9) delta = orthogonalUnit(axis);
        delta = normalize(delta);
        double[] u1 = positiveNormalize(add(multiply(axis, Math.cos(half)), multiply(delta, Math.sin(half))));
        double[] u2 = positiveNormalize(add(multiply(axis, Math.cos(half)), multiply(delta, -Math.sin(half))));
        return new double[] {u1[0], u1[1], u1[2], u2[0], u2[1], u2[2]};
    }

    private static FloatProcessor floatProcessor(int width, int height, double[] values) {
        float[] out = new float[values.length];
        for (int i = 0; i < values.length; i++) out[i] = (float) values[i];
        FloatProcessor fp = new FloatProcessor(width, height, out);
        fp.resetMinAndMax();
        return fp;
    }

    private static double maskedPercentile(double[] values, boolean[] mask, double pct) {
        double max = 0.0;
        int count = 0;
        for (int i = 0; i < values.length; i++) {
            if (!mask[i] || !Double.isFinite(values[i])) continue;
            max = Math.max(max, values[i]);
            count++;
        }
        if (count == 0 || max <= 0.0) return 0.0;
        int bins = 4096;
        int[] hist = new int[bins];
        for (int i = 0; i < values.length; i++) {
            if (!mask[i] || !Double.isFinite(values[i])) continue;
            int bin = (int) Math.floor(clamp(values[i] / max, 0.0, 1.0) * (bins - 1));
            hist[bin]++;
        }
        int target = (int) Math.ceil((pct / 100.0) * count);
        int cumulative = 0;
        for (int i = 0; i < bins; i++) {
            cumulative += hist[i];
            if (cumulative >= target) return max * i / (double) (bins - 1);
        }
        return max;
    }

    private static double percentileSorted(double[] sorted, double pct) {
        if (sorted.length == 0) return 0.0;
        double rank = (pct / 100.0) * (sorted.length - 1);
        int lo = (int) Math.floor(rank);
        int hi = (int) Math.ceil(rank);
        if (lo == hi) return sorted[lo];
        double frac = rank - lo;
        return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
    }

    private static double[] residualVector(double[] v1, double[] v2) {
        double[] cross = cross(v1, v2);
        if (norm(cross) < 1e-9) cross = orthogonalUnit(v1);
        return normalize(cross);
    }

    private static double[][] inverse(double[][] m) {
        double a = m[0][0], b = m[0][1], c = m[0][2];
        double d = m[1][0], e = m[1][1], f = m[1][2];
        double g = m[2][0], h = m[2][1], i = m[2][2];
        double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
        if (Math.abs(det) < 1e-12) det = det >= 0 ? 1e-12 : -1e-12;
        return new double[][] {
                {(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det},
                {(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det},
                {(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det},
        };
    }

    private static Eigen eigenSymmetric3x3(double[][] input) {
        double[][] a = new double[3][3];
        double[][] v = new double[][] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) a[r][c] = input[r][c];
        }
        for (int iter = 0; iter < 60; iter++) {
            int p = 0, q = 1;
            double max = Math.abs(a[p][q]);
            if (Math.abs(a[0][2]) > max) { p = 0; q = 2; max = Math.abs(a[0][2]); }
            if (Math.abs(a[1][2]) > max) { p = 1; q = 2; max = Math.abs(a[1][2]); }
            if (max < 1e-12) break;
            double theta = 0.5 * Math.atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
            double cos = Math.cos(theta);
            double sin = Math.sin(theta);
            for (int k = 0; k < 3; k++) {
                double apk = a[p][k], aqk = a[q][k];
                a[p][k] = cos * apk - sin * aqk;
                a[q][k] = sin * apk + cos * aqk;
            }
            for (int k = 0; k < 3; k++) {
                double akp = a[k][p], akq = a[k][q];
                a[k][p] = cos * akp - sin * akq;
                a[k][q] = sin * akp + cos * akq;
            }
            for (int k = 0; k < 3; k++) {
                double vkp = v[k][p], vkq = v[k][q];
                v[k][p] = cos * vkp - sin * vkq;
                v[k][q] = sin * vkp + cos * vkq;
            }
        }
        double[] vals = new double[] {a[0][0], a[1][1], a[2][2]};
        double[][] vecs = new double[][] {
                {v[0][0], v[1][0], v[2][0]},
                {v[0][1], v[1][1], v[2][1]},
                {v[0][2], v[1][2], v[2][2]},
        };
        for (int i = 0; i < 2; i++) {
            for (int j = i + 1; j < 3; j++) {
                if (vals[j] > vals[i]) {
                    double tv = vals[i]; vals[i] = vals[j]; vals[j] = tv;
                    double[] tmp = vecs[i]; vecs[i] = vecs[j]; vecs[j] = tmp;
                }
            }
        }
        return new Eigen(vals, vecs);
    }

    private static double[] parseVector(String text, double[] fallback) {
        if (text == null) return fallback;
        String[] parts = text.trim().split("[,;\\s]+");
        if (parts.length < 3) return fallback;
        try {
            return positiveNormalize(new double[] {
                    Double.parseDouble(parts[0]),
                    Double.parseDouble(parts[1]),
                    Double.parseDouble(parts[2])
            });
        } catch (NumberFormatException e) {
            return fallback;
        }
    }

    private static String normalizeMode(String mode) {
        if (mode == null) return "auto";
        String m = mode.toLowerCase().replaceAll("[^a-z0-9]+", "");
        if (m.contains("hdab")) return "hdab";
        if (m.equals("he") || m.contains("hepreset")) return "he";
        if (m.contains("custom")) return "custom";
        return "auto";
    }

    private static double[] normalize(double[] v) {
        double norm = norm(v);
        return new double[] {v[0] / Math.max(norm, 1e-12), v[1] / Math.max(norm, 1e-12), v[2] / Math.max(norm, 1e-12)};
    }

    private static double[] positiveNormalize(double[] v) {
        if (sum(v) < 0.0) v = multiply(v, -1.0);
        double[] out = new double[] {Math.max(1e-6, v[0]), Math.max(1e-6, v[1]), Math.max(1e-6, v[2])};
        return normalize(out);
    }

    private static double[] orthogonalUnit(double[] axis) {
        double[] ref = Math.abs(axis[0]) < 0.8 ? new double[] {1, 0, 0} : new double[] {0, 1, 0};
        double[] out = subtract(ref, multiply(axis, dot(ref, axis)));
        return normalize(out);
    }

    private static double norm(double[] v) {
        return Math.sqrt(dot(v, v));
    }

    private static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double sum(double[] a) {
        return a[0] + a[1] + a[2];
    }

    private static double[] add(double[] a, double[] b) {
        return new double[] {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }

    private static double[] subtract(double[] a, double[] b) {
        return new double[] {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    private static double[] multiply(double[] a, double k) {
        return new double[] {a[0] * k, a[1] * k, a[2] * k};
    }

    private static double[] cross(double[] a, double[] b) {
        return new double[] {
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
        };
    }

    private static double vectorAngleDegrees(double[] a, double[] b) {
        return Math.toDegrees(Math.acos(clamp(dot(a, b) / Math.max(norm(a) * norm(b), 1e-12), -1.0, 1.0)));
    }

    private static boolean allFinite(double[] v) {
        return Double.isFinite(v[0]) && Double.isFinite(v[1]) && Double.isFinite(v[2]);
    }

    private static double clamp(double x, double lo, double hi) {
        return Math.max(lo, Math.min(hi, x));
    }

    private static double asinh(double x) {
        return Math.log(x + Math.sqrt(x * x + 1.0));
    }

    private static int clamp255(double x) {
        return (int) Math.max(0, Math.min(255, Math.round(x)));
    }

    private static String sanitizeTitle(String title) {
        String s = title == null ? "image" : title.replaceAll("\\.[^.]+$", "");
        s = s.replaceAll("[^A-Za-z0-9_.-]+", "_");
        if (s.length() == 0) s = "image";
        return s;
    }

    public static class Params {
        public String mode = "Auto OD-PCA";
        public double tissueMinOD = 0.06;
        public double whitePercentile = 99.8;
        public double extremePercentile = 1.0;
        public double angleGain = 1.8;
        public double minDisplayAngleDegrees = 55.0;
        public double maxDisplayAngleDegrees = 118.0;
        public double localRadius = 24.0;
        public double localWeight = 0.70;
        public double chromaGain = 2.00;
        public double densityGamma = 0.72;
        public double scalePercentile = 95.0;
        public double displayStrength = 1.85;
        public String customVector1 = "0.650,0.704,0.286";
        public String customVector2 = "0.268,0.570,0.776";
        public boolean showOutputs = true;
        public boolean headless = false;
        public boolean fromMacro = false;
        public String outputDir = "";
        public String prefix = "";

        public static Params fromOptions(String arg) {
            Params p = new Params();
            String opts = Macro.getOptions();
            if (opts == null || opts.length() == 0) opts = arg;
            if (opts == null) opts = "";
            p.fromMacro = opts.length() > 0;
            p.mode = Macro.getValue(opts, "mode", p.mode);
            p.tissueMinOD = getDouble(opts, "tissue", getDouble(opts, "tissueMinOD", p.tissueMinOD));
            p.whitePercentile = getDouble(opts, "white", getDouble(opts, "whitePercentile", p.whitePercentile));
            p.extremePercentile = getDouble(opts, "extreme", getDouble(opts, "extremePercentile", p.extremePercentile));
            p.angleGain = getDouble(opts, "angleGain", p.angleGain);
            p.minDisplayAngleDegrees = getDouble(opts, "minAngle", getDouble(opts, "minDisplayAngle", p.minDisplayAngleDegrees));
            p.maxDisplayAngleDegrees = getDouble(opts, "maxAngle", getDouble(opts, "maxDisplayAngle", p.maxDisplayAngleDegrees));
            p.localRadius = getDouble(opts, "radius", getDouble(opts, "localRadius", p.localRadius));
            p.localWeight = getDouble(opts, "localWeight", p.localWeight);
            p.chromaGain = getDouble(opts, "chroma", getDouble(opts, "chromaGain", p.chromaGain));
            p.densityGamma = getDouble(opts, "gamma", getDouble(opts, "densityGamma", p.densityGamma));
            p.scalePercentile = getDouble(opts, "scale", getDouble(opts, "scalePercentile", p.scalePercentile));
            p.displayStrength = getDouble(opts, "strength", getDouble(opts, "displayStrength", p.displayStrength));
            p.customVector1 = Macro.getValue(opts, "vector1", p.customVector1);
            p.customVector2 = Macro.getValue(opts, "vector2", p.customVector2);
            p.showOutputs = getBoolean(opts, "show", getBoolean(opts, "showOutputs", p.showOutputs));
            p.headless = getBoolean(opts, "headless", p.headless);
            p.outputDir = Macro.getValue(opts, "output", Macro.getValue(opts, "outputDir", p.outputDir));
            p.prefix = Macro.getValue(opts, "prefix", p.prefix);
            return p;
        }

        private static double getDouble(String opts, String key, double fallback) {
            String raw = Macro.getValue(opts, key, null);
            if (raw == null) return fallback;
            try {
                return Double.parseDouble(raw);
            } catch (NumberFormatException e) {
                return fallback;
            }
        }

        private static boolean getBoolean(String opts, String key, boolean fallback) {
            String raw = Macro.getValue(opts, key, null);
            if (raw == null) return fallback;
            return raw.equalsIgnoreCase("true") || raw.equalsIgnoreCase("yes") || raw.equals("1");
        }
    }

    public static class Result {
        public final ImagePlus angleBoostedImage;
        public final ImagePlus stain1EnhancedImage;
        public final ImagePlus stain2EnhancedImage;
        public final ImagePlus stain1DensityImage;
        public final ImagePlus stain2DensityImage;
        public final ResultsTable measurements;

        public Result(ImagePlus angleBoostedImage, ImagePlus stain1EnhancedImage, ImagePlus stain2EnhancedImage,
                      ImagePlus stain1DensityImage, ImagePlus stain2DensityImage, ResultsTable measurements) {
            this.angleBoostedImage = angleBoostedImage;
            this.stain1EnhancedImage = stain1EnhancedImage;
            this.stain2EnhancedImage = stain2EnhancedImage;
            this.stain1DensityImage = stain1DensityImage;
            this.stain2DensityImage = stain2DensityImage;
            this.measurements = measurements;
        }

        public void save(String outputDir, String prefix) throws IOException {
            File dir = new File(outputDir);
            if (!dir.exists() && !dir.mkdirs()) {
                throw new IOException("Could not create " + outputDir);
            }
            saveTiff(angleBoostedImage, new File(dir, prefix + "_angle_boosted_overlay.tif"));
            saveTiff(stain1EnhancedImage, new File(dir, prefix + "_stain1_enhanced_color.tif"));
            saveTiff(stain2EnhancedImage, new File(dir, prefix + "_stain2_enhanced_color.tif"));
            saveTiff(stain1DensityImage, new File(dir, prefix + "_stain1_density.tif"));
            saveTiff(stain2DensityImage, new File(dir, prefix + "_stain2_density.tif"));
            measurements.save(new File(dir, prefix + "_angle_boost_report.csv").getAbsolutePath());
        }

        private static void saveTiff(ImagePlus imp, File file) throws IOException {
            if (!new FileSaver(imp).saveAsTiff(file.getAbsolutePath())) {
                throw new IOException("Could not save " + file.getAbsolutePath());
            }
        }
    }

    private static class StainPair {
        final double[] v1;
        final double[] v2;
        final String mode;

        StainPair(double[] v1, double[] v2, String mode) {
            this.v1 = normalize(v1);
            this.v2 = normalize(v2);
            this.mode = mode;
        }
    }

    private static class Eigen {
        final double[] values;
        final double[][] vectors;

        Eigen(double[] values, double[][] vectors) {
            this.values = values;
            this.vectors = vectors;
        }

        double[] vector(int index) {
            return normalize(vectors[index]);
        }
    }
}
