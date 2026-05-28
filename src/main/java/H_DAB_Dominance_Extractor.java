import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class H_DAB_Dominance_Extractor implements PlugIn {
    private static final double[] HEMA_REF = normalize(new double[] {0.650, 0.704, 0.286});
    private static final double[] DAB_REF = normalize(new double[] {0.268, 0.570, 0.776});

    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("H-DAB Dominance Extractor", "Open an RGB H-DAB image first.");
            return;
        }
        if (imp.getType() != ImagePlus.COLOR_RGB) {
            IJ.error("H-DAB Dominance Extractor", "This plugin expects an RGB image.");
            return;
        }

        Params params = Params.fromOptions(arg);
        if (!params.headless && !params.fromMacro && !showDialog(params)) {
            return;
        }

        Result result = analyze(imp, params);
        if (result == null) {
            IJ.error("H-DAB Dominance Extractor", "Could not analyze this image.");
            return;
        }

        if (params.showOutputs) {
            result.hemaColorImage.show();
            result.dabColorImage.show();
            result.residualColorImage.show();
            result.hemaImage.show();
            result.dabImage.show();
            result.residualImage.show();
            result.dominanceImage.show();
            result.maskImage.show();
            result.measurements.show("H-DAB Dominance Results");
        }

        if (params.outputDir.length() > 0) {
            try {
                result.save(params.outputDir, params.prefix.length() > 0 ? params.prefix : sanitizeTitle(imp.getTitle()));
            } catch (IOException e) {
                IJ.error("H-DAB Dominance Extractor", "Could not save outputs:\n" + e.getMessage());
            }
        }
    }

    public static Result analyze(ImagePlus imp, Params params) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        int n = width * height;
        int[] pixels = (int[]) imp.getProcessor().convertToRGB().getPixels();

        double[] hema = new double[n];
        double[] dab = new double[n];
        double[] residual = new double[n];
        double[] odNorm = new double[n];

        double[][] inv = inverse(stainMatrix());
        for (int i = 0; i < n; i++) {
            int rgb = pixels[i];
            double r = ((rgb >> 16) & 0xff);
            double g = ((rgb >> 8) & 0xff);
            double b = (rgb & 0xff);
            double odR = -Math.log((r + 1.0) / 256.0);
            double odG = -Math.log((g + 1.0) / 256.0);
            double odB = -Math.log((b + 1.0) / 256.0);
            odNorm[i] = Math.sqrt(odR * odR + odG * odG + odB * odB);
            double c0 = inv[0][0] * odR + inv[0][1] * odG + inv[0][2] * odB;
            double c1 = inv[1][0] * odR + inv[1][1] * odG + inv[1][2] * odB;
            double c2 = inv[2][0] * odR + inv[2][1] * odG + inv[2][2] * odB;
            hema[i] = Math.max(0.0, c0);
            dab[i] = Math.max(0.0, c1);
            residual[i] = Math.max(0.0, c2);
        }

        boolean[] tissue = new boolean[n];
        int tissueCount = 0;
        for (int i = 0; i < n; i++) {
            tissue[i] = odNorm[i] >= params.tissueMinOD;
            if (tissue[i]) tissueCount++;
        }
        if (tissueCount < 10) {
            return null;
        }

        double hScale = Math.max(percentile(hema, tissue, 99.0), 1e-9);
        double dScale = Math.max(percentile(dab, tissue, 99.0), 1e-9);
        double rScale = Math.max(percentile(residual, tissue, 99.0), 1e-9);

        double[] hNorm = normalizeByScale(hema, hScale);
        double[] dNorm = normalizeByScale(dab, dScale);
        double[] rNorm = normalizeByScale(residual, rScale);

        ChosenThresholds chosen = params.autoOptimize
                ? optimizeThresholds(tissue, hNorm, dNorm, rNorm, params)
                : new ChosenThresholds(params.alpha, params.dabMin, params.dominanceMin, params.residualMax, params.ratioMin, 0.0);

        double[] dominance = new double[n];
        byte[] mask = new byte[n];
        int positiveCount = 0;
        double sumDab = 0.0;
        double sumHema = 0.0;
        double sumResidual = 0.0;
        for (int i = 0; i < n; i++) {
            dominance[i] = dNorm[i] - chosen.alpha * hNorm[i];
            boolean positive = tissue[i]
                    && dNorm[i] >= chosen.dabMin
                    && dominance[i] >= chosen.dominanceMin
                    && dNorm[i] >= chosen.ratioMin * hNorm[i]
                    && rNorm[i] <= chosen.residualMax;
            if (positive) {
                mask[i] = (byte) 255;
                positiveCount++;
                sumDab += dab[i];
                sumHema += hema[i];
                sumResidual += residual[i];
            }
        }

        String base = sanitizeTitle(imp.getTitle());
        ImagePlus hemaColorImage = new ImagePlus(base + " - Hematoxylin color deconvolution", colorDeconvolvedProcessor(width, height, hema, HEMA_REF));
        ImagePlus dabColorImage = new ImagePlus(base + " - DAB color deconvolution", colorDeconvolvedProcessor(width, height, dab, DAB_REF));
        ImagePlus residualColorImage = new ImagePlus(base + " - Residual color deconvolution", colorDeconvolvedProcessor(width, height, residual, residualVector()));
        ImagePlus hemaImage = new ImagePlus(base + " - Hematoxylin concentration", floatProcessor(width, height, hema));
        ImagePlus dabImage = new ImagePlus(base + " - DAB concentration", floatProcessor(width, height, dab));
        ImagePlus residualImage = new ImagePlus(base + " - Residual concentration", floatProcessor(width, height, residual));
        ImagePlus dominanceImage = new ImagePlus(base + " - DAB dominance score", floatProcessor(width, height, dominance));
        ImagePlus maskImage = new ImagePlus(base + " - DAB positive mask", new ByteProcessor(width, height, mask));

        ResultsTable rt = new ResultsTable();
        rt.incrementCounter();
        rt.addValue("image", imp.getTitle());
        rt.addValue("width", width);
        rt.addValue("height", height);
        rt.addValue("tissue_pixels", tissueCount);
        rt.addValue("dab_positive_pixels", positiveCount);
        rt.addValue("dab_area_fraction", positiveCount / (double) tissueCount);
        rt.addValue("mean_dab_concentration", positiveCount > 0 ? sumDab / positiveCount : 0.0);
        rt.addValue("integrated_dab_concentration", sumDab);
        rt.addValue("mean_hema_in_positive", positiveCount > 0 ? sumHema / positiveCount : 0.0);
        rt.addValue("mean_residual_in_positive", positiveCount > 0 ? sumResidual / positiveCount : 0.0);
        rt.addValue("alpha", chosen.alpha);
        rt.addValue("dab_min_normalized", chosen.dabMin);
        rt.addValue("dominance_min", chosen.dominanceMin);
        rt.addValue("residual_max_normalized", chosen.residualMax);
        rt.addValue("dab_to_hema_ratio_min", chosen.ratioMin);
        rt.addValue("optimization_score", chosen.score);
        rt.addValue("hema_scale_p99", hScale);
        rt.addValue("dab_scale_p99", dScale);
        rt.addValue("residual_scale_p99", rScale);

        return new Result(hemaColorImage, dabColorImage, residualColorImage,
                hemaImage, dabImage, residualImage, dominanceImage, maskImage, rt, chosen);
    }

    private static boolean showDialog(Params params) {
        GenericDialog gd = new GenericDialog("H-DAB Dominance Extractor");
        gd.addNumericField("Alpha: DAB - alpha * Hematoxylin", params.alpha, 2);
        gd.addNumericField("Tissue minimum OD", params.tissueMinOD, 3);
        gd.addNumericField("DAB minimum, normalized", params.dabMin, 3);
        gd.addNumericField("Dominance minimum", params.dominanceMin, 3);
        gd.addNumericField("DAB/Hematoxylin ratio minimum", params.ratioMin, 2);
        gd.addNumericField("Residual maximum, normalized", params.residualMax, 3);
        gd.addCheckbox("Auto-optimize thresholds", params.autoOptimize);
        gd.addCheckbox("Show output images", params.showOutputs);
        gd.addStringField("Output directory (optional)", params.outputDir, 32);
        gd.addStringField("Output prefix (optional)", params.prefix, 24);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        params.alpha = gd.getNextNumber();
        params.tissueMinOD = gd.getNextNumber();
        params.dabMin = gd.getNextNumber();
        params.dominanceMin = gd.getNextNumber();
        params.ratioMin = gd.getNextNumber();
        params.residualMax = gd.getNextNumber();
        params.autoOptimize = gd.getNextBoolean();
        params.showOutputs = gd.getNextBoolean();
        params.outputDir = gd.getNextString();
        params.prefix = gd.getNextString();
        return true;
    }

    private static ChosenThresholds optimizeThresholds(boolean[] tissue, double[] hNorm, double[] dNorm, double[] rNorm, Params params) {
        double baseDab = clamp(otsuThreshold(dNorm, tissue) * 0.45, 0.02, 0.42);
        double residualMax = params.residualMax;
        boolean[] dabDominant = new boolean[dNorm.length];
        int dominantCount = 0;
        for (int i = 0; i < dNorm.length; i++) {
            dabDominant[i] = tissue[i]
                    && dNorm[i] >= params.ratioMin * hNorm[i]
                    && rNorm[i] <= residualMax;
            if (dabDominant[i]) dominantCount++;
        }
        if (dominantCount >= 50) {
            baseDab = Math.max(baseDab, percentile(dNorm, dabDominant, 85.0));
        }
        ChosenThresholds best = null;
        for (double alpha = 0.20; alpha <= 1.201; alpha += 0.05) {
            double[] dominance = new double[dNorm.length];
            boolean[] candidateArea = new boolean[dNorm.length];
            for (int i = 0; i < dNorm.length; i++) {
                dominance[i] = dNorm[i] - alpha * hNorm[i];
                candidateArea[i] = tissue[i] && dNorm[i] >= Math.max(0.015, baseDab * 0.35) && rNorm[i] <= residualMax;
            }
            double domMin = clamp(otsuThreshold(dominance, candidateArea) * 0.50, -0.08, 0.58);
            ChosenThresholds chosen = scoreThresholds(tissue, hNorm, dNorm, rNorm, alpha, baseDab, domMin, residualMax, params.ratioMin);
            if (best == null || chosen.score > best.score) {
                best = chosen;
            }
        }
        if (best == null) {
            return new ChosenThresholds(params.alpha, params.dabMin, params.dominanceMin, params.residualMax, params.ratioMin, 0.0);
        }
        return best;
    }

    private static ChosenThresholds scoreThresholds(boolean[] tissue, double[] hNorm, double[] dNorm, double[] rNorm,
                                                   double alpha, double dabMin, double dominanceMin, double residualMax,
                                                   double ratioMin) {
        int tissueCount = 0;
        int positiveCount = 0;
        double posD = 0.0, negD = 0.0, posH = 0.0, posR = 0.0, posDom = 0.0, negDom = 0.0;
        int negCount = 0;
        for (int i = 0; i < tissue.length; i++) {
            if (!tissue[i]) continue;
            tissueCount++;
            double dom = dNorm[i] - alpha * hNorm[i];
            boolean pos = dNorm[i] >= dabMin
                    && dom >= dominanceMin
                    && dNorm[i] >= ratioMin * hNorm[i]
                    && rNorm[i] <= residualMax;
            if (pos) {
                positiveCount++;
                posD += dNorm[i];
                posH += hNorm[i];
                posR += rNorm[i];
                posDom += dom;
            } else {
                negCount++;
                negD += dNorm[i];
                negDom += dom;
            }
        }
        if (positiveCount < 5 || tissueCount < 10) {
            return new ChosenThresholds(alpha, dabMin, dominanceMin, residualMax, ratioMin, -1e9);
        }
        double frac = positiveCount / (double) tissueCount;
        double fracPenalty = 0.0;
        if (frac < 0.002) fracPenalty += (0.002 - frac) * 40.0;
        if (frac > 0.45) fracPenalty += (frac - 0.45) * 6.0;
        double meanPosD = posD / positiveCount;
        double meanNegD = negCount > 0 ? negD / negCount : 0.0;
        double meanPosH = posH / positiveCount;
        double meanPosR = posR / positiveCount;
        double meanPosDom = posDom / positiveCount;
        double meanNegDom = negCount > 0 ? negDom / negCount : 0.0;
        double score = (meanPosD - meanNegD)
                + 0.55 * (meanPosDom - meanNegDom)
                - 0.22 * meanPosH
                - 0.18 * meanPosR
                - fracPenalty;
        return new ChosenThresholds(alpha, dabMin, dominanceMin, residualMax, ratioMin, score);
    }

    private static FloatProcessor floatProcessor(int width, int height, double[] values) {
        float[] out = new float[values.length];
        for (int i = 0; i < values.length; i++) out[i] = (float) values[i];
        FloatProcessor fp = new FloatProcessor(width, height, out);
        fp.resetMinAndMax();
        return fp;
    }

    private static ColorProcessor colorDeconvolvedProcessor(int width, int height, double[] amount, double[] vector) {
        int[] out = new int[amount.length];
        for (int i = 0; i < amount.length; i++) {
            int r = clamp255(255.0 * Math.exp(-amount[i] * vector[0]));
            int g = clamp255(255.0 * Math.exp(-amount[i] * vector[1]));
            int b = clamp255(255.0 * Math.exp(-amount[i] * vector[2]));
            out[i] = (r << 16) | (g << 8) | b;
        }
        return new ColorProcessor(width, height, out);
    }

    private static double[] normalizeByScale(double[] values, double scale) {
        double[] out = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = clamp(values[i] / scale, 0.0, 1.5);
        }
        return out;
    }

    private static double percentile(double[] values, boolean[] mask, double pct) {
        ArrayList<Double> vals = new ArrayList<Double>();
        for (int i = 0; i < values.length; i++) {
            if (mask[i] && Double.isFinite(values[i])) vals.add(values[i]);
        }
        if (vals.isEmpty()) return 0.0;
        double[] arr = new double[vals.size()];
        for (int i = 0; i < vals.size(); i++) arr[i] = vals.get(i);
        Arrays.sort(arr);
        double rank = (pct / 100.0) * (arr.length - 1);
        int lo = (int) Math.floor(rank);
        int hi = (int) Math.ceil(rank);
        if (lo == hi) return arr[lo];
        double frac = rank - lo;
        return arr[lo] * (1.0 - frac) + arr[hi] * frac;
    }

    private static double otsuThreshold(double[] values, boolean[] mask) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        int count = 0;
        for (int i = 0; i < values.length; i++) {
            if (!mask[i] || !Double.isFinite(values[i])) continue;
            min = Math.min(min, values[i]);
            max = Math.max(max, values[i]);
            count++;
        }
        if (count < 10 || !(max > min)) return 0.05;
        int bins = 256;
        double[] hist = new double[bins];
        for (int i = 0; i < values.length; i++) {
            if (!mask[i] || !Double.isFinite(values[i])) continue;
            int bin = (int) Math.floor((values[i] - min) / (max - min) * (bins - 1));
            if (bin < 0) bin = 0;
            if (bin >= bins) bin = bins - 1;
            hist[bin] += 1.0;
        }
        double total = count;
        double sum = 0.0;
        for (int i = 0; i < bins; i++) sum += i * hist[i];
        double sumB = 0.0;
        double wB = 0.0;
        double bestVar = -1.0;
        int best = 0;
        for (int i = 0; i < bins; i++) {
            wB += hist[i];
            if (wB <= 0.0) continue;
            double wF = total - wB;
            if (wF <= 0.0) break;
            sumB += i * hist[i];
            double mB = sumB / wB;
            double mF = (sum - sumB) / wF;
            double between = wB * wF * (mB - mF) * (mB - mF);
            if (between > bestVar) {
                bestVar = between;
                best = i;
            }
        }
        return min + (best / (double) (bins - 1)) * (max - min);
    }

    private static double[][] stainMatrix() {
        double[] residual = residualVector();
        return new double[][] {
                {HEMA_REF[0], DAB_REF[0], residual[0]},
                {HEMA_REF[1], DAB_REF[1], residual[1]},
                {HEMA_REF[2], DAB_REF[2], residual[2]},
        };
    }

    private static double[] residualVector() {
        double[] residual = new double[3];
        for (int i = 0; i < 3; i++) {
            residual[i] = Math.sqrt(Math.max(0.0, 1.0 - HEMA_REF[i] * HEMA_REF[i] - DAB_REF[i] * DAB_REF[i]));
        }
        return normalize(residual);
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

    private static double[] normalize(double[] v) {
        double norm = 0.0;
        for (double x : v) norm += x * x;
        norm = Math.sqrt(norm);
        double[] out = new double[v.length];
        for (int i = 0; i < v.length; i++) out[i] = v[i] / Math.max(norm, 1e-12);
        return out;
    }

    private static double clamp(double x, double lo, double hi) {
        return Math.max(lo, Math.min(hi, x));
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
        public double alpha = 0.70;
        public double tissueMinOD = 0.08;
        public double dabMin = 0.10;
        public double dominanceMin = 0.00;
        public double ratioMin = 1.15;
        public double residualMax = 1.05;
        public boolean autoOptimize = true;
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
            p.alpha = getDouble(opts, "alpha", p.alpha);
            p.tissueMinOD = getDouble(opts, "tissue", getDouble(opts, "tissueMinOD", p.tissueMinOD));
            p.dabMin = getDouble(opts, "dab", getDouble(opts, "dabMin", p.dabMin));
            p.dominanceMin = getDouble(opts, "dominance", getDouble(opts, "dominanceMin", p.dominanceMin));
            p.ratioMin = getDouble(opts, "ratio", getDouble(opts, "ratioMin", p.ratioMin));
            p.residualMax = getDouble(opts, "residual", getDouble(opts, "residualMax", p.residualMax));
            p.autoOptimize = getBoolean(opts, "auto", getBoolean(opts, "autoOptimize", p.autoOptimize));
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

    public static class ChosenThresholds {
        public final double alpha;
        public final double dabMin;
        public final double dominanceMin;
        public final double residualMax;
        public final double ratioMin;
        public final double score;

        public ChosenThresholds(double alpha, double dabMin, double dominanceMin, double residualMax, double ratioMin, double score) {
            this.alpha = alpha;
            this.dabMin = dabMin;
            this.dominanceMin = dominanceMin;
            this.residualMax = residualMax;
            this.ratioMin = ratioMin;
            this.score = score;
        }
    }

    public static class Result {
        public final ImagePlus hemaColorImage;
        public final ImagePlus dabColorImage;
        public final ImagePlus residualColorImage;
        public final ImagePlus hemaImage;
        public final ImagePlus dabImage;
        public final ImagePlus residualImage;
        public final ImagePlus dominanceImage;
        public final ImagePlus maskImage;
        public final ResultsTable measurements;
        public final ChosenThresholds thresholds;

        public Result(ImagePlus hemaColorImage, ImagePlus dabColorImage, ImagePlus residualColorImage,
                      ImagePlus hemaImage, ImagePlus dabImage, ImagePlus residualImage,
                      ImagePlus dominanceImage, ImagePlus maskImage, ResultsTable measurements,
                      ChosenThresholds thresholds) {
            this.hemaColorImage = hemaColorImage;
            this.dabColorImage = dabColorImage;
            this.residualColorImage = residualColorImage;
            this.hemaImage = hemaImage;
            this.dabImage = dabImage;
            this.residualImage = residualImage;
            this.dominanceImage = dominanceImage;
            this.maskImage = maskImage;
            this.measurements = measurements;
            this.thresholds = thresholds;
        }

        public void save(String outputDir, String prefix) throws IOException {
            File dir = new File(outputDir);
            if (!dir.exists() && !dir.mkdirs()) {
                throw new IOException("Could not create " + outputDir);
            }
            saveTiff(hemaColorImage, new File(dir, prefix + "_hematoxylin_color.tif"));
            saveTiff(dabColorImage, new File(dir, prefix + "_dab_color.tif"));
            saveTiff(residualColorImage, new File(dir, prefix + "_residual_color.tif"));
            saveTiff(hemaImage, new File(dir, prefix + "_hematoxylin.tif"));
            saveTiff(dabImage, new File(dir, prefix + "_dab.tif"));
            saveTiff(residualImage, new File(dir, prefix + "_residual.tif"));
            saveTiff(dominanceImage, new File(dir, prefix + "_dab_dominance.tif"));
            saveTiff(maskImage, new File(dir, prefix + "_dab_positive_mask.tif"));
            measurements.save(new File(dir, prefix + "_measurements.csv").getAbsolutePath());
        }

        private static void saveTiff(ImagePlus imp, File file) throws IOException {
            if (!new FileSaver(imp).saveAsTiff(file.getAbsolutePath())) {
                throw new IOException("Could not save " + file.getAbsolutePath());
            }
        }
    }
}
