import ij.ImagePlus;
import ij.process.ColorProcessor;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;

public class H_DAB_Plugin_SmokeTest {
    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            throw new IllegalArgumentException("Usage: H_DAB_Plugin_SmokeTest <input-rgb-image> <output-dir> [prefix]");
        }
        BufferedImage image = ImageIO.read(new File(args[0]));
        if (image == null) {
            throw new IllegalArgumentException("Could not read " + args[0]);
        }
        int width = image.getWidth();
        int height = image.getHeight();
        int[] pixels = new int[width * height];
        image.getRGB(0, 0, width, height, pixels, 0, width);
        ImagePlus imp = new ImagePlus(new File(args[0]).getName(), new ColorProcessor(width, height, pixels));

        H_DAB_Dominance_Extractor.Params params = new H_DAB_Dominance_Extractor.Params();
        params.autoOptimize = true;
        params.showOutputs = false;
        H_DAB_Dominance_Extractor.Result result = H_DAB_Dominance_Extractor.analyze(imp, params);
        if (result == null) {
            throw new IllegalStateException("No result returned");
        }
        String prefix = args.length >= 3 ? args[2] : "smoke";
        result.save(args[1], prefix);
        System.out.println("alpha=" + result.thresholds.alpha);
        System.out.println("dab_min=" + result.thresholds.dabMin);
        System.out.println("dominance_min=" + result.thresholds.dominanceMin);
        System.out.println("rows=" + result.measurements.size());
    }
}
