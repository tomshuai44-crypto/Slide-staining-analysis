import ij.IJ;
import ij.ImagePlus;

public class Angle_Booster_SmokeTest {
    public static void main(String[] args) throws Exception {
        String input = args[0];
        String output = args[1];
        String prefix = args.length > 2 ? args[2] : "smoke";
        String mode = args.length > 3 ? args[3] : "H-DAB";
        ImagePlus imp = IJ.openImage(input);
        if (imp == null) throw new IllegalArgumentException("Could not open " + input);
        Stain_Angle_Contrast_Booster.Params params = new Stain_Angle_Contrast_Booster.Params();
        params.mode = mode;
        params.showOutputs = false;
        Stain_Angle_Contrast_Booster.Result result = Stain_Angle_Contrast_Booster.analyze(imp, params);
        if (result == null) throw new IllegalStateException("No result");
        result.save(output, prefix);
        System.out.println("saved " + output);
        System.exit(0);
    }
}
