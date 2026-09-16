import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.lang.reflect.Field;
import java.util.Arrays;

/** Same runner works against the original JAR and the updated JAR. */
public class H_DAB_Tissue_QA {
    public static void main(String[] args) throws Exception {
        BufferedImage image=ImageIO.read(new File(args[0]));int w=image.getWidth(),h=image.getHeight();
        int[] rgb=image.getRGB(0,0,w,h,null,0,w);ImagePlus im=new ImagePlus(new File(args[0]).getName(),new ColorProcessor(w,h,rgb));
        H_DAB_Dominance_Extractor.Params p=new H_DAB_Dominance_Extractor.Params();p.showOutputs=false;
        if(args.length>3)p=H_DAB_Dominance_Extractor.Params.fromOptions(args[3]);
        long start=System.nanoTime();H_DAB_Dominance_Extractor.Result r=H_DAB_Dominance_Extractor.analyze(im,p);
        String dir=args[1],prefix=args[2];new File(dir).mkdirs();
        if(r==null)throw new IllegalStateException("No analysis result");
        ImagePlus tissue,overlay;
        try {
            tissue=(ImagePlus)r.getClass().getField("tissueMaskImage").get(r);
            overlay=(ImagePlus)r.getClass().getField("tissueOverlayImage").get(r);
            new FileSaver((ImagePlus)r.getClass().getField("tissueSeedImage").get(r)).saveAsTiff(dir+"/"+prefix+"_tissue_h_evidence.tif");
        } catch(NoSuchFieldException e) {
            byte[] mask=new byte[w*h];int[] over=rgb.clone();
            for(int i=0;i<rgb.length;i++) {
                double od2=0;for(int shift:new int[]{16,8,0}){double od=-Math.log((((rgb[i]>>shift)&255)+1.0)/256);od2+=od*od;}
                if(Math.sqrt(od2)>=p.tissueMinOD){mask[i]=(byte)255;int c=rgb[i];over[i]=((int)(.65*((c>>16)&255))<<16)|((int)(.65*((c>>8)&255)+.35*255)<<8)|(int)(.65*(c&255));}
            }
            tissue=new ImagePlus("legacy OD mask reconstructed from exact old predicate",new ByteProcessor(w,h,mask));
            overlay=new ImagePlus("legacy overlay",new ColorProcessor(w,h,over));
            int count=0;for(byte b:mask)if(b!=0)count++;
            if(count!=r.measurements.getValue("tissue_pixels",0))throw new AssertionError("Legacy mask disagrees with installed Java count");
        }
        new FileSaver(tissue).saveAsTiff(dir+"/"+prefix+"_tissue_mask.tif");
        new FileSaver(overlay).saveAsTiff(dir+"/"+prefix+"_tissue_overlay.tif");
        new FileSaver(r.maskImage).saveAsTiff(dir+"/"+prefix+"_dab_positive_mask.tif");
        r.measurements.setPrecision(9);r.measurements.addValue("qa_elapsed_seconds",(System.nanoTime()-start)/1e9);
        r.measurements.addValue("qa_raw_h_hash",Integer.toString(Arrays.hashCode((float[])r.hemaImage.getProcessor().getPixels())));
        r.measurements.addValue("qa_raw_d_hash",Integer.toString(Arrays.hashCode((float[])r.dabImage.getProcessor().getPixels())));
        r.measurements.addValue("qa_raw_r_hash",Integer.toString(Arrays.hashCode((float[])r.residualImage.getProcessor().getPixels())));
        r.measurements.save(dir+"/"+prefix+"_measurements.csv");
        System.out.println(prefix+" tissue="+r.measurements.getValue("tissue_pixels",0)+" positive="+r.measurements.getValue("dab_positive_pixels",0)+" seconds="+r.measurements.getValue("qa_elapsed_seconds",0));
    }
}
