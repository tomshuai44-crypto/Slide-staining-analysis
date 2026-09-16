import ij.ImagePlus;
import ij.process.ColorProcessor;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.PrintWriter;

/** Whole-image cyclic translations preserve the RGB histogram; compare only an
 * interior farther than the growth/closing/support scales from wrapped edges. */
public class H_DAB_Tissue_Stability {
    public static void main(String[] args) throws Exception {
        BufferedImage bi=ImageIO.read(new File(args[0]));int w=bi.getWidth(),h=bi.getHeight();
        int[] a=bi.getRGB(0,0,w,h,null,0,w);
        H_DAB_Dominance_Extractor.Params p=new H_DAB_Dominance_Extractor.Params();
        H_DAB_Tissue_Selector.Selection base=H_DAB_Tissue_Selector.select(new ImagePlus("base",new ColorProcessor(w,h,a)),p);
        try(PrintWriter out=new PrintWriter(args[1])) {
            out.println("dx,dy,seed,h_threshold,tissue_pixels,interior_pixels,aligned_interior_iou,disagreement_pixels");
            int[][] cases={{1,1},{7,11},{15,15},{0,0}};
            for(int k=0;k<cases.length;k++) {
                int dx=cases[k][0],dy=cases[k][1];int[] shifted=new int[a.length];
                for(int y=0;y<h;y++)for(int x=0;x<w;x++)shifted[((y+dy)%h)*w+(x+dx)%w]=a[y*w+x];
                if(k==3)p.tissueSeed=47;
                H_DAB_Tissue_Selector.Selection r=H_DAB_Tissue_Selector.select(new ImagePlus("shifted",new ColorProcessor(w,h,shifted)),p);
                int margin=p.tissueReachPixels+p.tissueSupportRadius+2*p.tissueCloseRadius+16;
                long intersection=0,union=0,diff=0,area=0;
                for(int y=margin;y<h-margin;y++)for(int x=margin;x<w-margin;x++) {
                    boolean b=base.mask[y*w+x],c=r.mask[(y+dy)*w+x+dx];
                    if(b&&c)intersection++;if(b||c)union++;if(b!=c)diff++;area++;
                }
                double iou=union==0?1:intersection/(double)union;
                out.println(dx+","+dy+","+p.tissueSeed+","+r.hThreshold+","+r.tissuePixels+","+area+","+iou+","+diff);
                System.out.println("offset "+dx+","+dy+" seed "+p.tissueSeed+" interior IoU="+iou+" difference="+diff);
                if(iou<.97)throw new AssertionError("Translation/seed stability below diagnostic bound 0.97");
            }
        }
    }
}
