import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import java.util.Arrays;
import java.util.Random;

/** Regression tests exercise the Java selector and complete measurement path. */
public final class H_DAB_Local_Tissue_Test {
    private static int checks=0;
    private static final int W=360,H=160;
    private static void check(boolean ok,String message) {
        if(!ok)throw new AssertionError(message);
        checks++;System.out.println("PASS LOCAL "+message);
    }
    private static H_DAB_Dominance_Extractor.Params params() {
        H_DAB_Dominance_Extractor.Params p=new H_DAB_Dominance_Extractor.Params();
        p.showOutputs=false;p.autoOptimize=false;p.tissueAutoWhite=false;
        p.tissueWhiteR=205;p.tissueWhiteG=195;p.tissueWhiteB=198;return p;
    }
    private static int pixel(double h,double d) {
        double hn=Math.sqrt(.650*.650+.704*.704+.286*.286),dn=Math.sqrt(.268*.268+.570*.570+.776*.776);
        double[] hv={.650/hn,.704/hn,.286/hn},dv={.268/dn,.570/dn,.776/dn},white={205,195,198};int result=0;
        for(int c=0;c<3;c++)result=(result<<8)|(int)Math.max(0,Math.min(255,Math.round((white[c]+1)*Math.exp(-h*hv[c]-d*dv[c])-1)));
        return result;
    }
    private static ImagePlus blank(int width,int height,int color) {
        int[] a=new int[width*height];Arrays.fill(a,color);return new ImagePlus("local_test",new ColorProcessor(width,height,a));
    }
    private static boolean truth(int x,int y) {return y>=40&&y<120&&x>=20&&x<330&&!(y>=65&&y<95&&x>=130&&x<180);}
    private static ImagePlus phantom() {
        ImagePlus im=blank(W,H,pixel(0,0));int[] a=(int[])im.getProcessor().getPixels();
        for(int y=0;y<H;y++)for(int x=0;x<W;x++)if(truth(x,y))a[y*W+x]=pixel(x<60?.5:0,x>=50?.45:0);
        return im;
    }
    public static void main(String[] args) {
        H_DAB_Dominance_Extractor.Params p=params();
        check(p.tissueModel.equals("local")&&p.tissueLocalK==8,"default model and reviewed k=8");
        ImagePlus im=phantom();boolean[] previous=null;
        H_DAB_Tissue_Selector.Selection selected=null;
        for(double k:new double[]{2,3,4,5,8}) {
            p.tissueLocalK=k;selected=H_DAB_Local_Tissue_Selector.select(im,p);
            int intersection=0,union=0;
            for(int y=0;y<H;y++)for(int x=0;x<W;x++) {
                boolean actual=selected.mask[y*W+x],truth=truth(x,y);
                if(actual&&truth)intersection++;if(actual||truth)union++;
            }
            double iou=intersection/(double)union;
            check(iou>.97,"k="+k+" analytic H/D extension shape IoU="+iou);
            check(selected.mask[80*W+315],"k="+k+" DAB-only extension >250 pixels beyond H retained");
            check(selected.mask[60*W+30],"k="+k+" DAB-negative H tissue retained");
            check(!selected.mask[80*W+150],"k="+k+" clear lumen preserved");
            if(previous!=null){boolean nested=true;for(int i=0;i<previous.length;i++)if(selected.mask[i]&&!previous[i])nested=false;check(nested,"nested masks as k increases");}
            previous=selected.mask;
        }
        H_DAB_Tissue_Selector.Selection repeat=H_DAB_Local_Tissue_Selector.select(im,p);
        check(Arrays.equals(selected.mask,repeat.mask)&&selected.diagnostics.equals(repeat.diagnostics),"fixed-seed reproducibility including local statistics");
        p.tissueReachPixels=0;
        check(Arrays.equals(selected.mask,H_DAB_Local_Tissue_Selector.select(im,p).mask),"previous H-only reach does not constrain local model");
        for(int color:new int[]{0xffffff,0x808080,0xcdc3c6,0x000000}) {
            check(H_DAB_Local_Tissue_Selector.select(blank(80,80,color),p).tissuePixels==0,"blank/no-H color "+Integer.toHexString(color)+" does not force a mask");
        }
        check(H_DAB_Local_Tissue_Selector.select(blank(80,80,pixel(0,.7)),p).tissuePixels==0,"DAB alone without reliable H cannot seed tissue");
        ImagePlus artifact=phantom();int[] aa=(int[])artifact.getProcessor().getPixels();
        for(int y=5;y<20;y++)for(int x=290;x<310;x++)aa[y*W+x]=pixel(0,.8);
        check(!H_DAB_Local_Tissue_Selector.select(artifact,p).mask[10*W+300],"isolated unseeded DAB island excluded");
        for(int y=5;y<20;y++)for(int x=290;x<310;x++)aa[y*W+x]=pixel(.7,.2);
        check(H_DAB_Local_Tissue_Selector.select(artifact,p).mask[10*W+300],"known limitation: seeded H-coloured artifact remains");
        ImagePlus noise=blank(W,H,pixel(0,0));int[] nn=(int[])noise.getProcessor().getPixels();Random random=new Random(98);
        for(int i=0;i<nn.length;i++)nn[i]=((205+random.nextInt(7)-3)<<16)|((195+random.nextInt(7)-3)<<8)|(198+random.nextInt(7)-3);
        for(int y=40;y<80;y++)for(int x=20;x<60;x++)nn[y*W+x]=pixel(.5,0);
        boolean[] nm=H_DAB_Local_Tissue_Selector.select(noise,p).mask;int far=0;
        for(int y=0;y<H;y++)for(int x=150;x<W;x++)if(nm[y*W+x])far++;
        check(far==0,"H patch does not flood independent neutral RGB noise");
        im=phantom();im.setRoi(new OvalRoi(5,25,230,110));selected=H_DAB_Local_Tissue_Selector.select(im,p);
        boolean inside=true;for(int y=0;y<H;y++)for(int x=0;x<W;x++)if(selected.mask[y*W+x]&&!im.getRoi().contains(x,y))inside=false;
        check(inside&&selected.tissuePixels>0,"nonrectangular area ROI clips local mask");
        int[] a=(int[])im.getProcessor().getPixels();for(int y=0;y<H;y++)for(int x=0;x<W;x++)if(!im.getRoi().contains(x,y))a[y*W+x]=pixel(1.2,2);
        repeat=H_DAB_Local_Tissue_Selector.select(im,p);
        check(Arrays.equals(selected.mask,repeat.mask)&&selected.diagnostics.equals(repeat.diagnostics),"pixels outside ROI cannot alter local mask or background thresholds");
        im=phantom();im.setRoi(new Roi(200,50,80,50));
        check(H_DAB_Local_Tissue_Selector.select(im,p).tissuePixels==0,"ROI cannot borrow an H seed outside its domain");
        H_DAB_Dominance_Extractor.Result empty=H_DAB_Dominance_Extractor.analyze(blank(7,7,pixel(.5,0)),p);
        check(empty.tissueSelection.status.equals("FAILED_LOCAL_BACKGROUND")&&Double.isNaN(empty.measurements.getValue("dab_area_fraction",0)),"insufficient local background explicitly fails with NaN denominator");
        im=phantom();p.tissueModel="local";
        H_DAB_Dominance_Extractor.Result result=H_DAB_Dominance_Extractor.analyze(im,p);
        p.tissueModel="h_only";p.tissueReachPixels=96;
        H_DAB_Dominance_Extractor.Result legacy=H_DAB_Dominance_Extractor.analyze(im,p);
        check(Arrays.equals((float[])result.hemaImage.getProcessor().getPixels(),(float[])legacy.hemaImage.getProcessor().getPixels())
                &&Arrays.equals((float[])result.dabImage.getProcessor().getPixels(),(float[])legacy.dabImage.getProcessor().getPixels())
                &&Arrays.equals((float[])result.residualImage.getProcessor().getPixels(),(float[])legacy.residualImage.getProcessor().getPixels()),"raw measurement H/D/residual maps unchanged across tissue models");
        byte[] mask=(byte[])result.tissueMaskImage.getProcessor().getPixels(),positive=(byte[])result.maskImage.getProcessor().getPixels();int count=0,pos=0;boolean subset=true;
        for(int i=0;i<mask.length;i++){if(mask[i]!=0)count++;if(positive[i]!=0){pos++;if(mask[i]==0)subset=false;}}
        check(subset&&count==result.measurements.getValue("tissue_pixels",0),"saved tissue mask equals denominator and contains all DAB positives");
        check(Math.abs(pos/(double)count-result.measurements.getValue("dab_area_fraction",0))<1e-12,"local tissue area is used by measurement fraction");
        check(result.measurements.getStringValue("tissue_algorithm",0).equals(H_DAB_Local_Tissue_Selector.VERSION)
                &&Double.isNaN(result.measurements.getValue("tissue_reach_px",0))
                &&result.measurements.getValue("tissue_local_k",0)==8,"CSV identifies local model, threshold and absence of reach cap");
        p=H_DAB_Dominance_Extractor.Params.fromOptions("tissuemodel=local tissuelocalk=5 tissuelocalminod=0.03");
        check(p.tissueModel.equals("local")&&p.tissueLocalK==5&&p.tissueLocalMinOD==.03,"local macro options parsed");
        p.tissueLocalK=Double.NaN;boolean rejected=false;
        try{H_DAB_Local_Tissue_Selector.select(im,p);}catch(IllegalArgumentException e){rejected=true;}
        check(rejected,"invalid local threshold rejected");
        System.out.println("ALL "+checks+" LOCAL CHECKS PASSED");
    }
}
