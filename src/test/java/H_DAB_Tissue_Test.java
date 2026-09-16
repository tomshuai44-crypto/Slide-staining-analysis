import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import java.util.Arrays;
import java.util.Random;

/** Run against the built/installed JAR, not a Python reimplementation. */
public class H_DAB_Tissue_Test {
    private static int checks=0;
    private static void check(boolean ok,String name){if(!ok)throw new AssertionError(name);checks++;System.out.println("PASS "+name);}
    private static int stain(double h,double d,double gray){
        double hn=Math.sqrt(.650*.650+.704*.704+.286*.286),dn=Math.sqrt(.268*.268+.570*.570+.776*.776);
        double[] hv={.650/hn,.704/hn,.286/hn},dv={.268/dn,.570/dn,.776/dn};int rgb=0;
        for(int c=0;c<3;c++)rgb=(rgb<<8)|Math.max(0,Math.min(255,(int)Math.round(256*Math.exp(-h*hv[c]-d*dv[c]-gray)-1)));
        return rgb;
    }
    private static ImagePlus filled(int w,int h,int c){int[] a=new int[w*h];Arrays.fill(a,c);return new ImagePlus("synthetic",new ColorProcessor(w,h,a));}
    private static H_DAB_Dominance_Extractor.Params params(){H_DAB_Dominance_Extractor.Params p=new H_DAB_Dominance_Extractor.Params();p.autoOptimize=false;p.showOutputs=false;p.tissueModel="h_only";return p;}
    private static H_DAB_Tissue_Selector.Selection select(ImagePlus im){return H_DAB_Tissue_Selector.select(im,params());}
    private static ImagePlus phantom(){
        ImagePlus im=filled(256,192,0xffffff);int[] a=(int[])im.getProcessor().getPixels();
        for(int y=0;y<192;y++)for(int x=0;x<256;x++) {
            boolean t=(x<96 && y>=16 && y<160)||(x>=176 && y>=112);
            boolean hole=x>=32&&x<72&&y>=64&&y<112;
            if(t&&!hole)a[y*256+x]=stain(.025,.04,.07);
            if(t&&!hole && x%16<5 && y%16<5)a[y*256+x]=stain(.55,.06,.03);
        }
        return im;
    }
    public static void main(String[] args) {
        for(int c:new int[]{0xffffff,0xf0f0f0,0x808080,0xcebfc2,0x000000}) {
            H_DAB_Tissue_Selector.Selection s=select(filled(128,128,c));check(s.tissuePixels==0,"background-only "+Integer.toHexString(c));
        }
        check(select(filled(128,128,stain(0,.65,0))).tissuePixels==0,"DAB alone never supplies tissue denominator");
        H_DAB_Dominance_Extractor.Params manual=params();manual.tissueAutoWhite=false;
        check(H_DAB_Tissue_Selector.select(filled(128,128,0x808080),manual).tissuePixels==0,"neutral darkness rejected even with fixed white reference");
        check(H_DAB_Tissue_Selector.select(filled(128,128,stain(0,.65,0)),manual).tissuePixels==0,"pure DAB rejected even with fixed white reference");
        check(H_DAB_Tissue_Selector.select(filled(128,128,stain(.4,0,0)),manual).tissuePixels==128*128,"known reference permits an entirely H-stained frame");
        check(select(filled(128,128,stain(.4,0,0))).tissuePixels==0,"uniform full-frame stain without a reference is not guessed");
        ImagePlus noise=filled(256,192,0xffffff);int[] na=(int[])noise.getProcessor().getPixels();
        for(int y=4;y<192;y+=11)for(int x=4;x<256;x+=11)na[y*256+x]=stain(.8,0,0);
        check(select(noise).tissuePixels==0,"isolated H pixels rejected");
        Random rnd=new Random(47);for(int i=0;i<na.length;i++){int c=0;for(int j=0;j<3;j++)c=(c<<8)|(200+rnd.nextInt(7)-3);na[i]=c;}
        check(select(noise).tissuePixels==0,"neutral RGB channel noise rejected");
        ImagePlus weak=filled(128,128,0xffffff);int[] wa=(int[])weak.getProcessor().getPixels();
        for(int y=20;y<100;y++)for(int x=20;x<100;x++)wa[y*128+x]=stain(.025,0,0);
        check(select(weak).status.equals("FAILED_NO_RELIABLE_H"),"weak H explicitly fails");
        ImagePlus im=phantom();H_DAB_Tissue_Selector.Selection s=select(im),again=select(im);
        check(Arrays.equals(s.mask,again.mask)&&s.hThreshold==again.hThreshold,"fixed-seed exact reproducibility");
        check(s.mask[30*256],"tissue touching image edge retained");
        check(s.mask[150*256+230],"disconnected peripheral tissue retained");
        check(!s.mask[96*256+128],"empty center does not seed tissue");
        check(!s.mask[88*256+52],"large internal lumen remains empty");
        check(s.mask[40*256+24]&&s.seeds[40*256+24]==0,"pale inter-nuclear tissue included");
        check(s.tissuePixels>4*s.seedPixels,"mask covers tissue beyond H nuclei");
        int intersection=0,union=0;
        for(int y=0;y<192;y++)for(int x=0;x<256;x++) {
            boolean truth=((x<96&&y>=16&&y<160)||(x>=176&&y>=112))&&!(x>=32&&x<72&&y>=64&&y<112);
            if(truth&&s.mask[y*256+x])intersection++;if(truth||s.mask[y*256+x])union++;
        }
        check(intersection/(double)union>=.95,"synthetic pale-tissue geometry IoU >= 0.95");
        System.out.println("PHANTOM_IOU="+intersection/(double)union+" TISSUE="+s.tissuePixels+" H_SUPPORT="+s.seedPixels);
        boolean translationStable=true;boolean[] reference=null;
        for(int offset=0;offset<16;offset++) {
            ImagePlus shifted=filled(384,320,0xffffff);int[] target=(int[])shifted.getProcessor().getPixels();
            int[] source=(int[])im.getProcessor().getPixels();
            for(int y=0;y<192;y++)System.arraycopy(source,y*256,target,(y+48+offset)*384+48+offset,256);
            H_DAB_Tissue_Selector.Selection shiftedMask=select(shifted);
            boolean[] aligned=new boolean[256*192];
            for(int y=0;y<192;y++)for(int x=0;x<256;x++)aligned[y*256+x]=shiftedMask.mask[(y+48+offset)*384+x+48+offset];
            if(offset==0)reference=aligned;else if(!Arrays.equals(reference,aligned))translationStable=false;
        }
        check(translationStable,"all 16 sampling-grid offsets give identical aligned phantom masks");
        H_DAB_Dominance_Extractor.Params changed=params();changed.tissueSeed=27;
        check(Arrays.equals(s.mask,H_DAB_Tissue_Selector.select(im,changed).mask),"representative phantom robust to sampling seed");
        im.setRoi(new OvalRoi(0,16,96,144));s=select(im);boolean outside=false;
        for(int y=0;y<192;y++)for(int x=0;x<256;x++)outside|=s.mask[y*256+x]&&!im.getRoi().contains(x,y);
        check(!outside&&s.tissuePixels>0,"nonrectangular area ROI clips mask");
        im.setRoi(new Roi(112,16,32,32));check(select(im).tissuePixels==0,"background-only ROI does not borrow external nuclei");
        H_DAB_Dominance_Extractor.Result empty=H_DAB_Dominance_Extractor.analyze(filled(32,32,0xffffff),params());
        check(empty!=null&&Double.isNaN(empty.measurements.getValue("dab_area_fraction",0)),"empty denominator exports NaN with status");
        im=phantom();H_DAB_Dominance_Extractor.Params p=params();H_DAB_Dominance_Extractor.Result r=H_DAB_Dominance_Extractor.analyze(im,p);
        byte[] tm=(byte[])r.tissueMaskImage.getProcessor().getPixels(),dm=(byte[])r.maskImage.getProcessor().getPixels();
        boolean subset=true,policy=true;int tc=0,dc=0;double sum=0;
        float[] hc=(float[])r.hemaImage.getProcessor().getPixels(),dd=(float[])r.dabImage.getProcessor().getPixels(),rr=(float[])r.residualImage.getProcessor().getPixels();
        double hs=r.measurements.getValue("hema_scale_p99",0),ds=r.measurements.getValue("dab_scale_p99",0),rs=r.measurements.getValue("residual_scale_p99",0);
        for(int i=0;i<tm.length;i++) {
            if(tm[i]!=0)tc++;if(dm[i]!=0){dc++;sum+=dd[i];if(tm[i]==0)subset=false;}
            double hn=Math.min(1.5,hc[i]/hs),dn=Math.min(1.5,dd[i]/ds),rn=Math.min(1.5,rr[i]/rs);
            boolean expected=tm[i]!=0&&dn>=p.dabMin&&dn-p.alpha*hn>=p.dominanceMin&&dn>=p.ratioMin*hn&&rn<=p.residualMax;
            if(expected!=(dm[i]!=0))policy=false;
        }
        check(subset&&tc==r.measurements.getValue("tissue_pixels",0),"saved mask equals measurement denominator; DAB is subset");
        check(dc>0&&dc<tc,"DAB regression fixture exercises both positive and negative decisions");
        check(policy,"unchanged strict DAB predicate applied inside new tissue mask");
        check(Math.abs(sum-r.measurements.getValue("integrated_dab_concentration",0))<1e-3,"raw DAB integral agrees with selected positives");
        check(Math.abs(dc/(double)tc-r.measurements.getValue("dab_area_fraction",0))<1e-12,"DAB fraction uses selected tissue");
        H_DAB_Dominance_Extractor.Params invalid=params();invalid.tissueMinOD=0;boolean thrown=false;
        try{H_DAB_Tissue_Selector.select(im,invalid);}catch(IllegalArgumentException e){thrown=true;}
        check(thrown,"zero tissue threshold rejected");
        for(String option:new String[]{"tissuereachpixels=NaN","tissuesupportradius=2.5","tissuehmin=abc"}) {
            boolean bad=false;try{H_DAB_Dominance_Extractor.Params.fromOptions(option);}catch(IllegalArgumentException e){bad=true;}
            check(bad,"malformed option rejected: "+option);
        }
        H_DAB_Dominance_Extractor.Params parsed=H_DAB_Dominance_Extractor.Params.fromOptions("tissueblocksize=24 tissuehmin=0.12 tissueautowhite=false tissueseed=42");
        check(parsed.tissueBlockSize==24&&parsed.tissueHMin==.12&&!parsed.tissueAutoWhite&&parsed.tissueSeed==42,"macro parameter parsing");
        System.out.println("ALL "+checks+" CHECKS PASSED");
    }
}
