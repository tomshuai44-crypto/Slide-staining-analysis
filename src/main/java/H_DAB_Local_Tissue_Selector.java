import ij.ImagePlus;
import ij.gui.Roi;
import java.util.Arrays;
import java.util.Random;

/** Native Java implementation of the reviewed local-transition experiment.
 * Both diagnostic stains use H/D/neutral OD. Raw measurement channels remain
 * in the extractor's original H/D/residual basis. No maximum growth radius.
 */
public final class H_DAB_Local_Tissue_Selector {
    public static final String VERSION = "local-hdab-1.0";
    private static final int STRIDE = 4;
    private static final int BACKGROUND_LIMIT = 100000;
    private static final double[] H = unit(.650,.704,.286), D = unit(.268,.570,.776);
    private static final double DET = (H[0]-H[2])*(D[1]-D[2])-(D[0]-D[2])*(H[1]-H[2]);
    private static final double HR=(D[1]-D[2])/DET, HG=-(D[0]-D[2])/DET;
    private static final double DR=-(H[1]-H[2])/DET, DG=(H[0]-H[2])/DET;

    public static H_DAB_Tissue_Selector.Selection select(ImagePlus imp, H_DAB_Dominance_Extractor.Params p) {
        validate(p);
        H_DAB_Tissue_Selector.Selection s=H_DAB_Tissue_Selector.selectSeeds(imp,p);
        s.algorithm=VERSION;
        s.diagnostics.put("k",p.tissueLocalK);
        s.diagnostics.put("material_min_od",p.tissueLocalMinOD);
        s.diagnostics.put("feature_stride_px",(double)STRIDE);
        s.diagnostics.put("transition_window_px",5.0);
        s.diagnostics.put("transition_fraction",.60);
        if(!"SEEDS_READY".equals(s.status))return s;
        int w=imp.getWidth(),h=imp.getHeight(),n=w*h,radius=p.tissueSupportRadius;
        int[] pixels=(int[])imp.getProcessor().convertToRGB().getPixels();
        Roi roi=imp.getRoi();
        boolean[] domain=new boolean[n],notDark=new boolean[n];
        float[] weights=roi==null?null:new float[n];
        float[][] values={new float[n],new float[n],new float[n]};
        float[] lum=new float[n],lumSquared=new float[n];
        float[] white={(float)s.whiteR,(float)s.whiteG,(float)s.whiteB};
        float[][] od=new float[3][256];
        for(int c=0;c<3;c++)for(int v=0;v<256;v++)od[c][v]=(float)-Math.log((v+1f)/(white[c]+1f));
        for(int y=0;y<h;y++)for(int x=0;x<w;x++) {
            int i=y*w+x;
            if(roi!=null&&!roi.contains(x,y))continue;
            domain[i]=true;if(weights!=null)weights[i]=1;
            int c=pixels[i],rr=(c>>16)&255,gg=(c>>8)&255,bb=c&255;
            notDark[i]=Math.min(rr,Math.min(gg,bb))>p.tissueDarkMax;
            float a=od[0][rr],b=od[1][gg],d=od[2][bb];
            values[0][i]=(float)Math.max(0,HR*((double)a-d)+HG*((double)b-d));
            values[1][i]=(float)Math.max(0,DR*((double)a-d)+DG*((double)b-d));
            a=Math.max(0,a);b=Math.max(0,b);d=Math.max(0,d);
            values[2][i]=(float)Math.sqrt(a*a+b*b+d*d);
            lum[i]=(rr+gg+bb)/3f;lumSquared[i]=lum[i]*lum[i];
        }
        float[] windowWeights=weights==null?null:uniform(weights,w,h,radius,true);
        float[] mean=maskedMean(lum,windowWeights,w,h,radius);
        float[] meanSquared=maskedMean(lumSquared,windowWeights,w,h,radius);
        float[] sd=new float[n];
        for(int i=0;i<n;i++)sd[i]=(float)Math.sqrt(Math.max(0,meanSquared[i]-mean[i]*mean[i]));
        double bright=quantile(collect(mean,domain),.90);
        boolean[] brightMask=new boolean[n];
        for(int i=0;i<n;i++)brightMask[i]=domain[i]&&mean[i]>=bright;
        double flat=quantile(collect(sd,brightMask),.25);
        int[] background=new int[n];int count=0;
        for(int i=0;i<n;i++)if(brightMask[i]&&sd[i]<=flat&&notDark[i])background[count++]=i;
        s.diagnostics.put("background_pool_pixels",(double)count);
        s.diagnostics.put("brightness_q90",bright);
        s.diagnostics.put("flat_sd_q25",flat);
        if(count<64){s.status="FAILED_LOCAL_BACKGROUND";return s;}
        // Java's fixed-seed partial shuffle replaces NumPy PCG64 choice. The
        // distribution rule is the same; native-port agreement is validated.
        Random random=new Random(p.tissueSeed);
        int sampleCount=Math.min(BACKGROUND_LIMIT,count);
        for(int i=0;i<sampleCount;i++) {
            int j=i+random.nextInt(count-i),swap=background[i];background[i]=background[j];background[j]=swap;
        }
        background=Arrays.copyOf(background,sampleCount);
        s.diagnostics.put("background_samples",(double)sampleCount);
        boolean[][] evidence=new boolean[3][];
        boolean[] texture=new boolean[n];
        String[] names={"h","d","v"};
        for(int stain=0;stain<3;stain++) {
            float[] a=values[stain];
            float[][] q=localQuantiles(a,domain,w,h,radius);
            float[] spread=new float[n];
            for(int i=0;i<n;i++)spread[i]=q[2][i]-q[0][i];
            double excess=Math.max(.04,quantile(sample(a,background),.99));
            float[] indicator=new float[n];
            for(int i=0;i<n;i++)indicator[i]=domain[i]&&a[i]>excess?1:0;
            float[] fraction=maskedMean(indicator,windowWeights,w,h,radius);
            String name=names[stain];
            double tm=threshold(q[1],background,.01,p.tissueLocalK,s,name+"_median");
            double ts=threshold(spread,background,.015,p.tissueLocalK,s,name+"_spread");
            double tf=threshold(fraction,background,.02,p.tissueLocalK,s,name+"_fraction");
            s.diagnostics.put(name+"_pixel_excess_cut",excess);
            evidence[stain]=new boolean[n];
            for(int i=0;i<n;i++) {
                evidence[stain][i]=q[1][i]>tm||(spread[i]>ts&&fraction[i]>tf);
                if(stain<2&&spread[i]>ts)texture[i]=true;
            }
        }
        boolean[] seedMask=new boolean[n];
        for(int i=0;i<n;i++)seedMask[i]=s.seeds[i]!=0;
        boolean[] nearbySeeds=morph(seedMask,w,h,radius,true);
        float[] local=new float[n];
        for(int i=0;i<n;i++)if(domain[i]&&(evidence[0][i]||evidence[1][i]||(evidence[2][i]&&texture[i])||nearbySeeds[i]))local[i]=1;
        float[] sustained=uniform(local,w,h,2,false);
        // ROI windows use only the domain; image edges retain the experiment's
        // constant-zero support convention. Seeds are always preserved.
        if(weights!=null) {
            float[] supportWeights=uniform(weights,w,h,2,false);
            for(int i=0;i<n;i++)if(supportWeights[i]>0)sustained[i]/=supportWeights[i];
        }
        boolean[] material=new boolean[n];
        for(int i=0;i<n;i++)material[i]=domain[i]&&notDark[i]&&values[2][i]>=p.tissueLocalMinOD;
        boolean[] closed=morph(morph(material,w,h,p.tissueCloseRadius,true),w,h,p.tissueCloseRadius,false);
        boolean[] allowed=new boolean[n];
        for(int i=0;i<n;i++)allowed[i]=domain[i]&&((sustained[i]>=.6f&&(material[i]||closed[i])&&notDark[i])||seedMask[i]);
        // Multi-source flood fill is exactly component reconstruction, not a
        // limited-distance dilation. Clear gaps and the ROI boundary block paths.
        int[] queue=new int[n];int tail=0;
        for(int i=0;i<n;i++)if(seedMask[i]){s.mask[i]=true;queue[tail++]=i;}
        for(int head=0;head<tail;head++) {
            int i=queue[head],x=i%w,y=i/w;
            for(int yy=Math.max(0,y-1);yy<=Math.min(h-1,y+1);yy++)for(int xx=Math.max(0,x-1);xx<=Math.min(w-1,x+1);xx++) {
                int j=yy*w+xx;if(allowed[j]&&!s.mask[j]){s.mask[j]=true;queue[tail++]=j;}
            }
        }
        s.tissuePixels=tail;
        s.status=tail>=10?"REVIEW_REQUIRED":"FAILED_INSUFFICIENT_TISSUE";
        if(tail<10){Arrays.fill(s.mask,false);s.tissuePixels=0;}
        return s;
    }

    private static double threshold(float[] feature,int[] background,double floor,double k,H_DAB_Tissue_Selector.Selection s,String key) {
        float[] sample=sample(feature,background);
        double center=quantile(sample,.5);
        float[] deviation=new float[sample.length];
        for(int i=0;i<sample.length;i++)deviation[i]=(float)Math.abs(sample[i]-center);
        double scale=Math.max(floor,1.4826*quantile(deviation,.5));
        double cutoff=center+k*scale;
        s.diagnostics.put(key+"_bg_median",center);
        s.diagnostics.put(key+"_bg_scale",scale);
        s.diagnostics.put(key+"_threshold",cutoff);
        return cutoff;
    }

    private static float[][] localQuantiles(float[] a,boolean[] domain,int w,int h,int r) {
        int nw=(w+STRIDE-1)/STRIDE,nh=(h+STRIDE-1)/STRIDE;
        float[][] grid=new float[3][nw*nh];
        float[] scratch=new float[(2*r+1)*(2*r+1)];
        for(int gy=0;gy<nh;gy++)for(int gx=0;gx<nw;gx++) {
            int count=0;
            for(int dy=-r;dy<=r;dy++)for(int dx=-r;dx<=r;dx++) {
                int i=reflect(gy*STRIDE+dy,h,false)*w+reflect(gx*STRIDE+dx,w,false);
                if(domain[i])scratch[count++]=a[i];
            }
            if(count==0)continue;
            Arrays.sort(scratch,0,count);
            for(int q=0;q<3;q++)grid[q][gy*nw+gx]=(float)sortedQuantile(scratch,count,q==0?.1:q==1?.5:.9);
        }
        float[][] result=new float[3][a.length];
        for(int y=0;y<h;y++)for(int x=0;x<w;x++) {
            double yy=Math.min(nh-1,y/(double)STRIDE),xx=Math.min(nw-1,x/(double)STRIDE);
            int y0=(int)yy,x0=(int)xx,y1=Math.min(nh-1,y0+1),x1=Math.min(nw-1,x0+1);
            double fy=yy-y0,fx=xx-x0;
            for(int q=0;q<3;q++) {
                float[] g=grid[q];
                result[q][y*w+x]=(float)((1-fy)*((1-fx)*g[y0*nw+x0]+fx*g[y0*nw+x1])+fy*((1-fx)*g[y1*nw+x0]+fx*g[y1*nw+x1]));
            }
        }
        return result;
    }

    private static float[] maskedMean(float[] a,float[] windowWeights,int w,int h,int r) {
        float[] out=uniform(a,w,h,r,true);
        if(windowWeights!=null)for(int i=0;i<out.length;i++)out[i]=windowWeights[i]>0?out[i]/windowWeights[i]:0;
        return out;
    }

    /** Separable uniform filter, vertical then horizontal, rounded to float at
     * each pass like SciPy's float32 uniform_filter. */
    private static float[] uniform(float[] a,int w,int h,int r,boolean reflected) {
        float[] tmp=new float[a.length],out=new float[a.length];int size=2*r+1;
        for(int x=0;x<w;x++) {
            double sum=0;for(int dy=-r;dy<=r;dy++)sum+=get(a,w,h,x,dy,reflected);
            for(int y=0;y<h;y++) {
                tmp[y*w+x]=(float)(sum/size);
                sum-=get(a,w,h,x,y-r,reflected);sum+=get(a,w,h,x,y+r+1,reflected);
            }
        }
        for(int y=0;y<h;y++) {
            double sum=0;for(int dx=-r;dx<=r;dx++)sum+=get(tmp,w,h,dx,y,reflected);
            for(int x=0;x<w;x++) {
                out[y*w+x]=(float)(sum/size);
                sum-=get(tmp,w,h,x-r,y,reflected);sum+=get(tmp,w,h,x+r+1,y,reflected);
            }
        }
        return out;
    }

    private static float get(float[] a,int w,int h,int x,int y,boolean reflected) {
        if(x>=0&&x<w&&y>=0&&y<h)return a[y*w+x];
        return reflected?a[reflect(y,h,true)*w+reflect(x,w,true)]:0;
    }

    // repeatEdge=true matches ndimage reflect; false matches numpy.pad reflect.
    private static int reflect(int p,int length,boolean repeatEdge) {
        if(length==1)return 0;
        int period=repeatEdge?2*length:2*(length-1);
        int q=p%period;if(q<0)q+=period;
        return q<length?q:(repeatEdge?period-1-q:period-q);
    }

    private static boolean[] morph(boolean[] a,int w,int h,int r,boolean dilation) {
        if(r==0)return a.clone();
        // Boolean min/max via integer horizontal then vertical counts. Erosion
        // treats outside-image pixels as false, matching the experiment.
        boolean[] tmp=new boolean[a.length],out=new boolean[a.length];int size=2*r+1;
        for(int y=0;y<h;y++) {
            int sum=0;for(int x=0;x<=Math.min(w-1,r);x++)if(a[y*w+x])sum++;
            for(int x=0;x<w;x++) {
                tmp[y*w+x]=dilation?sum>0:sum==size;
                if(x-r>=0&&a[y*w+x-r])sum--;if(x+r+1<w&&a[y*w+x+r+1])sum++;
            }
        }
        for(int x=0;x<w;x++) {
            int sum=0;for(int y=0;y<=Math.min(h-1,r);y++)if(tmp[y*w+x])sum++;
            for(int y=0;y<h;y++) {
                out[y*w+x]=dilation?sum>0:sum==size;
                if(y-r>=0&&tmp[(y-r)*w+x])sum--;if(y+r+1<h&&tmp[(y+r+1)*w+x])sum++;
            }
        }
        return out;
    }

    private static float[] collect(float[] a,boolean[] mask) {
        int n=0;for(boolean b:mask)if(b)n++;
        float[] out=new float[n];int j=0;for(int i=0;i<a.length;i++)if(mask[i])out[j++]=a[i];return out;
    }
    private static float[] sample(float[] a,int[] indexes) {
        float[] out=new float[indexes.length];for(int i=0;i<indexes.length;i++)out[i]=a[indexes[i]];return out;
    }
    private static double quantile(float[] a,double q) {Arrays.sort(a);return sortedQuantile(a,a.length,q);}
    private static double sortedQuantile(float[] a,int n,double q) {
        if(n==0)return 0;double k=q*(n-1);int lo=(int)k,hi=(int)Math.ceil(k);return a[lo]+(k-lo)*((double)a[hi]-a[lo]);
    }
    private static double[] unit(double a,double b,double c) {double n=Math.sqrt(a*a+b*b+c*c);return new double[]{a/n,b/n,c/n};}
    public static void validate(H_DAB_Dominance_Extractor.Params p) {
        if(!Double.isFinite(p.tissueLocalK)||p.tissueLocalK<=0||p.tissueLocalK>100
                ||!Double.isFinite(p.tissueLocalMinOD)||p.tissueLocalMinOD<=0
                ||p.tissueSupportRadius>32)
            throw new IllegalArgumentException("Local tissue selection requires k in (0,100], positive material OD, and support radius 1..32 pixels.");
    }
}
