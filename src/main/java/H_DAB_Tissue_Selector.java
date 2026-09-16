import ij.ImagePlus;
import ij.gui.Roi;
import java.util.Arrays;
import java.util.Random;
import java.util.LinkedHashMap;

/** H-guided tissue support, independent of the downstream DAB-positive gate. */
public final class H_DAB_Tissue_Selector {
    public static final String VERSION = "h-spatial-1.2";
    private static final double[] H = unit(.650, .704, .286);
    private static final double[] D = unit(.268, .570, .776);
    // Solve OD differences against H and D. A neutral OD offset cancels exactly.
    private static final double DET = (H[0]-H[2])*(D[1]-D[2])-(D[0]-D[2])*(H[1]-H[2]);
    private static final double HR = (D[1]-D[2])/DET;
    private static final double HG = -(D[0]-D[2])/DET;

    public static class Selection {
        public boolean[] mask;
        public byte[] seeds;
        public int domainPixels, tissuePixels, seedPixels, seedBlocks, sampledBlocks;
        public double whiteR, whiteG, whiteB, hThreshold, backgroundHSigma;
        public String status, referenceStatus;
        public String algorithm = VERSION;
        public final LinkedHashMap<String, Double> diagnostics = new LinkedHashMap<String, Double>();
    }

    public static Selection select(ImagePlus imp, H_DAB_Dominance_Extractor.Params p) {
        return selectInternal(imp, p, false);
    }

    public static Selection selectSeeds(ImagePlus imp, H_DAB_Dominance_Extractor.Params p) {
        return selectInternal(imp, p, true);
    }

    private static Selection selectInternal(ImagePlus imp, H_DAB_Dominance_Extractor.Params p, boolean seedsOnly) {
        validate(p);
        int w=imp.getWidth(), h=imp.getHeight(), n=w*h, bs=p.tissueBlockSize;
        int nx=(w+bs-1)/bs, ny=(h+bs-1)/bs, nb=nx*ny;
        int[] rgb=(int[])imp.getProcessor().convertToRGB().getPixels();
        Roi roi=imp.getRoi();
        if (roi!=null && !roi.isArea()) throw new IllegalArgumentException("Tissue selection requires an area ROI, not a line or point ROI.");
        boolean[] domain=new boolean[n];
        Selection s=new Selection(); s.mask=new boolean[n]; s.seeds=new byte[n];
        for(int y=0;y<h;y++) for(int x=0;x<w;x++) if(roi==null || roi.contains(x,y)) {domain[y*w+x]=true;s.domainPixels++;}
        if(s.domainPixels==0) {s.status="FAILED_EMPTY_ROI";s.referenceStatus="unavailable";return s;}

        // Sample the center first, then one random block in every spatial stratum.
        // Every block/pixel is subsequently scanned, so sampling never determines coverage.
        Random random=new Random(p.tissueSeed);
        int sx=Math.min(nx,64), sy=Math.min(ny,64);
        int[] blocks=new int[sx*sy+1]; int bc=0;
        blocks[bc++]=(ny/2)*nx+nx/2;
        for(int y=0;y<sy;y++) for(int x=0;x<sx;x++) {
            int x0=x*nx/sx,x1=(x+1)*nx/sx,y0=y*ny/sy,y1=(y+1)*ny/sy;
            int b=(y0+random.nextInt(y1-y0))*nx+x0+random.nextInt(x1-x0);
            if(b!=blocks[0]) blocks[bc++]=b;
        }
        double[][] means=new double[bc][3]; double[] lum=new double[bc],sd=new double[bc];
        int[] samples=new int[bc*32]; int ns=0, validBlocks=0;
        for(int k=0;k<bc;k++) {
            int bx=blocks[k]%nx*bs,by=blocks[k]/nx*bs, bw=Math.min(bs,w-bx),bh=Math.min(bs,h-by);
            int count=0;double sum=0,sum2=0;double[] m=new double[3];
            for(int j=0;j<32;j++) {
                int i=(by+random.nextInt(bh))*w+bx+random.nextInt(bw);
                if(!domain[i])continue;
                samples[ns++]=i;int c=rgb[i]; double r=(c>>16)&255,g=(c>>8)&255,b=c&255;
                m[0]+=r;m[1]+=g;m[2]+=b;double v=(r+g+b)/3;sum+=v;sum2+=v*v;count++;
            }
            if(count==0)continue;
            for(int c=0;c<3;c++) means[validBlocks][c]=m[c]/count;
            lum[validBlocks]=sum/count;sd[validBlocks]=Math.sqrt(Math.max(0,sum2/count-Math.pow(sum/count,2)));validBlocks++;
        }
        s.sampledBlocks=validBlocks;
        double[] white={p.tissueWhiteR,p.tissueWhiteG,p.tissueWhiteB};
        if(p.tissueAutoWhite && validBlocks>0) {
            double bright=quantile(lum,validBlocks,.90);double[] flat=new double[validBlocks];int nf=0;
            for(int k=0;k<validBlocks;k++)if(lum[k]>=bright)flat[nf++]=sd[k];
            double flatLimit=quantile(flat,nf,.25);double[][] refs=new double[3][validBlocks];int nr=0;
            for(int k=0;k<validBlocks;k++)if(lum[k]>=bright && sd[k]<=flatLimit) {
                for(int c=0;c<3;c++)refs[c][nr]=means[k][c];nr++;
            }
            for(int c=0;c<3;c++)white[c]=quantile(refs[c],nr,.5);
            s.referenceStatus=flatLimit<=8 && Math.min(white[0],Math.min(white[1],white[2]))>=128
                    ? "AUTO_ESTIMATED_REVIEW" : "AUTO_UNRELIABLE";
        } else s.referenceStatus=p.tissueAutoWhite ? "AUTO_UNRELIABLE" : "MANUAL_REFERENCE";
        s.whiteR=white[0];s.whiteG=white[1];s.whiteB=white[2];
        if(s.referenceStatus.equals("AUTO_UNRELIABLE")) {s.status="FAILED_REFERENCE";return s;}
        double[][] od=new double[3][256];
        for(int c=0;c<3;c++)for(int v=0;v<256;v++)od[c][v]=-Math.log((v+1)/(white[c]+1));
        // Estimate the density distribution using bounded Monte Carlo samples.
        double[] hs=new double[ns],bg=new double[ns];int nh=0,ng=0;
        for(int j=0;j<ns;j++) {
            int c=rgb[samples[j]],r=(c>>16)&255,g=(c>>8)&255,b=c&255;
            if(Math.min(r,Math.min(g,b))<=p.tissueDarkMax)continue;
            double hr=hSignal(od[0][r],od[1][g],od[2][b]);
            if(hr>=p.tissueHMin)hs[nh++]=hr;
            if(norm(od[0][r],od[1][g],od[2][b])<p.tissueMinOD)bg[ng++]=hr;
        }
        double med=quantile(bg,ng,.5);double[] dev=new double[ng];
        for(int i=0;i<ng;i++)dev[i]=Math.abs(bg[i]-med);
        s.backgroundHSigma=1.4826*quantile(dev,ng,.5);
        s.hThreshold=Math.max(p.tissueHMin, Math.max(.5*quantile(hs,nh,.90), med+4*s.backgroundHSigma));

        boolean[] strong=new boolean[n],candidate=new boolean[n],supported=new boolean[n];
        for(int i=0;i<n;i++) {
            if(!domain[i])continue;
            int c=rgb[i],r=(c>>16)&255,g=(c>>8)&255,b=c&255;
            if(Math.min(r,Math.min(g,b))<=p.tissueDarkMax)continue;
            double a=od[0][r],d=od[1][g],e=od[2][b];
            candidate[i]=norm(a,d,e)>=p.tissueMinOD;
            strong[i]=candidate[i] && hSignal(a,d,e)>=s.hThreshold;
        }
        // Require a local cluster; a single H-coloured pixel cannot seed tissue.
        for(int y=0;y<h;y++)for(int x=0;x<w;x++) {
            int i=y*w+x;if(!strong[i])continue;int count=0;
            for(int yy=Math.max(0,y-1);yy<=Math.min(h-1,y+1);yy++)
                for(int xx=Math.max(0,x-1);xx<=Math.min(w-1,x+1);xx++)if(strong[yy*w+xx])count++;
            supported[i]=count>=3;
        }
        // Sliding neighbourhoods replace hard tile decisions. Sampling tiles now
        // influence only reference/threshold estimation, never spatial boundaries.
        int radius=p.tissueSupportRadius;
        int[] counts=boxCounts(supported,w,h,radius);
        int[] domainCounts=roi==null ? null : boxCounts(domain,w,h,radius);
        boolean[] seededBlocks=new boolean[nb];
        int[] queue=new int[n];int tail=0;
        for(int y=0;y<h;y++)for(int x=0;x<w;x++) {
            int i=y*w+x;
            int localArea=domainCounts==null ? windowArea(x,y,w,h,radius) : domainCounts[i];
            if(supported[i] && counts[i]>=Math.max(p.tissueMinSeedPixels,Math.ceil(localArea*p.tissueSeedFraction))) {
                s.seeds[i]=(byte)255;s.seedPixels++;s.mask[i]=true;queue[tail++]=i;
                int block=(y/bs)*nx+x/bs;
                if(!seededBlocks[block]){seededBlocks[block]=true;s.seedBlocks++;}
            }
        }
        if(tail==0){s.status="FAILED_NO_RELIABLE_H";return s;}
        if(seedsOnly) {
            Arrays.fill(s.mask,false);
            s.status="SEEDS_READY";
            return s;
        }
        counts=boxCounts(candidate,w,h,radius);
        // Closing fills only small inter-pixel gaps. Large clear spaces stay empty.
        boolean[] allowed=close(candidate,domain,w,h,p.tissueCloseRadius);
        for(int y=0;y<h;y++)for(int x=0;x<w;x++) {
            int i=y*w+x,c=rgb[i];
            int localArea=domainCounts==null ? windowArea(x,y,w,h,radius) : domainCounts[i];
            allowed[i]=allowed[i] && domain[i]
                    && Math.min((c>>16)&255,Math.min((c>>8)&255,c&255))>p.tissueDarkMax
                    && (s.mask[i] || counts[i]>=Math.ceil(localArea*p.tissueVisibleFraction));
        }
        // Multi-source bounded chamfer-geodesic expansion. Cardinal steps cost
        // 3 and diagonals 4 (an octagonal approximation to Euclidean distance).
        // Unlike chessboard distance, this does not impose square reach contours.
        int maxDistance=3*p.tissueReachPixels;
        short[] distance=new short[n];Arrays.fill(distance,Short.MAX_VALUE);
        int[][] buckets=new int[maxDistance+1][];int[] sizes=new int[maxDistance+1];
        buckets[0]=Arrays.copyOf(queue,tail);sizes[0]=tail;
        for(int k=0;k<tail;k++)distance[queue[k]]=0;
        s.tissuePixels=tail;
        for(int d=0;d<=maxDistance;d++) {
            for(int k=0;k<sizes[d];k++) {
                int i=buckets[d][k];if(distance[i]!=d)continue;
                int x=i%w,y=i/w;
                for(int yy=Math.max(0,y-1);yy<=Math.min(h-1,y+1);yy++)for(int xx=Math.max(0,x-1);xx<=Math.min(w-1,x+1);xx++) {
                    if(xx==x&&yy==y)continue;
                    int j=yy*w+xx, next=d+((xx==x||yy==y)?3:4);
                    if(allowed[j]&&next<=maxDistance&&next<distance[j]) {
                        distance[j]=(short)next;
                        if(!s.mask[j]){s.mask[j]=true;s.tissuePixels++;}
                        if(buckets[next]==null)buckets[next]=new int[256];
                        if(sizes[next]==buckets[next].length)buckets[next]=Arrays.copyOf(buckets[next],2*sizes[next]);
                        buckets[next][sizes[next]++]=j;
                    }
                }
            }
            buckets[d]=null;
        }
        s.status=s.tissuePixels>=10 ? "REVIEW_REQUIRED" : "FAILED_INSUFFICIENT_TISSUE";
        if(s.tissuePixels<10){Arrays.fill(s.mask,false);s.tissuePixels=0;}
        return s;
    }

    private static int windowArea(int x,int y,int w,int h,int r) {
        return (Math.min(w-1,x+r)-Math.max(0,x-r)+1)*(Math.min(h-1,y+r)-Math.max(0,y-r)+1);
    }
    private static int[] boxCounts(boolean[] a,int w,int h,int r) {
        int[] tmp=new int[a.length],out=new int[a.length];
        for(int y=0;y<h;y++) {
            int sum=0;for(int x=0;x<=Math.min(w-1,r);x++)if(a[y*w+x])sum++;
            for(int x=0;x<w;x++) {
                tmp[y*w+x]=sum;
                if(x-r>=0 && a[y*w+x-r])sum--;if(x+r+1<w && a[y*w+x+r+1])sum++;
            }
        }
        for(int x=0;x<w;x++) {
            int sum=0;for(int y=0;y<=Math.min(h-1,r);y++)sum+=tmp[y*w+x];
            for(int y=0;y<h;y++) {
                out[y*w+x]=sum;
                if(y-r>=0)sum-=tmp[(y-r)*w+x];if(y+r+1<h)sum+=tmp[(y+r+1)*w+x];
            }
        }
        return out;
    }
    private static boolean[] close(boolean[] mask,boolean[] domain,int w,int h,int r) {
        if(r==0)return mask;
        boolean[] out=boxMorph(boxMorph(mask,w,h,r,true),w,h,r,false);
        for(int i=0;i<out.length;i++)out[i]=domain[i]&&(out[i]||mask[i]);
        return out;
    }
    // Separable square morphology, truncated at the image boundary (no edge erosion).
    private static boolean[] boxMorph(boolean[] a,int w,int h,int r,boolean dilation) {
        boolean[] tmp=new boolean[a.length],out=new boolean[a.length];
        for(int y=0;y<h;y++) {
            int sum=0;for(int x=0;x<=Math.min(w-1,r);x++)if(a[y*w+x])sum++;
            for(int x=0;x<w;x++) {
                int count=Math.min(w-1,x+r)-Math.max(0,x-r)+1;tmp[y*w+x]=dilation ? sum>0 : sum==count;
                if(x-r>=0 && a[y*w+x-r])sum--;if(x+r+1<w && a[y*w+x+r+1])sum++;
            }
        }
        for(int x=0;x<w;x++) {
            int sum=0;for(int y=0;y<=Math.min(h-1,r);y++)if(tmp[y*w+x])sum++;
            for(int y=0;y<h;y++) {
                int count=Math.min(h-1,y+r)-Math.max(0,y-r)+1;out[y*w+x]=dilation ? sum>0 : sum==count;
                if(y-r>=0 && tmp[(y-r)*w+x])sum--;if(y+r+1<h && tmp[(y+r+1)*w+x])sum++;
            }
        }
        return out;
    }
    private static double hSignal(double r,double g,double b){return HR*(r-b)+HG*(g-b);}
    private static double norm(double r,double g,double b){r=Math.max(0,r);g=Math.max(0,g);b=Math.max(0,b);return Math.sqrt(r*r+g*g+b*b);}
    private static double[] unit(double a,double b,double c){double n=Math.sqrt(a*a+b*b+c*c);return new double[]{a/n,b/n,c/n};}
    private static double quantile(double[] a,int n,double q){if(n==0)return 0;double[] b=Arrays.copyOf(a,n);Arrays.sort(b);double k=q*(n-1);int lo=(int)k,hi=(int)Math.ceil(k);return b[lo]+(k-lo)*(b[hi]-b[lo]);}
    public static void validate(H_DAB_Dominance_Extractor.Params p) {
        if(!Double.isFinite(p.tissueMinOD)||p.tissueMinOD<=0 || !Double.isFinite(p.tissueHMin)||p.tissueHMin<=0
            || p.tissueBlockSize<2||p.tissueBlockSize>512 || p.tissueReachPixels<0||p.tissueReachPixels>4096
            || p.tissueSupportRadius<1||p.tissueSupportRadius>256
            || p.tissueCloseRadius<0||p.tissueCloseRadius>32 || p.tissueMinSeedPixels<3
            || !Double.isFinite(p.tissueSeedFraction)||p.tissueSeedFraction<=0||p.tissueSeedFraction>1
            || !Double.isFinite(p.tissueVisibleFraction)||p.tissueVisibleFraction<=0||p.tissueVisibleFraction>1
            || p.tissueDarkMax<0||p.tissueDarkMax>127
            || !Double.isFinite(p.tissueWhiteR)||p.tissueWhiteR<128||p.tissueWhiteR>255
            || !Double.isFinite(p.tissueWhiteG)||p.tissueWhiteG<128||p.tissueWhiteG>255
            || !Double.isFinite(p.tissueWhiteB)||p.tissueWhiteB<128||p.tissueWhiteB>255)
            throw new IllegalArgumentException("Invalid tissue parameters: positive OD/H thresholds, sampling block 2..512 px, support radius 1..256 px, reach 0..4096 px, closing 0..32 px, seed pixels >=3, fractions (0,1], dark cutoff 0..127, reference RGB 128..255 required.");
    }
}
