import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.util.ArrayList;

public class SimulationRunY4 {

    private ArrayList<Yeast4> y1;
    private ArrayList<Yeast4> y2;

    private ArrayList<Yeast4> yFused;

    private int n0_1;
    private int n0_2;

    private int type1;
    private int type2;

    private int L;

    private String label[] = new String[2];

    MersenneTwister rd;

    ImagePlus imp;
    ImageProcessor ip;

    private double[][] cost;

    private double[] n1;
    private double[] n2;
    private double[] nf;

    ArrayList<double[]>[] matingDistances;


    SimulationRunY4(int tp1, int tp2){
        rd = new MersenneTwister();
        y1 = new ArrayList<Yeast4>();
        y2 = new ArrayList<Yeast4>();
        yFused = new ArrayList<Yeast4>();

        label[Yeast4.LOCAL]="L";
        label[Yeast4.GLOBAL]="G";

        cost = new double[2][4];

        n1 = new double[4];
        n2 = new double[4];
        nf = new double[4];

        matingDistances = new ArrayList[4];
        for(int j=0;j<4;j++){
            matingDistances[j]=new ArrayList<double[]>();
        }

        type1 = tp1;  // if tp1 is GLOBAL, species 1 is a global !!sensor!!
        type2 = tp2;

    }

    public void run(double density, double ratio, int size, boolean print, boolean diploDivide){

        y1 = new ArrayList<Yeast4>();
        y2 = new ArrayList<Yeast4>();

        L = size;

        n0_1 = (int)(ratio * density *L*L*4/Math.PI);
        n0_2 = (int)((1-ratio) * density *L*L*4/Math.PI);

        initiate(y1,n0_1,type1,type2);
        initiate(y2,n0_2,type2,type1);
        equilibrate();
        computeInitialDistances();

        if(print) {
            imp = IJ.createImage("Simu_d-"+density+"_r-"+ratio+"_Y1-b-"+label[type1]+"_Y2-r-"+label[type2], "RGB", 4*L, 4*L, 5);
            ip = imp.getProcessor();

        }

        for(int k=0;k<4; k++){

            computeStimulation();

            if(print){
                write(k+1);
            }

            cost[0][k] = computeCost(y1); // compute the number of stimulated before division
            cost[1][k] = computeCost(y2);

            nf[k]-= yFused.size(); // differential at each time step

            n1[k] = y1.size(); // compute the number of cells before division
            n2[k] = y2.size();

            computePutativeYeasts(); // it aims

            if(diploDivide){
                diploDivisionRound(yFused);
            }

            fusionDivisionRound(k);

            nf[k] = yFused.size(); // instantaneous number of diploids

        }

        if(print){ // final stage after t=4
            write(5);
        }

    }

    public double[][] getCost(){
        return cost;
    }
    public double[] getN1(){
        return n1;
    }
    public double[] getN2(){
        return n2;
    }
    public double[] getNF(){
        return nf;
    }
    public ArrayList<double[]>[] getMatingDistances(){
        return matingDistances;
    }

    private void initiate(ArrayList<Yeast4> y, int n, int typeSense,int typeEmit ){

          for(int i=0; i<n; i++){
               Yeast4 yy = new Yeast4(rd,L,typeSense,typeEmit);
               y.add(yy);
          }

    }

    private void diploDivisionRound(ArrayList<Yeast4> y){

        int N = y.size();

        for(int i = 0; i<N; i++){

            y.add(y.get(i).copyYeast());

        }

    }

    private double[][] getPositions(ArrayList<Yeast4> y){

        double[][] out = new double[y.size()][2];

        for(int i = 0; i<y.size(); i++){
            out[i] = y.get(i).getPosition();
        }

        return out;
    }

    private void computeStimulation(){

        for(Yeast4 yy:y1) {
            yy.resetStim( rd, y2, y1.size());
        }
        for(Yeast4 yy:y2) {
            yy.resetStim( rd, y1, y2.size());
        }
        for(Yeast4 yy:y1) {
            yy.setDL();
        }
        for(Yeast4 yy:y2) {
            yy.setDL();
        }

    }

    private double computeCost(ArrayList<Yeast4> y){
        double out=0;
        for(Yeast4 yy:y){
            out+=yy.getCost();
        }
        return out;
    }

    private void computePutativeYeasts(){

        for(Yeast4 yc:y1){
            yc.startDaughter(rd,y2);
        }
        for(Yeast4 yc:y2){
            yc.startDaughter(rd,y1);
        }

    }

    private void fusionDivisionRound(int time){

        double[][] pos2 = getPositions(y2);

        ArrayList<Integer>[][] twosInProx = new ArrayList[L][L];
        for(int i=0;i<L;i++){
            for(int j=0; j<L; j++){
                twosInProx[i][j] = new ArrayList<Integer>();
            }
        }

        for(int i=0;i<y2.size();i++){

            int xx = (int) pos2[i][0];
            int yy = (int) pos2[i][1];

            for (int xxx = -3; xxx < 4; xxx++) {
                for (int yyy = -3; yyy < 4; yyy++) {

                    int xl = (xx + xxx + L) % L;
                    int yl = (yy + yyy + L) % L;

                    twosInProx[xl][yl].add(i);

                }
            }

        }

        boolean[] hit1 = new boolean[y1.size()];
        boolean[] hit2 = new boolean[y2.size()];

        boolean[] mated2 = new boolean[y2.size()];
        boolean[] mated1 = new boolean[y1.size()];

        double nxtX,nxtY,xb,yb,xx,yy;
        double d2_m1m2,d2_b1m2,d2_m1b2,d2_b1b2;

        ArrayList<Distance> distances = new ArrayList<Distance>();

        for(int i1=0;i1<y1.size();i1++){

            xx = y1.get(i1).getX();
            yy = y1.get(i1).getY();

            // future position of the bud
            xb = y1.get(i1).getDaughter().getX();
            yb = y1.get(i1).getDaughter().getY();

            ArrayList<Integer> ngb2 = twosInProx[((int)xx+L)%L][((int)yy+L)%L];

            for(int i2:ngb2){

                double r2[] = pos2[i2];

                nxtX = y2.get(i2).getDaughter().getX();
                nxtY = y2.get(i2).getDaughter().getY();

                d2_m1m2 = sqrDbc( (r2[0] - xx) , (r2[1] - yy) ); // check mother1 mother2 contact
                d2_b1m2 = sqrDbc( (r2[0] - xb) , (r2[1] - yb) ); // check bud1 mother2 contact
                d2_m1b2 = sqrDbc( (nxtX - xx) , (nxtY - yy) );     // check mother1 bud2 contact
                d2_b1b2 = sqrDbc( (nxtX - xb) , (nxtY - yb) );     // check bud1 bud2 contact

                boolean found = d2_m1m2<=1.0 || d2_b1m2<=1.0 || d2_m1b2<=1.0 || d2_b1b2<=1.0;

                if(found){

                    hit1[i1] = true;
                    hit2[i2] = true;
                    distances.add(new Distance(d2_m1m2,i1,i2));

                }

            }
        }

        while(!distances.isEmpty()){
            Distance d0 = distances.get(0);
            boolean test = (!hit1[d0.getI1()] || !hit2[d0.getI2()]);
            while(test){
                distances.remove(0);
                test = !distances.isEmpty();
                if(test) {
                    d0 = distances.get(0);
                    test = (!hit1[d0.getI1()] || !hit2[d0.getI2()]);
                }
            }
            if(!distances.isEmpty()) {

                for(int i=distances.size()-1;i>=0;i--) {
                    Distance d = distances.get(i);
                    if (!hit1[d.getI1()] || !hit2[d.getI2()]) {
                        distances.remove(i);
                    } else if (d.isSmaller(d0.getD())) {
                        d0 = d;
                    }
                }

                hit1[d0.getI1()] = false;
                hit2[d0.getI2()] = false;

                mated1[d0.getI1()] = true;
                mated2[d0.getI2()] = true;

                double[] nD = new double[3];

                double dxo = y1.get(d0.getI1()).kins.get(0).getX() - y2.get(d0.getI2()).kins.get(0).getX() ;
                double dyo = y1.get(d0.getI1()).kins.get(0).getY() - y2.get(d0.getI2()).kins.get(0).getY() ;

                nD[0]=Math.sqrt(sqrDbc( dxo,  dyo ));
                nD[1]=y1.get(d0.getI1()).getFirstClosestD();
                nD[2]=y2.get(d0.getI2()).getFirstClosestD();
                matingDistances[time].add(nD);

                distances.remove(d0);
            }

        }

        fuse(y2,mated2);
        fuse(y1,mated1);

    }

    private void fuse(ArrayList<Yeast4> y, boolean[] mated){
        for(int i=y.size()-1; i>=0; i--){
            if(mated[i]){
                yFused.add(y.remove(i));
            }
            else {

                y.add(y.get(i).getDaughter());
            }

        }
    }

    private void write(int i){

        imp.setSlice(i);
        double xl,yl;
        int color;


        for(int k=0; k<y1.size(); k++){

            xl = 4*y1.get(k).getX();
            yl = 4*y1.get(k).getY();

            color = (int)(255 - 105*(y1.get(k).getCost()));

            for(int xx = (int)xl-2; xx<=(int)xl+2; xx++){
                for(int yy = (int)yl-2; yy<=(int)yl+2; yy++){

                    if((xx-xl)*(xx-xl) + (yy-yl)*(yy-yl) <=2){
                        ip.putPixel((xx+4*L)%(4*L),(yy+4*L)%(4*L),new int[]{0,0,color});
                    }

                }
            }

        }

        for(int k=0; k<y2.size(); k++){

            xl = 4*y2.get(k).getX();
            yl = 4*y2.get(k).getY();

            color = (int)(255 - 105*(y2.get(k).getCost()));

            for(int xx = (int)xl-2; xx<=(int)xl+2; xx++){
                for(int yy = (int)yl-2; yy<=(int)yl+2; yy++){

                    if((xx-xl)*(xx-xl) + (yy-yl)*(yy-yl) <=2){
                        ip.putPixel((xx+4*L)%(4*L),(yy+4*L)%(4*L),new int[]{color,0,0});
                    }

                }
            }
        }

        for(int k=0; k<yFused.size(); k++){

            xl = 4*yFused.get(k).getX();
            yl = 4*yFused.get(k).getY();

            for(int xx = (int)xl-2; xx<=(int)xl+2; xx++){
                for(int yy = (int)yl-2; yy<=(int)yl+2; yy++){

                    if((xx-xl)*(xx-xl) + (yy-yl)*(yy-yl) <=2){
                        ip.putPixel((xx+4*L)%(4*L),(yy+4*L)%(4*L),new int[]{0,255,0});
                    }

                }
            }
        }

    }

    public void saveImage(String dir, int n){

        IJ.save(imp,dir+imp.getTitle()+"_minDL-"+(1+Yeast4.getMinDL())+"_"+"_maxDL-"+(1+Yeast4.getMaxDL())+"_"+n+".tif");

    }

    private void equilibrate(){
        int n_1 = y1.size(), n_2 = y2.size();
        double[][] pos = new double[n_1+n_2][2];
        int nt = n_1+n_2;

        for(int i = 0; i<n_1; i++){
            pos[i] = y1.get(i).getPosition();

        }
        for(int i = 0; i<n_2; i++){
            pos[n_1+i] = y2.get(i).getPosition();
        }

        double[][] force = new double[nt][2];
        int cnt=0;
        while(computeForce(pos,force)>=1e-4*(nt) && cnt<10000){
            for(int i = 0; i<nt; i++){
                pos[i][0] = pos[i][0]+force[i][0];
                pos[i][1] = pos[i][1]+force[i][1];
            }
            cnt++;
        }

        for(int i = 0; i<n_1; i++){
            y1.get(i).setPosition(pos[i]);

        }
        for(int i = 0; i<n_2; i++){
            y2.get(i).setPosition(pos[n_1+i]);
        }

    }

    private double computeForce(double[][] pos, double[][] force){
        double sum = 0;
        double d,sqrtd,th;
        int n = pos.length;
        double[] r,r2;
        for(int i=0; i<n; i++){
            force[i][0]=0;
            force[i][1]=0;
            r=pos[i];
            for(int j=0; j<i; j++){
                r2 = pos[j];
                d=(r[0]-r2[0])*(r[0]-r2[0]) + (r[1]-r2[1])*(r[1]-r2[1]);
                if(d<1 && d!=0.0){
                    sqrtd = Math.sqrt(d);
                    force[i][0] += 0.1*(1-d) * (r[0]-r2[0])/sqrtd;
                    force[i][1] += 0.1*(1-d) * (r[1]-r2[1])/sqrtd;
                    force[j][0] += 0.1*(1-d) * (r2[0]-r[0])/sqrtd;
                    force[j][1] += 0.1*(1-d) * (r2[1]-r[1])/sqrtd;
                    sum+= 0.1*(1-d);
                }
                else if(d==0.0){
                    th = 2*Math.PI*rd.nextDouble();
                    d = 0.1*Math.cos(th);
                    force[i][0] += d;
                    force[j][0] -= d;
                    d = 0.1*Math.sin(th);
                    force[i][1] += d;
                    force[j][1] -= d;
                    sum+= 0.1;
                }

            }

        }

        return sum;
    }

    private void computeInitialDistances(){

        double xx,yy;
        double d;

        for(Yeast4 yy2:y2){
            yy2.setFirstClosestD(2*L);
        }

        for(Yeast4 yy1:y1){

            xx = yy1.getX();
            yy = yy1.getY();

            yy1.setFirstClosestD(2*L);

            for(Yeast4 yy2:y2){

                double r2[] = yy2.getPosition();

                d = Math.sqrt( sqrDbc(r2[0] - xx,r2[1] - yy) ); // check mother1 mother2 contact

                if(d<yy1.getFirstClosestD()){
                    yy1.setFirstClosestD(d);
                }
                if(d<yy2.getFirstClosestD()){
                    yy2.setFirstClosestD(d);
                }

            }
        }

        for(Yeast4 yy1:y1){
            yy1.setStim(rd,yy1.getFirstClosestD()/Yeast4.getMinDL(),(double)y2.size()/(y1.size()+y2.size()));
        }
        for(Yeast4 yy2:y2){
            yy2.setStim(rd,yy2.getFirstClosestD()/Yeast4.getMinDL(),(double)y1.size()/(y1.size()+y2.size()));
        }


    }

    private double sqrDbc(double dx, double dy){
        double dxbc = bc(dx);
        double dybc = bc(dy);
        return dxbc*dxbc + dybc * dybc;
    }

    private double bc(double x){ // return between -L/2 and L/2
        return x - Math.round(x/L)*L;
    }
}
