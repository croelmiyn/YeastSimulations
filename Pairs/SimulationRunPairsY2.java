import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.util.ArrayList;

public class SimulationRunPairsY2 {

    private ArrayList<Yeast2> y1;
    private ArrayList<Yeast2> y2;

    private ArrayList<Yeast2> yFused;

    MersenneTwister rd;

    private double[] fusionSuccess;

    SimulationRunPairsY2(){
        rd = new MersenneTwister();
        y1 = new ArrayList<Yeast2>();
        y2 = new ArrayList<Yeast2>();
        yFused = new ArrayList<Yeast2>();
        fusionSuccess = new double[4];
    }

    public void run(double distance, int type1, int type2, boolean globalSwitch, int nRepeats){

        for(int cnt=0;cnt<nRepeats;cnt++) {

            y1 = new ArrayList<Yeast2>();
            y2 = new ArrayList<Yeast2>();
            yFused = new ArrayList<Yeast2>();

            initiate(y1, 0, type1);
            initiate(y2, distance, type2);

            for (int k = 0; k < 4; k++) {

                if(k==2 && globalSwitch){
                    if(type1==Yeast2.GLOBAL)
                        for(Yeast2 yy1:y1)
                            yy1.setDivMode(Yeast2.AXIAL);
                    if(type2==Yeast2.GLOBAL)
                        for(Yeast2 yy2:y2)
                            yy2.setDivMode(Yeast2.AXIAL);
                }

                computePutativeYeasts(); // it aims

                fusionDivisionRound(k);

                fusionSuccess[k] += (yFused.size() > 0) ? 1 : 0;

            }
        }
        for (int k = 0; k < 4; k++) {
            fusionSuccess[k]/=nRepeats;
        }


    }

    public double[] getFusionSuccess(){

        return fusionSuccess;


    }

    private void initiate(ArrayList<Yeast2> y, double x, int type ){

        Yeast2 yy = new Yeast2(rd,0.0,type);
        yy.setX(x);
        y.add(yy);

    }
    private double[][] getPositions(ArrayList<Yeast2> y){

        double[][] out = new double[y.size()][2];

        for(int i = 0; i<y.size(); i++){
            out[i] = y.get(i).getPosition();
        }

        return out;
    }

    private void computePutativeYeasts(){

        for(Yeast2 yc:y1){
            yc.startDaughter(rd,y2);
        }
        for(Yeast2 yc:y2){
            yc.startDaughter(rd,y1);
        }

    }

    private void fusionDivisionRound(int time){

        double[][] pos2 = getPositions(y2);

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

            for(int i2=0;i2<y2.size();i2++){

                double r2[] = pos2[i2];

                nxtX = y2.get(i2).getDaughter().getX();
                nxtY = y2.get(i2).getDaughter().getY();

                d2_m1m2 = sqrD( (r2[0] - xx) , (r2[1] - yy) ); // check mother1 mother2 contact
                d2_b1m2 = sqrD( (r2[0] - xb) , (r2[1] - yb) ); // check bud1 mother2 contact
                d2_m1b2 = sqrD( (nxtX - xx) , (nxtY - yy) );     // check mother1 bud2 contact
                d2_b1b2 = sqrD( (nxtX - xb) , (nxtY - yb) );     // check bud1 bud2 contact

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

                distances.remove(d0);
            }

        }

        fuse(y2,mated2);
        fuse(y1,mated1);

    }

    private void fuse(ArrayList<Yeast2> y, boolean[] mated){
        for(int i=y.size()-1; i>=0; i--){
            if(mated[i]){
                yFused.add(y.remove(i));
            }
            else {
                y.add(y.get(i).getDaughter());
            }

        }
    }

    private double sqrD(double dx, double dy){ // bc not necess
        return dx*dx + dy * dy;
    }
}
