import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.util.*;

public class Yeast2 {

    public static final int LOCAL = 0;
    public static final int GLOBAL = 1;
    public static final int DIPLOID = 2;

    public static final int AXIAL = 0; // PI/3
    public static final int POLAR = 1; // PI

    private static double maxD;

    private double x,y;

    private double oldDivisionSite;

    private int fs;

    private Yeast2 daughter;

    private int yeastType;
    private int divMode;

    public ArrayList<Yeast2> kins;

    private double firstClosestD;

    Yeast2(MersenneTwister rd, double L, int type){  // make a primary yeast

        x = L*rd.nextDouble();
        y = L*rd.nextDouble();
        oldDivisionSite = 2*Math.PI*rd.nextDouble();
        fs = (rd.nextDouble()<0.5)?-1:1;

        setYeastType(type);

        kins = new ArrayList<Yeast2>();
        kins.add(this);

    }
    Yeast2(double xx, double yy, double ods, int futs, int type){ // make a daughter yeast
        x=xx;
        y=yy;
        oldDivisionSite = ods;
        fs = futs;
        setYeastType(type);
    }


    void startDaughter(MersenneTwister rd, ArrayList<Yeast2> y2s) {

        double nds = oldDivisionSite;
        boolean neighbor = false;
        if (yeastType == LOCAL) {
            double d = maxD * maxD;
            for (Yeast2 yt : y2s) {
                double[] r = yt.getPosition();
                double dl = (r[0] - x) * (r[0] - x) + (r[1] - y) * (r[1] - y);

                if (dl <= d) {
                    neighbor = true;
                    d = dl;
                    nds = Math.atan2((r[1] - y), (r[0] - x));
                }

            }
        }
        if (!neighbor) {
            switch (divMode) {
                case AXIAL:
                    nds = axialSearch();
                    break;
                case POLAR:
                    nds = polarSearch();
                    break;
            }
        }

        oldDivisionSite = nds;
        double xn = x+Math.cos(nds);
        double yn = y+Math.sin(nds);
        Yeast2 y2 = new Yeast2(xn,yn,nds-Math.PI,((rd.nextDouble()<0.5)?-1:1),yeastType);
        kins.add(y2);
        y2.kins = this.kins;

        daughter = y2;
    }

    Yeast2 copyYeast(){
        return this;
    }

    double getX(){
        return x;
    }

    double getY(){
        return y;
    }

    int getYeastType(){
        return yeastType;
    }

    Yeast2 getDaughter(){return daughter;}

    void setX(double xx){
        this.x=xx;
    }

    void setY(double yy){
        this.y=yy;
    }

    void setDivMode(int ndm){
        divMode=ndm;
    }

    void setYeastType(int type){
        yeastType = type;
        switch(type){
            case LOCAL:divMode=AXIAL;break;
            case GLOBAL:divMode=POLAR;break;
        }
    }

    void setPosition(double[] r){
        x = r[0];
        y = r[1];
    }

    double[] getPosition(){
        double r[] = new double[2];
        r[0] = x;
        r[1] = y;
        return r;
    }

    private double axialSearch(){
        double[] nextAngle = new double[5];
        nextAngle[0]=Math.PI/3;
        nextAngle[1]=-Math.PI/3;
        nextAngle[2]=2*Math.PI/3;
        nextAngle[3]=-2*Math.PI/3;
        nextAngle[4]=Math.PI;

        return searchNext(nextAngle);
    }

    private double polarSearch(){
        double[] nextAngle = new double[5];
        nextAngle[0]=Math.PI;
        nextAngle[1]=-2*Math.PI/3;
        nextAngle[2]=2*Math.PI/3;
        nextAngle[3]=-Math.PI/3;
        nextAngle[4]=Math.PI/3;

        return searchNext(nextAngle);
    }


    private double searchNext(double[] na){
        double xn,yn;
        double nds = oldDivisionSite;
        boolean found = false;
        int i=0;
        while(!found){
            nds = oldDivisionSite+fs*na[i];
            i++;
            xn = x+Math.cos(nds);
            yn = y+Math.sin(nds);
            found = true;
            for(Yeast2 y2:kins){
                double xx = y2.getX();
                double yy = y2.getY();

                found = found && ((xx-xn)*(xx-xn)+(yy-yn)*(yy-yn)>=0.99);
            }
            found = (i==5) || found;
        }

        return nds;
    }

    public static void setMaxD(double mD){
        maxD = mD;
    }
    public static double getMaxD(){
        return maxD;
    }

    void setFirstClosestD(double a){
        firstClosestD = a;
    }
    double getFirstClosestD(){
        return firstClosestD;
    }
}
