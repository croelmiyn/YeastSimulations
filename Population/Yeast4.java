import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.util.ArrayList;

public class Yeast4 {

    public static final int LOCAL = 0;
    public static final int GLOBAL = 1;

    public static final int AXIAL = 0; // PI/3
    public static final int POLAR = 1; // PI

    private static double DLmax;
    private static double DLmin;

    private double x,y;

    private double oldDivisionSite;

    private int fs;

    private Yeast4 daughter;

    private int yeastTypeEmit;
    private int yeastTypeSense;
    private int divMode;

    private boolean stim;

    private double DL;

    private double p_t;

    public ArrayList<Yeast4> kins;

    private double firstClosestD;

    private double L;

    Yeast4(MersenneTwister rd, double L, int typeSense, int typeEmit){  // make a primary yeast

        x = L*rd.nextDouble();
        y = L*rd.nextDouble();
        oldDivisionSite = 2*Math.PI*rd.nextDouble();
        fs = (rd.nextDouble()<0.5)?-1:1;

        setYeastType(typeSense,typeEmit,false,0);

        kins = new ArrayList<Yeast4>();
        kins.add(this);

        this.L=L;

    }
    Yeast4(double xx, double yy, double ods, int futs, int typeSense, int typeEmit, boolean st, double pt){ // make a daughter yeast
        x=xx;
        y=yy;
        oldDivisionSite = ods;
        fs = futs;
        setYeastType(typeSense,typeEmit,st,pt);
    }


    void startDaughter(MersenneTwister rd, ArrayList<Yeast4> y2s) {

        double nds = oldDivisionSite;
        boolean neighbor = false;
        if (yeastTypeSense == LOCAL && stim) {
            double d = 1.0; // Perfect aim only if distance<DL and stimulation ON
            for (Yeast4 yt : y2s) {
                double[] r = yt.getPosition();
                double dl = (sqrtBc((r[0] - x),(r[1] - y))-1)/yt.getDL();

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
        Yeast4 y2 = new Yeast4(xn,yn,nds-Math.PI,((rd.nextDouble()<0.5)?-1:1),yeastTypeSense,yeastTypeEmit,stim,p_t);
        y2.setL(L);
        kins.add(y2);
        y2.kins = this.kins;

        daughter = y2;
    }

    Yeast4 copyYeast(){
        return this;
    }

    double getX(){
        return x;
    }

    double getY(){
        return y;
    }

    double getDL(){ return DL; }

    int getYeastTypeSense(){
        return yeastTypeSense;
    }
    int getYeastTypeEmit(){
        return yeastTypeEmit;
    }

    Yeast4 getDaughter(){return daughter;}

    void setX(double xx){
        this.x=xx;
    }

    void setY(double yy){
        this.y=yy;
    }

    void setL(double l){ L=l; }

    void setDivMode(){
        divMode= (yeastTypeSense==GLOBAL && stim)?POLAR:AXIAL;
    }

    void setDL(){
        if(yeastTypeEmit==LOCAL) DL=DLmin + p_t * (DLmax-DLmin);
    }

    void setYeastType(int typeSense, int typeEmit, boolean st, double pt){
        yeastTypeSense = typeSense;
        yeastTypeEmit = typeEmit;

        stim = st;
        p_t=pt;
        setDL();
        setDivMode();
    }

    void setStim(MersenneTwister rd, double d,double ratio){
        switch(yeastTypeSense){
            case LOCAL:
                {p_t = (d<1.0)?1.0:0.0;}
                break;
            case GLOBAL:
                p_t = ratio;
                break;
        }
        stim = (rd.nextDouble()<=p_t);
        if(stim){
            if(yeastTypeSense==GLOBAL){
                divMode=POLAR;
            }
        }
        setDL();
    }

    void resetStim(MersenneTwister rd, ArrayList<Yeast4> y2s, int nPeers){
        double p_tp1;
        switch(yeastTypeSense){
            case LOCAL:
                double d = L/DLmin;
                for(Yeast4 yt:y2s){
                    double[] r = yt.getPosition();
                    double dl = (sqrtBc((r[0] - x),(r[1] - y) )-1)/yt.getDL();
                    if (dl <= d) {
                        d = dl;
                    }
                }
                {p_tp1 = (d<1.0)?1.0:0.0;}
                break;
            case GLOBAL:
                p_tp1 = (double)y2s.size()/(y2s.size()+nPeers);
                break;
            default:
                p_tp1=p_t;
        }
        if(!stim && p_tp1>p_t){
            stim = (rd.nextDouble()<=((p_tp1-p_t)/(1-p_t)));
            if(stim){
                if(yeastTypeSense==GLOBAL){
                    divMode=POLAR;
                }
            }
        }
        if(stim && p_tp1<p_t) {
            stim = !(rd.nextDouble() <= (1-p_tp1/p_t));
            if(!stim){
                if(yeastTypeSense==GLOBAL){
                    divMode=AXIAL;
                }
            }
        }
        p_t = p_tp1;

    }

    // Cost is stimulation
    double getCost(){ return p_t; }

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
            for(Yeast4 y2:kins){
                double xx = y2.getX();
                double yy = y2.getY();

                found = found && ((xx-xn)*(xx-xn)+(yy-yn)*(yy-yn)>=0.99);
            }
            found = (i==5) || found;
        }

        return nds;
    }

    public static void setMaxDL(double mD){
        DLmax = mD;
    }
    public static double getMaxDL(){
        return DLmax;
    }
    public static void setMinDL(double mD){
        DLmin = mD;
    }
    public static double getMinDL(){
        return DLmin;
    }

    void setFirstClosestD(double a){
        firstClosestD = a;
    }
    double getFirstClosestD(){
        return firstClosestD;
    }

    private double sqrtBc(double dx, double dy){
        double dxbc = bc(dx);
        double dybc = bc(dy);
        return Math.sqrt(dxbc*dxbc + dybc * dybc);
    }

    private double bc(double x){ // return between -L/2 and L/2
        return x - Math.round(x/L)*L;
    }
}
