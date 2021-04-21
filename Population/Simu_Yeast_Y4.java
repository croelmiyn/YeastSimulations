import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import mpi.rc.IJ.IJutilities.ConcurrencyUtils;

import java.io.FileWriter;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import static java.lang.Thread.sleep;

public class Simu_Yeast_Y4 implements PlugIn {

    private int typ1;
    private int typ2;
    private int L;

    private double density;
    private double ratio;

    private double maxDist;
    private double minDist;

    private int Nrepeat;

    // outputs
    private double[][] meanCost;
    private double[][] stdCost;

    private double[] meanN1;
    private double[] stdN1;

    private double[] meanN2;
    private double[] stdN2;

    private double[] meanNF;
    private double[] stdNF;

    ArrayList<double[]>[] distances;

    // lock
    private final Object lock = new Object();
    private boolean canceled;
    private boolean parallel;

    // file
    private String saveFile;
    private String dir;
    private boolean write;
    private boolean diploDivide;

    public void run(String arg){

        initiate();

        if(!canceled) {

            if(parallel) {
                compute();
            }
            else {
                computeSerial();
            }
        }

    }

    private void compute(){

        meanCost = new double[2][4];
        stdCost  = new double[2][4];

        meanN1 = new double[4];
        stdN1  = new double[4];

        meanN2 = new double[4];
        stdN2  = new double[4];

        meanNF = new double[4];
        stdNF  = new double[4];

        distances = new ArrayList[4];
        for(int j=0;j<4;j++){
            distances[j]=new ArrayList<double[]>();
        }

        // parallelization block
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        Future<?>[] futures = new Future[nthreads];

        int p = (int) Math.ceil(Nrepeat / (double)nthreads);
        if (p < 1) {
            p = 1;
        }

        final double d = density;
        final double r = ratio;
        final int Lf = L;
        final int typ1f = typ1;
        final int typ2f = typ2;
        //IJ.log("toto");
        // running the simulations in parallel
        for (int l = 0; l < nthreads; l++) {

            final int firstI = l * p;

            final int lastI = ( (Nrepeat) < (firstI + p) )? Nrepeat : (firstI + p);

            futures[l] = ConcurrencyUtils.submit(new Runnable() {
                public void run() {
                    try {
                        sleep(firstI); // to make sure they all start at a different millisecond
                    } catch (InterruptedException e) {

                    }
                    for (int i = firstI; i < lastI; i++) {

                        simu(d, r, Lf, typ1f, typ2f, i);

                    }

                }
            });
        }
        try {
            ConcurrencyUtils.waitForCompletion(futures);
        } catch (InterruptedException ex) {
            IJ.log("Interruption Exception");
            IJ.log("Message --> " + ex.getMessage());
            IJ.log("Cause --> " + ex.getCause());
            IJ.log("LocMessage --> " + ex.getLocalizedMessage());
        } catch (ExecutionException ex) {
            IJ.log("Execution Exception");
            IJ.log("Message --> " + ex.getMessage());
            IJ.log("Cause --> " + ex.getCause());
            IJ.log("LocMessage --> " + ex.getLocalizedMessage());
        }

        for(int j=0; j<4; j++) {

            for(int k=0;k<2;k++) {
                meanCost[k][j] /= Nrepeat;
                stdCost[k][j] /= Nrepeat;
                stdCost[k][j] -= meanCost[k][j] * meanCost[k][j];
            }

            meanN1[j] /= Nrepeat;
            stdN1[j] /= Nrepeat;
            stdN1[j] -= meanN1[j] * meanN1[j];

            meanN2[j] /= Nrepeat;
            stdN2[j] /= Nrepeat;
            stdN2[j] -= meanN2[j] * meanN2[j];

            meanNF[j] /= Nrepeat;
            stdNF[j] /= Nrepeat;
            stdNF[j] -= meanNF[j] * meanNF[j];
        }

        //IJ.log("" + density + "\t" + ratio + "\t" + meanSuccess + "\t" + stdSuccess);
        save();

    }

    private void computeSerial(){

        meanCost = new double[2][4];
        stdCost  = new double[2][4];

        meanN1 = new double[4];
        stdN1  = new double[4];

        meanN2 = new double[4];
        stdN2  = new double[4];

        meanNF = new double[4];
        stdNF  = new double[4];

        distances = new ArrayList[4];
        for(int j=0;j<4;j++){
            distances[j]=new ArrayList<double[]>();
        }

        final double d = density;
        final double r = ratio;
        final int Lf = L;
        final int typ1f = typ1;
        final int typ2f = typ2;
        IJ.log("toto");

        for (int l = 0; l < Nrepeat; l++) {

            simu(d, r, Lf, typ1f, typ2f, l);

        }

        for(int j=0; j<4; j++) {

            for(int k=0;k<2;k++) {
                meanCost[k][j] /= Nrepeat;
                stdCost[k][j] /= Nrepeat;
                stdCost[k][j] -= meanCost[k][j] * meanCost[k][j];
            }

            meanN1[j] /= Nrepeat;
            stdN1[j] /= Nrepeat;
            stdN1[j] -= meanN1[j] * meanN1[j];

            meanN2[j] /= Nrepeat;
            stdN2[j] /= Nrepeat;
            stdN2[j] -= meanN2[j] * meanN2[j];

            meanNF[j] /= Nrepeat;
            stdNF[j] /= Nrepeat;
            stdNF[j] -= meanNF[j] * meanNF[j];
        }

        //IJ.log("" + density + "\t" + ratio + "\t" + meanSuccess + "\t" + stdSuccess);
        save();

    }

    private void initiate(){

        canceled = false;

        String yeastType[] = new String[2];
        yeastType[Yeast4.LOCAL]="Local";
        yeastType[Yeast4.GLOBAL]="Global";

        L=128;

        density = 0.05;
        ratio = 0.5;

        maxDist = 3;
        minDist = 2;

        Nrepeat = 10;

        parallel = false;
        diploDivide = false;

        GenericDialog gd = new GenericDialog("parameters");

        gd.addChoice("Yeast_type_1 ",yeastType,"");
        gd.addChoice("Yeast_type_2 ",yeastType,"");

        gd.addNumericField("Size (px)", L, 0);
        gd.addNumericField("density ", density, 2);
        gd.addNumericField("ratio ", ratio, 2);
        gd.addNumericField("min_distance ", minDist, 2);
        gd.addNumericField("max_distance ", maxDist, 2);
        gd.addNumericField("nb_repeats ", Nrepeat, 0);

        gd.addCheckbox("print", false);

        gd.addCheckbox("parallel_computing", parallel);
        gd.addCheckbox("diploid_division", diploDivide);

        gd.showDialog();
        if(gd.wasCanceled()){
            canceled=true;
            return;
        }

        typ1 = gd.getNextChoiceIndex();
        typ2 = gd.getNextChoiceIndex();

        L = (int)  gd.getNextNumber();

        density = gd.getNextNumber();
        ratio =   gd.getNextNumber();
        minDist = gd.getNextNumber();
        maxDist = gd.getNextNumber();
        if(maxDist<minDist){IJ.log("wrong distance input"); canceled=true;return;}

        Nrepeat = (int) gd.getNextNumber();

        write = gd.getNextBoolean();
        parallel = gd.getNextBoolean();
        diploDivide = gd.getNextBoolean();

        Yeast4.setMinDL(minDist-1);
        Yeast4.setMaxDL(maxDist-1);

        String yt[] = new String[2];
        yt[Yeast4.LOCAL]="L";
        yt[Yeast4.GLOBAL]="G";
//*
        saveFile = "MatingSuccess_"+yt[typ1]+yt[typ2]+"_d-"+density+"_r-"+ratio+"_L-"+L+"_minD-"+minDist+"_maxD-"+maxDist;
        SaveDialog sd = new SaveDialog("Track_File",saveFile,".txt");
        dir = sd.getDirectory();
        saveFile = sd.getFileName();
//*/
    }

    private void simu(double density, double ratio, int size, int type1, int type2, int n){

        SimulationRunY4 sr;
        boolean wt, dD;

        synchronized(lock){

            sr = new SimulationRunY4(type1,type2);
            wt = write;
            dD = diploDivide;
        }

        sr.run(density,ratio,size,wt,dD);

        if(wt){
            sr.saveImage(dir,n);
        }

        double[][] cost = sr.getCost();
        double[] n1 = sr.getN1();
        double[] n2 = sr.getN2();
        double[] nf = sr.getNF();
        ArrayList<double[]>[] dist = sr.getMatingDistances();

        synchronized (lock) {

            for(int j=0; j<4; j++) {

                for(int k=0; k<2; k++) {
                    meanCost[k][j] +=cost[k][j];
                    stdCost[k][j] += cost[k][j] * cost[k][j];
                }

                meanN1[j] += n1[j];
                stdN1[j] += n1[j] * n1[j];

                meanN2[j] += n2[j];
                stdN2[j] += n2[j] * n2[j];

                meanNF[j] += nf[j];
                stdNF[j] += nf[j] * nf[j];

                distances[j].addAll(dist[j]);
            }
        }

    }

    private void save(){

        String lineSep = System.getProperty("line.separator");
        try {
            FileWriter file = new FileWriter(dir + saveFile);

            for (int i = 0; i < 4; i++){
                String buffer = "" + density + "\t" + ratio ;
                buffer += "\t" + meanN1[i] + "\t" + Math.sqrt(stdN1[i]) + "\t" + meanN2[i] + "\t" + Math.sqrt(stdN2[i]) ;
                buffer += "\t" + meanCost[0][i] + "\t" + Math.sqrt(stdCost[0][i]) + "\t" + meanCost[1][i] + "\t" + Math.sqrt(stdCost[1][i]) ;
                buffer += "\t" + meanNF[i] + "\t" + Math.sqrt(stdNF[i]);
                buffer += lineSep;
                file.write(buffer);
            }

            file.close();
        } catch (Exception e){
            IJ.log("Error saveYeast --> "+e.getMessage());
            IJ.log("Error saveYeast --> "+e.getCause());
            IJ.log("Error saveYeast --> "+e.getLocalizedMessage());
        }
        IJ.showStatus("Done");


        try {
            FileWriter file = new FileWriter(dir + "Distances_" + saveFile);

            String buffer = ""+ distances[0].size();
            for (int i = 1; i < 4; i++) {
                buffer +=  "\t" + distances[i].size();
            }
            buffer+=lineSep;
            file.write(buffer);
            for (int i = 0; i < 4; i++){
                for(double[] d:distances[i]){
                    buffer = "" + d[0] + "\t" + d[1] + "\t" + d[2] + lineSep;
                    file.write(buffer);
                }

            }

            file.close();
        } catch (Exception e){
            IJ.log("Error saveYeast --> "+e.getMessage());
            IJ.log("Error saveYeast --> "+e.getCause());
            IJ.log("Error saveYeast --> "+e.getLocalizedMessage());
        }
        IJ.showStatus("Done");

    }
}
