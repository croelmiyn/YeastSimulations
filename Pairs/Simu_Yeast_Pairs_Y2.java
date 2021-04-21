import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;
import mpi.rc.IJ.IJutilities.ConcurrencyUtils;

import java.io.FileWriter;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class Simu_Yeast_Pairs_Y2 implements PlugIn {

    private int typ1;
    private int typ2;
    private boolean globalSwitch;

    private double maxDist;
    private double dx;

    private int Nrepeat;

    // outputs
    private double[][] meanSuccess;

    // lock
    private final Object lock = new Object();
    private boolean canceled;
    private boolean parallel;

    // file
    private String saveFile;
    private String dir;

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

        int nX = (int) (6.0/dx);

        meanSuccess = new double[nX][4];

        // parallelization block
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        Future<?>[] futures = new Future[nthreads];

        int p = (int) Math.ceil(nX / (double)nthreads);
        if (p < 1) {
            p = 1;
        }

        final int nR = Nrepeat;
        final int typ1f = typ1;
        final int typ2f = typ2;
        final boolean gsf = globalSwitch;
        final double dxf = dx;
        //IJ.log("toto");
        // running the simulations in parallel
        for (int l = 0; l < nthreads; l++) {

            final int firstI = l * p;

            final int lastI = ( (nX) < (firstI + p) )? nX : (firstI + p);

            futures[l] = ConcurrencyUtils.submit(new Runnable() {
                public void run() {

                    for (int i = firstI; i < lastI; i++) {

                        double[] out = simu(1.0+i*dxf, typ1f, typ2f, gsf, nR);

                        synchronized (lock) {
                            for(int j=0;j<4;j++) {
                                meanSuccess[i][j] = out[j];
                            }
                        }

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

        //IJ.log("" + density + "\t" + ratio + "\t" + meanSuccess + "\t" + stdSuccess);
        save();

    }

    private void computeSerial(){

        int nX = (int) (6.0/dx);

        meanSuccess = new double[nX][4];

        final int nR = Nrepeat;
        final int typ1f = typ1;
        final int typ2f = typ2;
        final boolean gsf = globalSwitch;
        final double dxf = dx;

        for (int l = 0; l < nX; l++) {

            double[] out = simu(1.0+l*dxf, typ1f, typ2f, gsf, nR);

            meanSuccess[l] = out;

        }

        //IJ.log("" + density + "\t" + ratio + "\t" + meanSuccess + "\t" + stdSuccess);
        save();

    }

    private void initiate(){

        canceled = false;

        String yeastType[] = new String[2];
        yeastType[0]="Local";
        yeastType[1]="Global";

        globalSwitch = true;

        maxDist = 2;
        dx=0.10;

        Nrepeat = 1000;

        parallel = false;


        GenericDialog gd = new GenericDialog("parameters");

        gd.addChoice("Yeast_type_1 ",yeastType,"");
        gd.addChoice("Yeast_type_2 ",yeastType,"");
        gd.addCheckbox("Global_2AxialRounds", globalSwitch);

        gd.addNumericField("max_distance ", maxDist, 2);
        gd.addNumericField("resolution ", dx, 2);
        gd.addNumericField("nb_repeats ", Nrepeat, 0);

        gd.addCheckbox("parallel_computing", parallel);

        gd.showDialog();
        if(gd.wasCanceled()){
            canceled=true;
            return;
        }

        typ1 = gd.getNextChoiceIndex();
        typ2 = gd.getNextChoiceIndex();
        globalSwitch = gd.getNextBoolean();

        maxDist = gd.getNextNumber();
        dx = gd.getNextNumber();

        Nrepeat = (int) gd.getNextNumber();

        parallel = gd.getNextBoolean();

        Yeast2.setMaxD(maxDist);
//*
        saveFile = "MatingSuccess_Pairs_"+yeastType[typ1]+"-"+yeastType[typ2]+"_maxD-"+maxDist;
        SaveDialog sd = new SaveDialog("Track_File",saveFile,".txt");
        dir = sd.getDirectory();
        saveFile = sd.getFileName();
//*/
    }

    private double[] simu(double distance, int type1, int type2, boolean gs, int nR){

        SimulationRunPairsY2 sr;

        sr = new SimulationRunPairsY2();

        sr.run(distance,type1,type2,gs,nR);

        return sr.getFusionSuccess();

    }

    private void save(){

        String lineSep = System.getProperty("line.separator");
        try {
            FileWriter file = new FileWriter(dir + saveFile);
            String buffer;

            for(int x=0;x<meanSuccess.length;x++) {
                buffer = "" + (1.0+x*dx);
                for (int i = 0; i < 4; i++) {
                    buffer += "\t"+ meanSuccess[x][i] ;
                }
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

    }
}
