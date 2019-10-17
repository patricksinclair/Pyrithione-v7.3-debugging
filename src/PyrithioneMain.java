public class PyrithioneMain {
    public static void main(String[] args){

        long startTime = System.currentTimeMillis();

        double[] tau_vals = {0.005, 0.01, 0.02, 0.05, 0.1};
        int nreps = 32;

        for(int i = 0; i < tau_vals.length; i++){
            BioSystem.debugReplications(tau_vals[i], nreps);
            BioSystem.debugDeaths(tau_vals[i], nreps);
            BioSystem.debugImmigration(tau_vals[i], nreps);
            BioSystem.debugDeterioration(tau_vals[i], nreps);
        }

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);

    }
}
