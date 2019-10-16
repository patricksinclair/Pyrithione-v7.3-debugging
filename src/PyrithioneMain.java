public class PyrithioneMain {
    public static void main(String[] args){

        double[] tau_vals = {0.005, 0.01, 0.02, 0.05, 0.1};

        for(int i = 0; i < tau_vals.length; i++){
            BioSystem.debugReplications(tau_vals[i]);
            BioSystem.debugDeaths(tau_vals[i]);
            BioSystem.debugImmigration(tau_vals[i]);
            BioSystem.debugDeterioration(tau_vals[i]);
        }

    }
}