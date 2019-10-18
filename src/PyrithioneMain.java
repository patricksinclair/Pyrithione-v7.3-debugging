import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.AbstractRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Random;

public class PyrithioneMain {
    public static void main(String[] args){

        //        double tau = 0.01;
        //        double attempts = 100./tau;
        //        Random rand = new Random();
        //        //double x = new PoissonDistribution(rand, 0.002*tau).probability(1)*attempts;
        //        //System.out.println(x);
        //        int counter = 0;
        //
        //        for(int b = 0; b < 120; b++){
        //            for(int i = 0; i < (int)attempts; i++){
        //                PoissonDistribution poss = new PoissonDistribution(0.002*tau);
        //                poss.reseedRandomGenerator(rand.nextLong());
        //                int j = poss.sample();
        //                counter+=j;
        //                //System.out.println(j);
        //            }
        //
        //        }
        //        System.out.println("c: "+counter);
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
