import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.ArrayList;
import java.util.Random;
import java.util.stream.IntStream;

public class BioSystem {

    private Random rand = new Random();

    private double alpha, c_max; //steepness and max val of antimicrobial concn
    private double scale, sigma; //mic distb shape parameters
    private ArrayList<Microhabitat> microhabitats;
    private double time_elapsed, exit_time; //exit time is the time it took for the biofilm to reach the thickness limit, if it did
    private int immigration_index;

    private double deterioration_rate = 0.002;
    private double biofilm_threshold = 0.6;
    private double immigration_rate = 0.8;
    private double tau;
    private double delta_x = 5.;
    private int thickness_limit = 6; //this is how big the system can get before we exit. should reduce overall simulation duration
    private int n_detachments = 0, n_deaths = 0, n_replications = 0, n_immigrations = 0;
    private int[] migrations_out, migrations_in;

    private BioSystem(int initial_pop, double tau_step){
        this.tau = tau_step;
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.immigration_index = 0;
        this.microhabitats = new ArrayList<>();

        microhabitats.add(new Microhabitat(0., this.biofilm_threshold));
        microhabitats.get(0).setSurface();
        microhabitats.get(0).setImmigration_zone(true);
        microhabitats.get(0).addARandomBacterium_x_N(initial_pop);
    }


    private BioSystem(int initial_pop, double tau_step, int L){
        //this is used to create a system of more than one microhabitat
        //currently it's just used for migraiton testing, hence the high initial pops

        this.tau = tau_step;
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.immigration_index = 0;
        this.migrations_out = new int[L];
        this.migrations_in = new int[L];

        for(int i = 0; i < L; i++){
            microhabitats.add(new Microhabitat(0., this.biofilm_threshold));
            microhabitats.get(i).addARandomBacterium_x_N(initial_pop);
        }
        microhabitats.get(0).setSurface();
        microhabitats.get(L-1).setImmigration_zone(true);
    }



    private int getN_detachments(){
        return n_detachments;
    }

    private int getN_deaths(){
        return n_deaths;
    }

    private int getN_replications(){
        return n_replications;
    }

    private int getN_immigrations(){
        return n_immigrations;
    }

    private double getTimeElapsed(){
        return time_elapsed;
    }

    private double getExit_time(){
        return exit_time;
    }

    private int getSystemSize(){
        return microhabitats.size();
    }

    private int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }


    private void immigrate(int mh_index, int n_immigrants){
        microhabitats.get(mh_index).addARandomBacterium_x_N(n_immigrants);
    }


    public void migrate(int mh_index, int bac_index){

        double migrating_bac = microhabitats.get(mh_index).getPopulation().get(bac_index);
        microhabitats.get(mh_index).removeABacterium(bac_index);

        if(microhabitats.get(mh_index).isSurface()) {
            microhabitats.get(mh_index + 1).addABacterium(migrating_bac);
            migrations_in[mh_index+1]++;
        } else if(microhabitats.get(mh_index).isImmigration_zone()) {
            microhabitats.get(mh_index - 1).addABacterium(migrating_bac);
            migrations_in[mh_index-1]++;
        } else {
            if(rand.nextBoolean()) {
                microhabitats.get(mh_index + 1).addABacterium(migrating_bac);
                migrations_in[mh_index+1]++;
            } else {
                microhabitats.get(mh_index - 1).addABacterium(migrating_bac);
                migrations_in[mh_index-1]++;
            }
        }
    }




    public void performAction_replication(){
        //this is a heavily reduced version of the usual method, designed  for debugging
        //this method only allows for replication to occur.
        double tau_step = tau;
        int system_size = microhabitats.size();
        int[][] replication_allocations;
        int[] original_popsizes;

        whileloop:
        while(true){
            replication_allocations = new int[system_size][];
            original_popsizes = new int[system_size];

            //calculate events
            for(int mh_index = 0; mh_index < system_size; mh_index++){

                int mh_pop = microhabitats.get(mh_index).getN();

                int[] n_replications = new int[mh_pop];
                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    double g_rate = microhabitats.get(mh_index).replicationRate(bac_index);
                    if(g_rate == 0.) {
                        n_replications[bac_index] = 0;
                    } else {
                        n_replications[bac_index] = new PoissonDistribution(g_rate*tau_step).sample();
                    }
                }

                replication_allocations[mh_index] = n_replications;
                original_popsizes[mh_index] = mh_pop;
            }
            break;
        }

        //carry out the events
        for(int mh_index = 0; mh_index < system_size; mh_index++){
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                microhabitats.get(mh_index).replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);
                n_replications += replication_allocations[mh_index][bac_index];
            }
        }

        time_elapsed += tau_step;
    }



    public void performAction_death(){
        //only allows for death
        double tau_step = tau;
        int system_size = microhabitats.size();

        int[][] death_allocations;
        int[] original_popsizes;

        whileloop:
        while(true){
            death_allocations = new int[system_size][];
            original_popsizes = new int[system_size];

            for(int mh_index = 0; mh_index < system_size; mh_index++){
                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_deaths = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    double d_rate = Math.abs(microhabitats.get(mh_index).deathRate(bac_index));
                    if(d_rate == 0.) {
                        n_deaths[bac_index] = 0;
                    } else {
                        n_deaths[bac_index] = new PoissonDistribution(d_rate*tau_step).sample();

                        if(n_deaths[bac_index] > 1) {
                            tau_step /= 2.;
                            continue whileloop;
                        }
                    }
                }
                death_allocations[mh_index] = n_deaths;
                original_popsizes[mh_index] = mh_pop;
            }
            break;
        }


        for(int mh_index = 0; mh_index < system_size; mh_index++) {
            for(int bac_index = original_popsizes[mh_index] - 1; bac_index >= 0; bac_index--) {

                if(death_allocations[mh_index][bac_index] != 0) {
                    microhabitats.get(mh_index).removeABacterium(bac_index);
                    n_deaths++;
                }
            }
        }

        time_elapsed += tau_step;
    }


    public void performAction_deterioration(){
        //only allows deterioration
        double tau_step = tau;
        int system_size = microhabitats.size();
        int[] detachment_allocations;
        int[] original_popsizes;

        whileloop:
        while(true){
            detachment_allocations = new int[microhabitats.get(immigration_index).getN()];
            original_popsizes = new int[system_size];

            for(int mh_index = 0; mh_index < system_size; mh_index++){

                int mh_pop = microhabitats.get(mh_index).getN();

                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    if(mh_index == immigration_index) {
                        detachment_allocations[bac_index] = new PoissonDistribution(deterioration_rate*tau_step).sample();

                        if(detachment_allocations[bac_index] > 1) {
                            tau_step /= 2.;
                            continue whileloop;
                        }
                    }
                }
                original_popsizes[mh_index] = mh_pop;
            }
            break;
        }


        for(int mh_index = 0; mh_index < system_size; mh_index++){
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(mh_index == immigration_index) {
                    if(detachment_allocations[bac_index] != 0) {
                        microhabitats.get(mh_index).removeABacterium(bac_index);
                        n_detachments++;
                    }
                }
            }
        }

        time_elapsed+=tau_step;
    }




    public void performAction_immigration(){
        //only allows for immigration
        double tau_step = tau;
        int n_immigrants;


        n_immigrants = new PoissonDistribution(immigration_rate*tau_step).sample();
        immigrate(immigration_index, n_immigrants);
        n_immigrations += n_immigrants;
        time_elapsed += tau_step;
    }



    public void performAction_migration(){
        //only allows migration between microhabitats
        double tau_step = tau;
        int system_size = microhabitats.size();
        int[][] migration_allocations;
        int[] original_popsizes;

        whileloop:
        while(true){
            migration_allocations = new int[system_size][];
            original_popsizes = new int[system_size];

            for(int mh_index = 0; mh_index < system_size; mh_index++){

                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_migrations = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++){

                    n_migrations[bac_index] = new PoissonDistribution(microhabitats.get(mh_index).migrate_rate()*tau_step).sample();

                    if(n_migrations[bac_index] > 1) {
                        tau_step /= 2.;
                        continue whileloop;
                    }
                }

                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = mh_pop;
            }
            break;
        }


        for(int mh_index = 0; mh_index < system_size; mh_index++){
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){
                if(system_size > 1) {
                    if(migration_allocations[mh_index][bac_index] != 0) migrate(mh_index, bac_index);
                    migrations_out[mh_index]++;
                }
            }
        }
        time_elapsed += tau_step;
    }



    public static void debugReplications(double tau_val){
        long startTime = System.currentTimeMillis();
        //method used to debug the replication events
        //run 1 microhab for 100 hours, do 16 reps and average them
        int initial_pop = 5;
        double duration = 100.;
        int nreps = 16;
        int nmeasurements = 50;
        double interval = duration/nmeasurements;

        String directoryName = "debugging";
        String filename = String.format("replication_debugging_t=%.2f_initial_pop=%d_tau=%.3f", duration, initial_pop, tau_val);

        Databox[] databoxes = new Databox[nreps];
        //todo make sure correct method here
        IntStream.range(0, nreps).parallel().forEach(i ->
                databoxes[i] = debugReplications_subroutine(i, nmeasurements, duration, initial_pop, tau_val));


        Databox avg_databox = Databox.averagedResults(databoxes);
        Toolbox.writeDataboxToFile(directoryName, filename, avg_databox);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }

    public static Databox debugReplications_subroutine(int i, int nMeasurements, double duration, int initial_pop, double tau_step){

        double interval = duration/(double)nMeasurements;
        boolean alreadyRecorded = false;
        int sampleCounter = 0;
        double[] times = new double[nMeasurements+1];
        double[] event_counts = new double[nMeasurements+1];
        double[] pop_sizes = new double[nMeasurements+1];

        BioSystem bs = new BioSystem(initial_pop, tau_step);

        while(bs.time_elapsed <= duration+0.2*interval){

            if(bs.time_elapsed%interval < tau_step && !alreadyRecorded){
                System.out.println("replication -- "+"\ttau: "+tau_step+"\trep: "+i+"\tt: "+bs.time_elapsed+"\tpop_size: "+bs.getTotalN());
                times[sampleCounter] = bs.time_elapsed;
                //todo - get correct counter here
                event_counts[sampleCounter] = bs.n_replications;
                pop_sizes[sampleCounter] = bs.getTotalN();

                alreadyRecorded = true;
                sampleCounter++;
            }
            if(bs.time_elapsed%interval > 0.1) alreadyRecorded = false;
            //todo - make sure the correct method used here
            bs.performAction_replication();
        }
        return new Databox(tau_step, times, event_counts, pop_sizes);
    }



    public static void debugDeaths(double tau_val){
        long startTime = System.currentTimeMillis();
        //method used to debug the replication events
        //run 1 microhab for 100 hours, do 10 reps and average them
        //todo make sure correct pop used
        int initial_pop = 120;
        double duration = 100.;
        int nreps = 16;
        int nmeasurements = 50;
        double interval = duration/nmeasurements;

        String directoryName = "debugging";
        //todo - change name here
        String filename = String.format("death_debugging_t=%.2f_initial_pop=%d_tau=%.3f", duration, initial_pop, tau_val);

        Databox[] databoxes = new Databox[nreps];
        //todo - make sure correct method here
        IntStream.range(0, nreps).parallel().forEach(i ->
                databoxes[i] = debugDeaths_subroutine(i, nmeasurements, duration, initial_pop, tau_val));


        Databox avg_databox = Databox.averagedResults(databoxes);
        Toolbox.writeDataboxToFile(directoryName, filename, avg_databox);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }

    public static Databox debugDeaths_subroutine(int i, int nMeasurements, double duration, int initial_pop, double tau_step){

        double interval = duration/(double)nMeasurements;
        boolean alreadyRecorded = false;
        int sampleCounter = 0;
        double[] times = new double[nMeasurements+1];
        double[] event_counts = new double[nMeasurements+1];
        double[] pop_sizes = new double[nMeasurements+1];

        BioSystem bs = new BioSystem(initial_pop, tau_step);

        while(bs.time_elapsed <= duration+0.2*interval){

            if(bs.time_elapsed%interval < tau_step && !alreadyRecorded){
                System.out.println("death -- "+"\ttau: "+tau_step+"\trep: "+i+"\tt: "+bs.time_elapsed+"\tpop_size: "+bs.getTotalN());
                times[sampleCounter] = bs.time_elapsed;
                //todo - get correct counter here
                event_counts[sampleCounter] = bs.n_deaths;
                pop_sizes[sampleCounter] = bs.getTotalN();

                alreadyRecorded = true;
                sampleCounter++;
            }
            if(bs.time_elapsed%interval > 0.1) alreadyRecorded = false;
            //todo - make sure the correct method used here
            bs.performAction_death();
        }
        return new Databox(tau_step, times, event_counts, pop_sizes);
    }


    public static void debugImmigration(double tau_val){
        long startTime = System.currentTimeMillis();
        //method used to debug the replication events
        //run 1 microhab for 100 hours, do 10 reps and average them
        //todo make sure correct pop used
        int initial_pop = 5;
        double duration = 100.;
        int nreps = 16;
        int nmeasurements = 50;
        double interval = duration/nmeasurements;

        String directoryName = "debugging";
        //todo - change name here
        String filename = String.format("immigration_debugging_t=%.2f_initial_pop=%d_tau=%.3f", duration, initial_pop, tau_val);

        Databox[] databoxes = new Databox[nreps];
        //todo - make sure correct method here
        IntStream.range(0, nreps).parallel().forEach(i ->
                databoxes[i] = debugImmigration_subroutine(i, nmeasurements, duration, initial_pop, tau_val));


        Databox avg_databox = Databox.averagedResults(databoxes);
        Toolbox.writeDataboxToFile(directoryName, filename, avg_databox);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }

    public static Databox debugImmigration_subroutine(int i, int nMeasurements, double duration, int initial_pop, double tau_step){

        double interval = duration/(double)nMeasurements;
        boolean alreadyRecorded = false;
        int sampleCounter = 0;
        double[] times = new double[nMeasurements+1];
        double[] event_counts = new double[nMeasurements+1];
        double[] pop_sizes = new double[nMeasurements+1];

        BioSystem bs = new BioSystem(initial_pop, tau_step);

        while(bs.time_elapsed <= duration+0.2*interval){

            if(bs.time_elapsed%interval < tau_step && !alreadyRecorded){
                System.out.println("immigration -- "+"\ttau: "+tau_step+"\trep: "+i+"\tt: "+bs.time_elapsed+"\tpop_size: "+bs.getTotalN());
                times[sampleCounter] = bs.time_elapsed;
                //todo - get correct counter here
                event_counts[sampleCounter] = bs.n_immigrations;
                pop_sizes[sampleCounter] = bs.getTotalN();

                alreadyRecorded = true;
                sampleCounter++;
            }
            if(bs.time_elapsed%interval > 0.1) alreadyRecorded = false;
            //todo - make sure the correct method used here
            bs.performAction_immigration();
        }
        return new Databox(tau_step, times, event_counts, pop_sizes);
    }


    public static void debugDeterioration(double tau_val){
        long startTime = System.currentTimeMillis();
        //method used to debug the replication events
        //run 1 microhab for 100 hours, do 10 reps and average them
        //todo make sure correct pop used
        int initial_pop = 120;
        double duration = 100.;
        int nreps = 16;
        int nmeasurements = 50;
        double interval = duration/nmeasurements;

        String directoryName = "debugging";
        //todo - change name here
        String filename = String.format("deterioration_debugging_t=%.2f_initial_pop=%d_tau=%.3f", duration, initial_pop, tau_val);

        Databox[] databoxes = new Databox[nreps];
        //todo - make sure correct method here
        IntStream.range(0, nreps).parallel().forEach(i ->
                databoxes[i] = debugDeterioration_subroutine(i, nmeasurements, duration, initial_pop, tau_val));


        Databox avg_databox = Databox.averagedResults(databoxes);
        Toolbox.writeDataboxToFile(directoryName, filename, avg_databox);

        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }

    public static Databox debugDeterioration_subroutine(int i, int nMeasurements, double duration, int initial_pop, double tau_step){

        double interval = duration/(double)nMeasurements;
        boolean alreadyRecorded = false;
        int sampleCounter = 0;
        double[] times = new double[nMeasurements+1];
        double[] event_counts = new double[nMeasurements+1];
        double[] pop_sizes = new double[nMeasurements+1];

        BioSystem bs = new BioSystem(initial_pop, tau_step);

        while(bs.time_elapsed <= duration+0.2*interval){

            if(bs.time_elapsed%interval < tau_step && !alreadyRecorded){
                System.out.println("deterioration -- "+"\ttau: "+tau_step+"\trep: "+i+"\tt: "+bs.time_elapsed+"\tpop_size: "+bs.getTotalN());
                times[sampleCounter] = bs.time_elapsed;
                //todo - get correct counter here
                event_counts[sampleCounter] = bs.n_detachments;
                pop_sizes[sampleCounter] = bs.getTotalN();

                alreadyRecorded = true;
                sampleCounter++;
            }
            if(bs.time_elapsed%interval > 0.1) alreadyRecorded = false;
            //todo - make sure the correct method used here
            bs.performAction_deterioration();
        }
        return new Databox(tau_step, times, event_counts, pop_sizes);
    }



















}
