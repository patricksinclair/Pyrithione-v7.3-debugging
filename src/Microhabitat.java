import java.util.ArrayList;

class Microhabitat {

    private double c;
    private ArrayList<Double> population;

    private int K = 120; //karrying kapacity
    private boolean surface = false, biofilm_region, immigration_zone = false;
    private double max_gRate = 0.083; //max growth rate =  2/day
    private double uniform_dRate = 0.018; //all bacteria have this death rate
    double biofilm_threshold; //fraction occupied needed to transition to biofilm
    double b = 0.2; //migration rate

    public Microhabitat(double c, double biofilm_threshold){
        this.c = c;
        this.population = new ArrayList<>(K);
        this.biofilm_threshold = biofilm_threshold;
        this.biofilm_region = false;
    }

    int getN(){return population.size();}
    boolean isSurface(){return surface;}
    boolean isBiofilm_region(){return biofilm_region;}
    boolean isImmigration_zone(){return immigration_zone;}
    ArrayList<Double> getPopulation(){return population;}


    void setSurface(){this.surface = true;}
    void setBiofilm_region(){this.biofilm_region = true;}
    void setImmigration_zone(boolean immigration_zone){this.immigration_zone = immigration_zone;}


    private double fractionFull(){
        return getN()/(double)K;
    }

    boolean atBiofilmThreshold(){
        //double threshold_density = 0.8;
        return fractionFull() >= biofilm_threshold;}

    double migrate_rate(){
        //returns 0.5*b for the microhabitat next to the ship hull, to account for the inability to migrate into the hull
        //also for the microhabitat that's the biofilm edge
        double b = 0.2;
        return (surface || immigration_zone) ? 0.5*b : b;
    }

    private double beta(int index){
        return population.get(index);
    }

    private double phi_c(int index){
        double cB = c/beta(index);
        return 1. - (6.*cB*cB)/(5. + cB*cB);
    }


    double replicationRate(int index){
        double phi_c_scaled = max_gRate*phi_c(index);
        return phi_c_scaled > 0. ? phi_c_scaled*(1. - getN()/(double)K) : 0.;
        //double dRate = phi_c_scaled < 0. ? (phi_c_scaled - uniform_dRate) : -1.*uniform_dRate; //this is -ve to reflect death

        //return gRate;
    }


    double deathRate(int index){
        double phi_c_scaled = max_gRate*phi_c(index);
        return phi_c_scaled < 0. ? (phi_c_scaled - uniform_dRate) : -1.*uniform_dRate; //this is -ve to reflect death
    }


    void addARandomBacterium_x_N(int n_bacteria){
        for(int i = 0; i < n_bacteria; i++){
            //just use a single species for now
            population.add(5.);
        }
    }

    void replicateABacterium_x_N(int index, int nReps){
        for(int i = 0; i < nReps; i++){
            population.add(population.get(index));
        }
    }

    void addABacterium(double MIC){population.add(MIC);}

    void removeABacterium(int index){population.remove(index);}

}
