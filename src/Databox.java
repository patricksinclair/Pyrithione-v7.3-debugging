public class Databox {

    private double tau;
    private double[] times;
    private double[] event_counters;
    private double[] pop_sizes;
    private double[][] event_counters_matrix;
    private double[][] pop_sizes_matrix;



    public Databox(double tau, double[] times, double[] event_counters, double[] pop_sizes){
        this.tau = tau;
        this.times = times;
        this.event_counters = event_counters;
        this.pop_sizes = pop_sizes;
    }

    public Databox(double tau, double[] times, double[][] event_counters_matrix, double[][] pop_sizes_matrix){
        this.tau = tau;
        this.times = times;
        this.event_counters_matrix = event_counters_matrix;
        this.pop_sizes_matrix = pop_sizes_matrix;
    }



    public double getTau(){
        return tau;
    }

    public double[] getTimes(){
        return times;
    }

    public double[] getEvent_counters(){
        return event_counters;
    }

    public double[] getPop_sizes(){
        return pop_sizes;
    }

    public double[][] getEvent_counters_matrix(){
        return event_counters_matrix;
    }

    public double[][] getPop_sizes_matrix(){
        return pop_sizes_matrix;
    }


    public static Databox averagedResults(Databox[] databoxes){

        double[] avg_event_counters = new double[databoxes[0].getEvent_counters().length];
        double[] avg_pop_sizes = new double[databoxes[0].getEvent_counters().length];

        for(int i = 0; i < databoxes[0].getEvent_counters().length; i++){
            for(int db = 0; db < databoxes.length; db++){

                avg_event_counters[i] += databoxes[db].getEvent_counters()[i]/(double)databoxes.length;
                avg_pop_sizes[i] += databoxes[db].getPop_sizes()[i]/(double)databoxes.length;
            }

        }

        return new Databox(databoxes[0].tau, databoxes[0].times, avg_event_counters, avg_pop_sizes);
    }



    public static Databox averagedResultsMatrix(Databox[] databoxes){

        double[][] avg_event_counters_matrix = new double[databoxes[0].getEvent_counters_matrix().length][databoxes[0].getEvent_counters_matrix()[0].length];
        double[][] avg_pop_sizes_matrix = new double[databoxes[0].getPop_sizes_matrix().length][databoxes[0].getPop_sizes_matrix()[0].length];

        for(int i = 0; i < databoxes[0].getEvent_counters_matrix().length; i++){
            for(int j = 0; j < databoxes[0].getEvent_counters_matrix()[0].length; j++){

                for(int db = 0; db < databoxes.length; db++){
                    avg_event_counters_matrix[i][j] += databoxes[db].getEvent_counters_matrix()[i][j]/(double)databoxes.length;
                    avg_pop_sizes_matrix[i][j] += databoxes[db].getPop_sizes_matrix()[i][j]/(double)databoxes.length;
                }
            }
        }

        return new Databox(databoxes[0].tau, databoxes[0].times, avg_event_counters_matrix, avg_pop_sizes_matrix);
    }







}
