import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

public class Toolbox {

    private static int largestHeaderLength(String[] headers){
        int biggun = 0;
        for(String s : headers){
            if(s.length() > biggun) biggun = s.length();
        }
        return biggun;
    }

    public static void writeDataboxToFile(String directoryName, String filename, Databox db){

        File directory = new File(directoryName);
        if(!directory.exists()) directory.mkdirs();

        File file = new File(directoryName+"/"+filename+".txt");

        try{

            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
            String[] headers = {"#times", "event_counts", "pop_size"};
            int string_length = Math.max(12, Toolbox.largestHeaderLength(headers)+3);

            String head_start = "#"+headers[0]+",";
            String file_header = String.format("%-"+string_length+"s", head_start);
            for(int i = 1; i < headers.length-1; i++){
                String heado = headers[i]+",";
                file_header += String.format("%-"+string_length+"s", heado);
            }
            String heado = headers[headers.length-1];
            file_header += String.format("%-"+string_length+"s", heado);
            bw.write(file_header);
            bw.newLine();


            for(int t = 0; t < db.getTimes().length; t++){

                String output = "";
                String t_val = String.format("%.3E,", db.getTimes()[t]);
                output += String.format("%-"+string_length+"s", t_val);
                String event_val = String.format("%.3E,", db.getEvent_counters()[t]);
                output += String.format("%-"+string_length+"s", event_val);
                String pop_val =  String.format("%.3E", db.getPop_sizes()[t]);
                output += String.format("%-"+string_length+"s", pop_val);


                bw.write(output);
                bw.newLine();
            }

            bw.close();

        }catch (IOException e){}

    }























    public static String millisToShortDHMS(long duration) {
        String res = "";
        long days  = TimeUnit.MILLISECONDS.toDays(duration);
        long hours = TimeUnit.MILLISECONDS.toHours(duration)
                - TimeUnit.DAYS.toHours(TimeUnit.MILLISECONDS.toDays(duration));
        long minutes = TimeUnit.MILLISECONDS.toMinutes(duration)
                - TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(duration));
        long seconds = TimeUnit.MILLISECONDS.toSeconds(duration)
                - TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(duration));
        if (days == 0) {
            res = String.format("%02d:%02d:%02d", hours, minutes, seconds);
        }
        else {
            res = String.format("%dd%02d:%02d:%02d", days, hours, minutes, seconds);
        }
        return res;
    }
}
