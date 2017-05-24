import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

/*
 * This program,
 * completed the data of the equal width of the discrete,
 * completed the conversion of csv data, each line in the data is a property, each column is a record,
 * the data is transposed after completion.
 */
public class preTreatment {

 @SuppressWarnings("unchecked")
public static void main(String[] args) {
	
  int l = Integer.valueOf(args[0]);		//columns number of input file
  int n = Integer.valueOf(args[1]);		//dispersion
  String f = args[2];					//the path of input file
  String fout = args[3];				//the path of output file
  try {
	  long startTime = System.currentTimeMillis(); 
	  File file = new File(fout);

	  //If file doesn't exists,create it.
	  if (!file.exists()) {
		  file.createNewFile();
	  }

	  try {
		  //the path of input file
		  BufferedReader reader = new BufferedReader(new FileReader(f));
          String line = null;
          @SuppressWarnings("rawtypes")
          Vector Data[] = new Vector[l - 1];
          reader.readLine();
          System.out.println("Start.");
          FileWriter fw = new FileWriter(file.getAbsoluteFile());
    	  BufferedWriter bw = new BufferedWriter(fw);
          while((line=reader.readLine())!=null){
        	    //remove quotation marks
        	  	line = line.replace("\"", "");
				//Divide with a comma
              	String item[] = line.split(",");           	
				//'i' represents the number of columns, each column is a record;
				//'count' indicates the number of genes involved in the association,
				//Purpose:Looking for the association between genes                         	
				for(int i = 0;i < item.length-1;i++){
					double d = Double.parseDouble(item[i+1]);
					//The performance of all transcriptions ranges between [0,1], separated by letters '*'
					//According to the level of performance,we construct a matrix that can be used for analysis.
					if(Data[i] == null){
              			Data[i] = new Vector<String>();
              		}
					int r = (int)(d*n);
					for(int k = 0;k <= n;k++){
						if(r == k && r !=n){
							Data[i].addAll(Arrays.asList(item[0]+"*" + k));
							break;
						}
						if(r == n){
							Data[i].addAll(Arrays.asList(item[0]+"*" + (r-1)));
							break;
						}
					}
				}
          	}
          for(int i = 0;i < l - 1;i++){
      		  //System.out.println(Data[i].size());
      		  bw.write(Data[i].toString().replace("[", "").replace("]", "").replace(" ", "") + "\r\n");      		  
      	  }
    	  bw.close();
          } catch (Exception e) {
          e.printStackTrace();
      }
	  long endTime = System.currentTimeMillis();
	  System.out.println("Use " + (endTime - startTime)/1000.0 + "s.");
	  System.out.println("Done.");

  	} catch (IOException e) {
  		e.printStackTrace();
  	}
 }
}