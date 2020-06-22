
import java.io.*;
import java.util.*;

class format1
{

  public static void main(String args[])
  {
     try
     {
         BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/selected."+args[0]+".-21.features")));
         PrintWriter p = new PrintWriter("../SVM/"+args[0]+".formatted");

         String data = "";

         while((data = r1.readLine()) != null)
         {
            StringTokenizer st = new StringTokenizer(data);
            String s = "1 ";
 
            for(int i=1; i<=20; i++)
             s += i+":"+st.nextToken()+" ";
             
            p.println(s);
         }   
            
        r1.close();
        p.close();     
         
     
     }

     catch(Exception e)
     {
         e.printStackTrace();
     }
  }

}

                      
   


