//calculates the MFEI2 and MFEI3.

import java.io.*;
import java.util.*;

class mfe23
{

  public static void main(String args[])
  {
     try
     {
         BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".data2")));
         
         PrintWriter p = new PrintWriter("../data/"+args[0]+".mfe23");

         String data = null;
                   
         r1.readLine(); //to remove the header        
 
         while((data = r1.readLine()) != null)
         {
                StringTokenizer st = new StringTokenizer(data);
                String s = st.nextToken()+ "  "; //ID
                  
                double mfe = Double.parseDouble(st.nextToken());

                int len = Integer.parseInt(st.nextToken());   
       
                int loops = Integer.parseInt(st.nextToken());                         
                    
                double stems = Double.parseDouble(st.nextToken());

                double Nmfe = mfe/len;

                double mfe2 = Nmfe/stems; 
                double mfe3 = Nmfe/loops;           

                p.println(s+"   "+mfe2+"  "+mfe3);
                
                //System.out.println(s);
          }
            

          p.close();
          r1.close();

      }

      catch(Exception e)
      {
         System.out.println(e.toString());
       }
   }

}    







