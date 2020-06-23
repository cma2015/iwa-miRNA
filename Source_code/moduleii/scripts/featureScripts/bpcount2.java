//calculates the base pairs related fetaures.

import java.io.*;
import java.util.*;

class bpcount2
{

   public static void main(String args[])
   {
       try 
        {
           BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".bp1")));
           BufferedReader r2 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".data2")));

           PrintWriter p1 = new PrintWriter("../data/"+args[0]+".bp2");

           r2.readLine(); 

           String data = null;
     
           int j = 0;

           while((data=r1.readLine()) != null)
           {
              StringTokenizer st = new StringTokenizer(data);
              String id = st.nextToken(); 
              st.nextToken(); 
              int tot = Integer.parseInt(st.nextToken());
              st.nextToken(); st.nextToken(); st.nextToken();
              double pau = Double.parseDouble(st.nextToken());
              double pcg = Double.parseDouble(st.nextToken()); 
              double pgu = Double.parseDouble(st.nextToken());               
              
              String data2 = r2.readLine();
              StringTokenizer t = new StringTokenizer(data2);
              t.nextToken(); t.nextToken(); t.nextToken(); t.nextToken(); 
              int stems = Integer.parseInt(t.nextToken());

              String s = id+"  "+stems+"  "+(1.0*tot/stems)+"  "+(pau/stems)+"  "+(pcg/stems)+"  "+(pgu/stems);

              p1.println(s);

              j++;

           }

          r1.close(); r2.close(); p1.close();

          //System.out.println(j);


         }

         catch(Exception e)
         {
           e.printStackTrace();
         }
     }
}   
              













                
                 