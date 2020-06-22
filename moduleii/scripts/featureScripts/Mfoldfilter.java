//Filters Mfold related features 

import java.io.*;
import java.util.*;

class Mfoldfilter
{

   public static void main(String args[])
   {
       try 
        {
           BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".mfold")));
           PrintWriter p1 = new PrintWriter("../data/"+args[0]+".mfold2");
 
           String data = null;

           int i = 0;

           boolean f = false;
                                 
           while((data=r1.readLine()) != null)
           {
               
              if(data.indexOf("dG	dH	dS	Tm") >= 0)
                 f = true;
              if(f == true)
              { 
                 i++;
                 p1.println(data);
              } 
        
           }
      
           //System.out.println("i = "+i);
           r1.close(); p1.close();                
        }

        catch(Exception e)
        {
           e.printStackTrace();
         }   
      }
} 

                   
