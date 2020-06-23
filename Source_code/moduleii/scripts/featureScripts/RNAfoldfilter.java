//Filters the RNAfold-related features

import java.io.*;
import java.util.*;

class RNAfoldfilter
{

   public static void main(String args[])
   {
       try 
        {
           BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".RNAfold1")));
           PrintWriter p1 = new PrintWriter("../data/"+args[0]+".RNAfold2");
 
           String data = null;

           int i = 0;
                                 
           while((data=r1.readLine()) != null)
           {
              
              if(data.indexOf('>') >= 0)
              {
                i++;
                p1.print(i+"  ");
				//System.out.println(data); 
		  } 
 
              		  
		  else if((data.indexOf(',') >= 0 || data.indexOf('}') >= 0 || data.indexOf('{') >= 0 || data.indexOf('(') >= 0 || data.indexOf(')') >= 0) && data.indexOf('=') < 0)
              {
                
                StringTokenizer st = new StringTokenizer(data);
                st.nextToken(); 
                String a = st.nextToken();
                String b = "";
                if (a.length() > 1)
                   b = a;
                else
                   b = st.nextToken();
                 
                int x1 = b.indexOf('-');
                              
                String c = b.substring(x1,(b.length()-1));                   
                p1.print(c+"  ");
		    //System.out.println(c);
              } 
              else if(data.indexOf("=") > 0)
              {   }

              else if(data.indexOf("frequency") >= 0)
              {
                StringTokenizer st = new StringTokenizer(data);
                for(int j=1; j<=6; j++)
                  st.nextToken(); 
                String d1 = st.nextToken();
                st.nextToken(); st.nextToken();
                String d2 = st.nextToken();
                StringTokenizer st1 = new StringTokenizer(d1,";");

                String freq = st1.nextToken();
                double div = Double.parseDouble(d2)/2;

	          p1.println(freq+"  "+div); 
              }                           
              
           }
          
           //System.out.println("i ="+i);
           r1.close(); p1.close();   
           
        }

        catch(Exception e)
        {
           e.printStackTrace();
         }   
      }
} 

                   
