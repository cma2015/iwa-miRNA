//Filters all the features calculated into a one file and the best subset of fetures into another file.


import java.io.*;
import java.util.*;

class filter_all
{

  public static void main(String args[])
  {
     try
     {
         BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".data1")));
         BufferedReader r2 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".data2")));
         BufferedReader r3 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".zdata")));
         BufferedReader r4 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".mfe14")));
         BufferedReader r5 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".mfe23")));
         BufferedReader r6 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".bp1")));
         BufferedReader r7 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".bp2")));
         BufferedReader r8 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".mfold2")));
         BufferedReader r9 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".RNAfold2")));

         PrintWriter p1 = new PrintWriter("../data/all."+args[0]+"-48.features");
         PrintWriter p2 = new PrintWriter("../data/selected."+args[0]+".-21.features");

         r1.readLine();
         r2.readLine();
         r3.readLine();
         r8.readLine();

         String d1 = null;
     
         int j = 0;

         while((d1=r1.readLine()) != null)
         {
            j++;
            //.data1 file
            StringTokenizer t1 = new StringTokenizer(d1);
            for(int i=1; i<=28; i++)
               t1.nextToken();
            String cg = t1.nextToken();
            
            t1.nextToken();

            String di = "";
            for(int i=1; i<=16; i++)
              di += t1.nextToken()+"  ";
            
            t1.nextToken();
            String dP = t1.nextToken();  
            t1.nextToken();
            String dG = t1.nextToken();
            t1.nextToken();
            String dQ = t1.nextToken();
            t1.nextToken();
            String dD = t1.nextToken();

            //.data2 file
            StringTokenizer t2 = new StringTokenizer(r2.readLine());
            t2.nextToken(); t2.nextToken();
            String len = t2.nextToken();
            for(int i=1; i<=5; i++)
                t2.nextToken();
            String dF = t2.nextToken(); 

            //.zdata file
            StringTokenizer t3 = new StringTokenizer(r3.readLine());
            for(int i=1; i<=4; i++)
               t3.nextToken(); 
            String zG = t3.nextToken();
            for(int i=1; i<=5; i++)
               t3.nextToken(); 
            String zQ = t3.nextToken();
            for(int i=1; i<=5; i++)
               t3.nextToken(); 
            String zD = t3.nextToken(); 
            for(int i=1; i<=5; i++)
               t3.nextToken(); 
            String zP = t3.nextToken(); 
            for(int i=1; i<=5; i++)
               t3.nextToken(); 
            String zF = t3.nextToken();

            //.mfe14 file
            StringTokenizer t4 = new StringTokenizer(r4.readLine());
            t4.nextToken();
            String mfe1 = t4.nextToken();
            String mfe4 = t4.nextToken();

            //.mfe23 file
            StringTokenizer t5 = new StringTokenizer(r5.readLine());
            t5.nextToken();
            String mfe2 = t5.nextToken();
            String mfe3 = t5.nextToken();

            //.bp1 file
            StringTokenizer t6 = new StringTokenizer(r6.readLine());
            for(int i=1; i<=9; i++)
               t6.nextToken(); 
            String auL = t6.nextToken(); // au/L
            String gcL = t6.nextToken();
            String guL = t6.nextToken(); 
            String pau = t6.nextToken(); // %au/L
            String pgc = t6.nextToken();
            String pgu = t6.nextToken(); 

            //.bp2 file
            StringTokenizer t7 = new StringTokenizer(r7.readLine());
            t7.nextToken(); t7.nextToken();
            String bpStems = t7.nextToken();
            String auStems = t7.nextToken();
            String gcStems = t7.nextToken();
            String guStems = t7.nextToken();

            //.mfold2 file
            StringTokenizer t8 = new StringTokenizer(r8.readLine());
            t8.nextToken(); 
            String dH = t8.nextToken();
            String dS = t8.nextToken();
            String Tm = t8.nextToken();
            double dHL = Double.parseDouble(dH)/Integer.parseInt(len);
            double dSL = Double.parseDouble(dS)/Integer.parseInt(len);
            double TmL = Double.parseDouble(Tm)/Integer.parseInt(len);

            //.RNAfold-based
            StringTokenizer t9 = new StringTokenizer(r9.readLine());
            t9.nextToken(); 
            String mfe = t9.nextToken();
            String efe = t9.nextToken();   
            String freq = t9.nextToken();   
            String div = t9.nextToken();

            double diff = (Double.parseDouble(mfe)-Double.parseDouble(efe))/Integer.parseInt(len);
            double nefe = Double.parseDouble(efe)/Integer.parseInt(len);
                          
            String s = cg+"  "+di+"  "+mfe1+"  "+mfe2+"  "+dG+"  "+dP+"  "+dQ+"  "+dD+"  "+dF;
                   s += "  "+zG+"  "+zP+"  "+zQ+"  "+zD+"  "+zF+"  "+mfe3+"  "+mfe4+"  "+nefe+"  "+freq+"  "+div+"  "+diff;
                   s += "  "+dH+"  "+dHL+"  "+dS+"  "+dSL+"  "+Tm+"  "+TmL+"  "+auL+"  "+gcL+"  "+guL;
                   s += "  "+bpStems+"  "+auStems+"  "+gcStems+"  "+guStems;

            p1.println(s);

 
            String selected = cg+"  "+mfe1+"  "+mfe2+"  "+dG+"  "+dQ+"  "+dF+"  "+zD+"  "+nefe+"  "+diff;
                   selected += "  "+dS+"  "+dSL+"  "+bpStems+"  "+auL+"  "+gcL+"  "+guL+"  "+auStems+"  "+gcStems+"  "+guStems;
                   selected += "  "+mfe3+"  "+mfe4+"  "+div;     

            p2.println(selected);


        }
      
        r1.close();r2.close();r3.close();r4.close();r5.close();r6.close();r7.close();r8.close();r9.close();
        p1.close();
        p2.close();
     
     }

     catch(Exception e)
     {
         e.printStackTrace();
     }
  }

}

                      
   


