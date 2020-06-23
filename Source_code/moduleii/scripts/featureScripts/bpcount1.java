//calculates the base pairs related fearures


import java.io.*;
import java.util.*;

class bpcount1
{

   public static void main(String args[])
   {
       try 
        {
           BufferedReader r1 = new BufferedReader(new InputStreamReader(new FileInputStream("../data/"+args[0]+".fold")));
           PrintWriter p1 = new PrintWriter("../data/"+args[0]+".bp1");
 
           String data = null;
     
           int j = 0;

           while((data=r1.readLine()) != null)
           {
             if(data.indexOf('>') >= 0)
             {
               j++;
                  
               //System.out.println(j+"");

               String s = r1.readLine();
               StringTokenizer st = new StringTokenizer(r1.readLine());
               String f = st.nextToken();

               char[] seq = s.toCharArray();
               char[] fold = f.toCharArray();

               //System.out.println(s);
               //System.out.println(f+"\n");
              
               int au = 0;
               int gc = 0;
               int gu = 0;
 
               Vector stk = new Vector();
               Vector qu = new Vector(); 
                  
               int i= 0;
               int tot = 0;

               while(i < fold.length)
               {
 
                  while(fold[i] == '(')
                  {
                    stk.add(0,new Character(seq[i]));
                    i++;
                    tot++;
                  }
                   
                  if(fold[i] == ')')
                  {
                     char c1 = ((Character)stk.remove(0)).charValue();
                     char c2 = seq[i];
                     //System.out.println(c1+" - "+c2);
                 
                     if((c1 == 'A' && c2 == 'U') || (c1 == 'U' && c2 == 'A') || (c1 == 'a' && c2 == 'u') || (c1 == 'u' && c2 == 'a'))
                       au++;
                     if((c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G') || (c1 == 'g' && c2 == 'c') || (c1 == 'c' && c2 == 'g'))
                       gc++;
                     if((c1 == 'G' && c2 == 'U') || (c1 == 'U' && c2 == 'G') || (c1 == 'g' && c2 == 'u') || (c1 == 'u' && c2 == 'g'))
                       gu++;

                     i++;
                  }
           
                  else if(fold[i] == '.')
                   i++;
 
                }                     
  
                String s1 = j+"   "+seq.length+"    "+tot+"    "+au+"     "+gc+"     "+gu+"     "+(au*100.0/tot)+"     "+(gc*100.0/tot)+"     "+(gu*100.0/tot)+
                              "   "+(1.0*au/seq.length)+"   "+(1.0*gc/seq.length)+"  "+(1.0*gu/seq.length)+"  "+
                                   (au*100.0/(tot))+"  "+(gc*100.0/(tot))+"  "+(gu*100.0/(tot));                          
                  
                p1.println(s1);


            }

           } 
            r1.close(); p1.close();

          }

          catch(Exception e)
          {
             e.printStackTrace();
           }
      }
}  
          