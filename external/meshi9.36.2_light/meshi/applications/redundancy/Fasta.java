package meshi.applications.redundancy;
import java.io.*;
import java.util.*;

/***
    Class for static methods handling fasta sequences
 ***/
public class Fasta {
    public static final int    DEBUG_LEVEL = 1; // debug printouts
    public static final int    PRINT_LINE_LENGTH = 80;
    
    /*******************************************************/
    /***** ***** ***** AUXILIARY FUNCTIONS ***** ***** *****/
    /*******************************************************/
    // Print fasta sequence to the screen
    public static void printFasta(String fasta, String header)
        {
            System.out.println(header);
            printFasta(fasta); 
        }
        
    public static void printFasta(String fasta)
        {
            //System.out.println("Reading FASTA sequence ===>");
            int i = 0;

            for (i=0; i<fasta.length(); i+=PRINT_LINE_LENGTH)
            {
                System.out.println(fasta.substring(i,Math.min(fasta.length(), i+PRINT_LINE_LENGTH)));
            }
        }

    // check if the sequence is in FASTA format
    public static boolean isFasta(String str)
        {
            byte[] b = str.getBytes();

            for (int i=0; i<b.length; i++)
                if (!(b[i]>='A' && b[i]<='Z'))
                    return false;

            if (DEBUG_LEVEL>=2)
                System.out.println("TESTED => isFasta=true for "+str);
            
            return true;
        }

    /***
        Clean the fragment from non A-Z characters - for querying the database
     ***/
    public static String cleanFasta(String fasta)
        {
            StringBuffer sb = new StringBuffer();
            int len = fasta.length();
            byte[] b = fasta.getBytes();

            for (int i=0; i<len; i++)
                if (b[i]>='A' && b[i]<='Z')
                    sb.append(fasta.substring(i,i+1));
            
            return sb.toString();
        }
    
}
