package programs;


import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.SecondaryStructure;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.parameters.AtomType;
import meshi.util.*;
import meshi.util.mysql.*;
import meshi.energy.rapdf.*;
import meshi.parameters.ResidueType;

import java.io.*;
import java.util.*;
import java.text.*;

/***
    Data: 06/03/2011
    Creator: N. Vanetik
    This programs scans through the given PDB database (all *.pdb files in the given directory)
    and moves files with <= (number of residues) into the target directory.
    The program receives following parameters:
        <source directory> <target directory> <number of residues> 
***/
public class PdbDatabaseSeparation{

    public static final String PDB_EXTENSION = ".pdb";
    public static final  String PDB_SEQRES_RECORD = "SEQRES";

    // constructor that performs the separation
    public PdbDatabaseSeparation(String srcDir, String trgDir, int numOfResidues)
    {
        File src, trg;
        src = new File(srcDir);
        trg = new File(trgDir);

	if (!src.exists() || !src.isDirectory())
	    {
		System.out.println(srcDir+" does not exists or is not a directory");
		System.exit(0);
	    }
       if (trg.exists() && !trg.isDirectory())
	    {
		System.out.println(trgDir+" is not a directory");
		System.exit(0);
	    }

       if (!trg.exists())
	   {
	       trg.mkdir();
               System.out.println("Created directory "+trgDir);
           }
       File srcFiles[] = src.listFiles();
       String fasta;
       System.out.println("Moving files from "+srcDir+" to "+trgDir);
       for (int i=0; i<srcFiles.length; i++)
	   if (srcFiles[i].getName().endsWith(PDB_EXTENSION))
	   {
	       fasta = readPdb(srcFiles[i]);
	       if (fasta.length()<=numOfResidues)
		   {
                      boolean success = srcFiles[i].renameTo(new File(trgDir, srcFiles[i].getName()));
                      if (!success) {
			  System.out.println("Error moving file "+srcFiles[i].getName());
                        }
		      else
                          System.out.println("Moved file "+srcFiles[i].getName()+" with "+fasta.length()+" residues");
		   }
	   }
    }

    /***
        Read FASTA sequence from a pdb file
    ***/
    private String readPdb(File toRead)
        {
            try
            {
                int lineCounter = 0, j;
                String line, token, oneLetterCode;
                StringTokenizer st;
                StringBuffer fasta = new StringBuffer();
                int total = 0;
                LineNumberReader currentLineReader = new LineNumberReader(new FileReader(toRead));
                
                while ((line = currentLineReader.readLine()) != null)
                        {
                            /*if (lineCounter==0)
                            {
                                // write header
                                fasta.append(FASTA_COMMENT+line.toLowerCase()+"\n");
                                }*/
                                
                            if (line.startsWith(PDB_SEQRES_RECORD))
                            {
                                st = new StringTokenizer(line, " \n\r\t");
                                //System.out.println("=========> Reading line "+lineCounter);
                                
                                // skip SEQRES <record num> type <amount> data
                                for (j=0; j<5 && st.hasMoreTokens(); j++)
                                    st.nextToken();

                                if (j<5)
                                    break;
                                // read 3-letter residue code and translate into 1-letter residue code
                                
                                while (st.hasMoreTokens())
                                {
                                    token = st.nextToken(); total++;
                                    oneLetterCode = ResidueType.typeOneLetter(token);
                                    if (oneLetterCode==null)
                                        continue;
                                    //System.out.println("=========> Found residue #"+total+" "+token+" translated to "
                                    //                   +oneLetterCode);
                                    fasta.append(oneLetterCode);
                                }
                            }
                            lineCounter++;
                        }

                        currentLineReader.close();
                        currentLineReader = null;
                        fasta.append("\n");
                        return fasta.toString();
            }
            catch (IOException e1)
                {
                   System.out.println("Exception "+e1.toString()); 
                }
            return null;
        }
    
 

    /******************************************************************/
    /*                          THE MAIN METHOD                       */
    /******************************************************************/
    public static void main(String[] args) throws IOException{
        
        if (args==null || args.length<3)
        {
            System.out.println("FORMAT: PdbDatabaseSeparation <source directory> <target directory> <number of residues>");
            System.exit(1);
        } 
        try {
	    PdbDatabaseSeparation rmc = new PdbDatabaseSeparation(args[0], args[1], Integer.parseInt(args[2]));
	}
	catch (NumberFormatException e)
	    {
		e.printStackTrace();
	    }
    }
}