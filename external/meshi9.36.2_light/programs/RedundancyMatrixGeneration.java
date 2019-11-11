package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.SecondaryStructure;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.parameters.AtomType;
import meshi.applications.redundancy.*;
import meshi.util.mysql.*;

import java.io.*;
import java.util.*;

/***
    Data: 26/12/2010
    Creator: N. Vanetik
    This programs computes the redundancy score of each residue in a protein for a given database. T
    his matrix is stored in a file.
    The program receives following parameters:
        <directory that holds PDB files> <full path to blast> <alignment lower bound 1..100> <evalue lower bound 0..1> <output file - optional>
***/
public class RedundancyMatrixGeneration {
    public static final int    DEBUG_LEVEL = 1; // debug printouts
    public static final int    SAVE_FRAGMENT_BOUND  = 0;
    public static final int    MAX_DOUBLES_ON_LINE = 6;
    public static final String FASTA_FORMAT = "fasta";
    public static final String PDB_FORMAT = "pdb";
    
    public static final int    STORE_IN_MEMORY = 0;
    public static final int    STORE_IN_FILE   = 1;
    public static final int    STORE_IN_DATABASE  = 2;
    private int                whereToStoreResults = STORE_IN_DATABASE;

    public static final int    NO_FIT = -1;
    
    public static final boolean STORE_MATRIX_IN_MEMORY = false;
    
    private HashMap         fragmentCountMap;            // hash map to store fragments already found
    private String          dbName          = null;          // database to derive BLAST queries from
    private String          blastDbName     = null;          // database to run BLAST reset
    private String          matrixFileName  = null;          // the file to store the results in
    private FileOutputStream matrixFile     = null;
    private int             currentResultNumber = 0; // the number of current file for which the redundancy measure is computed 
    
    // This list holds redunancy weights for each sequence in the database
    private ArrayList redunancyMatrix = new ArrayList();
    
    /*************************************************************/
    /****************     BLAST VARIABLES       ******************/
    /*************************************************************/
    private static final int DEFAULT_FRAGMENT_LENGTH = 100;
    private static final int DEFAULT_ALIGNMENT_BOUND = 100;
    private String           TEMP_QUERY_FILE_NAME = "tmpquery.txt";
    private String           TEMP_RESULTS_FILE_NAME = "tmpresults.txt";
    private static final int BLAST_OUTPUT_FORMAT = 7;
    private static final int DEFAULT_EVALUE = 1;

    // username, password and DB name for mysql access
    private   String filesTable     = null;
    private   String fragmentsTable = null;
    private   String redundancyTable = null;

    private static final String FRAGMENTS_TABLE_ATTRIBUTES = "(fragment varchar(1000), cnt double)";
    private static final String FILES_TABLE_ATTRIBUTES = "(filename varchar(1000))";
    private static final String REDUNDANCY_TABLE_ATTRIBUTES = "(filename varchar(100), residue integer, redu double)";

    private int fragmentLength = -1;
    private int alignmentBound = -1;
    private double highestEvalue    = -1;
    private String blastExecutable  = null;
    private boolean useMysql = true;

    private MysqlHandler mysql_updater = null;
    private DataReader dr = null;

    private CommandList commands = null;

    /*** constructor ***/
    public RedundancyMatrixGeneration(CommandList commands)
        {
            fragmentCountMap = new HashMap();
            this.mysql_updater = new MysqlHandler(commands);
            this.commands = commands;

            this.dbName = commands.firstWord("pdbdatabase").secondWord();
            System.out.println("Found dbname="+this.dbName);
            this.blastExecutable = commands.firstWord("blastpath").secondWord();
            System.out.println("Found blast path="+this.blastExecutable);
            this.blastDbName = commands.firstWord("blastdatabase").secondWord();
            System.out.println("Found blast database="+this.blastDbName);

            try
            {
              this.fragmentLength = Integer.parseInt(commands.firstWord("fragmentlength").secondWord());
              if (fragmentLength <= 0)
                   fragmentLength = DEFAULT_FRAGMENT_LENGTH;
              System.out.println("Found param="+fragmentLength);
              this.alignmentBound = Integer.parseInt(commands.firstWord("fragmentalignment").secondWord());
              if (alignmentBound <= 0)
                  alignmentBound = DEFAULT_ALIGNMENT_BOUND;
              System.out.println("Found param="+alignmentBound);
              this.highestEvalue = Double.parseDouble(commands.firstWord("evaluebound").secondWord());
              if (highestEvalue<=0)
                  highestEvalue = DEFAULT_EVALUE;
              System.out.println("Found param="+highestEvalue);
              matrixFileName = commands.firstWord("redudancyfile").secondWord();
              
              TEMP_QUERY_FILE_NAME = commands.firstWord("blastqueryfile").secondWord();
              TEMP_RESULTS_FILE_NAME = commands.firstWord("blastresultsfile").secondWord();
              System.out.println("Found param="+TEMP_QUERY_FILE_NAME);
            
              if (blastExecutable == null)
              {
                System.out.println("BLAST executable path is undefined, use -blast <path> ...");
                System.exit(0);
              }

              if (DEBUG_LEVEL>=2)
                  System.out.println("Computing redundancy for fragment length "
                                     +fragmentLength+" alignment bound "+alignmentBound+" evalue "+highestEvalue);
              /*** create DB files ***/
              boolean rc = false;
              if (commands.keyExists("redundancytable"))
              {
                  this.redundancyTable = commands.firstWord("redundancytable").secondWord();
                  rc |= mysql_updater.createNonExistingTable(this.redundancyTable,
                      REDUNDANCY_TABLE_ATTRIBUTES);
              }
              if (commands.keyExists("filestable"))
              {
                  this.filesTable = commands.firstWord("filestable").secondWord();
                  rc |= mysql_updater.createNonExistingTable(this.filesTable,
                      FILES_TABLE_ATTRIBUTES);
              }
              if (commands.keyExists("fragmentstable"))
              {
                  this.fragmentsTable = commands.firstWord("fragmentstable").secondWord();
                  rc |= mysql_updater.createNonExistingTable(this.fragmentsTable,
                      FRAGMENTS_TABLE_ATTRIBUTES);
              }
              
              openMatrixFile();
              computeMatrix();
              mysql_updater.close();
              closeMatrixFile();
            }
            catch (NumberFormatException e)
            {
                e.printStackTrace();
                System.exit(0);
            }            
        }

    /************************************************************************/
    /*****************  METHODS TO HANDLE THE RESULTS FILE ******************/
    /************************************************************************/
    private void openMatrixFile()
        {
            try
            {
            if (matrixFileName != null)
                matrixFile = new FileOutputStream(matrixFileName, true);
            }
            catch (IOException e)
            {
                e.printStackTrace();
            }
        }

    private void closeMatrixFile()
        {
            try
            {
            if (matrixFile != null)
                matrixFile.close();
            }
            catch (IOException e)
            {
                e.printStackTrace();
            }
        }


    /************************************************************************/
    /*********************  MATRIX COMPUTATION METHODS **********************/
    /************************************************************************/
    
    /***
        Compute redundancy measure for the query file(s) vs given database
     ***/
    public void computeMatrix()
        {
            // read all fasta lines from a file and compute their redundancy vs the database
            this.dr = new DataReader(this.dbName, PDB_FORMAT);
            String fasta, currentFileName;

            // get all FASTA sequences from the data reader
            redunancyMatrix.clear();
            
            while ((fasta = dr.getNext()) != null)
            {
                if (DEBUG_LEVEL>=2)
                    System.out.println("=====> read fasta = "+fasta);
                currentFileName = dr.getCurrentFileName();
                
                if (DEBUG_LEVEL>=1)
                    System.out.print("=====> checking for file = "+currentFileName);

                int rc = mysql_updater.selectCount(this.filesTable, "filename='"+currentFileName+"'");
                if (rc<=0)
                {
                   if (DEBUG_LEVEL>=1)
                    System.out.println("... reading file = "+currentFileName); 
                   computeRedundancy(fasta);
                   mysql_updater.insertIntoTable(this.filesTable, "(filename)", "('"+currentFileName+"')");
                }
                else
                {
                    if (DEBUG_LEVEL>=1)
                       System.out.println();
                }
            }
        }

    /***
        A getter
     ***/
    public ArrayList getMatrix() { return redunancyMatrix; }

    /***
        Print the entire matrix
     ***/
    public void printMatrix()
        {
            double [] W;
            System.out.println("=========> redundancy matrix:");
            for (int i=0; i<redunancyMatrix.size(); i++)
            {
                W = (double[])redunancyMatrix.get(i);
                System.out.print("r["+i+"]=");
                for (int j=0; j<W.length; j++)
                {
                    System.out.print(" "+W[j]);
                    if (j > 0 && j % MAX_DOUBLES_ON_LINE==0)
                      System.out.println();  
                }
                System.out.println();
            }
        }

    /*********************************************************************************/
    /****************         HANDLING REDUNDANT FRAGMENTS      **********************/
    /*********************************************************************************/
    private void storeFragment(String fragment, int count)
        {
            if (DEBUG_LEVEL>=3)
                System.out.println("+++++++++> storing fragment "+fragment+" with count "+count); 
            if (useMysql)
            {
                String parameters = "(fragment, cnt)";
                String values = "('"+fragment+"',"+count+")";
                mysql_updater.insertIntoTable(this.fragmentsTable, parameters, values);
            }
            else
              fragmentCountMap.put(fragment, new Integer(count));
                
        }

    private int countOfFragment(String fragment)
        {
            if (DEBUG_LEVEL>=3)
                System.out.println("+++++++++> looking for fragment "+fragment); 
            if (useMysql)
            {
                //int count = mysql_updater.find(fragment);
               String whereCondition = "fragment = '"+Fasta.cleanFasta(fragment)+"'";
               int count = (int)mysql_updater.selectDouble(this.fragmentsTable, whereCondition, "cnt");
               return count;
            }
            else
            {
                Object count;
                if ((count = fragmentCountMap.get(fragment)) != null)
                    return ((Integer)count).intValue();
                else
                    return -1;
            }
            //return -1;
        }

    
    /*********************************************************************************/
    /****************                 COMPUTE REDUNDANCY        **********************/
    /*********************************************************************************/
    /***
        Match given fasta sequence with blast and store results in the local matrix
     ***/
    public void computeRedundancy(String fasta)
        {
            int len = fasta.length();
            File queryTempFile;
            FileOutputStream fos  = null;
            Process p, p2;
            String commandString, line, token, fragment;
            StringTokenizer st;
            InputStream in;
            int fit = NO_FIT;
            double alignment, evalue;
            LineNumberReader lnr;
            Object storedFit;
            int effectiveLength;
            int[] weightedRC = new int[len]; // array to collect the RC of each residue
            int[] weights    = new int[len]; // array to collect the weights
            int factor = fragmentLength/2;
            int weight;
            double logweight;
            int previousFit = NO_FIT;

            Arrays.fill(weightedRC, 0);
            Arrays.fill(weights, 0);
            
            // pre-build command string parts
            commandString = this.blastExecutable+" -db "+this.blastDbName+" -query "+TEMP_QUERY_FILE_NAME
                +" -evalue "+highestEvalue+" -outfmt "+BLAST_OUTPUT_FORMAT+" -out "+TEMP_RESULTS_FILE_NAME;
            //countString = "wc -l "+TEMP_RESULTS_FILE_NAME;

            if (DEBUG_LEVEL>=3)
                System.out.println("++++++++++++> fasta="+fasta);

            queryTempFile = new File(TEMP_QUERY_FILE_NAME);
            try
            {
               if (!queryTempFile.exists())
                  queryTempFile.createNewFile();
            }
            catch (IOException e)
                {
                    e.printStackTrace();
                }

            boolean mustCallBlast = true;
            
            for (int i=0; i<len; i++)
            {
                if (i+fragmentLength>=len && i>0) // too-short fragments - stop unless this is the entire sequence
                    break;

                fragment = fasta.substring(i, Math.min(i+fragmentLength, len));
                if (DEBUG_LEVEL>=3)
                {
                    Fasta.printFasta(fragment, "Matching fragment ===>");
                    System.out.println("++++++++++++> fragment="+fragment);
                }
                // check if we have found this fragment already
                try
                {
                    previousFit = fit;    
                    fit = countOfFragment(fragment);

                    /*** with high probability, the fragments that fitted [x,x+100]
                         also fit fragment [x+1,x+101] for low evalues ***/
                    if (fit<0 && previousFit != NO_FIT && !mustCallBlast)
                    {
                        fit = previousFit;
                        previousFit = NO_FIT;
                        mustCallBlast = true;
                        if (DEBUG_LEVEL>=3)
                            System.out.println("++++++++++++> saving blast call for i="+i);
                    }
                    
                    if (DEBUG_LEVEL>=3)
                      System.out.println("++++++++++++> found fit="+fit);

                    /*** otherwise, use blast to compute the fit ***/
                    if (fit<0) // actually match the fragment
                    {
                     mustCallBlast = false;   
                                     // store the fragment in a temporary file
                     fos = new FileOutputStream(queryTempFile);
                     fos.write(fragment.getBytes());
                     fos.close();

                     // call BLAST as a new process
                     if (DEBUG_LEVEL>=3)
                         System.out.println("++++++++++++> "+commandString);
                     p = Runtime.getRuntime().exec(commandString);

                     // wait for it to end
                     p.waitFor();
                     
                     // read results from a temporary file
                     //p2 = Runtime.getRuntime().exec(countString);
                     //in = p2.getInputStream();
                     //p2.waitFor();
                     

                     // count suitable matches
                     // FileReader fl = new FileReader(TEMP_RESULTS_FILE_NAME);
                     lnr = new LineNumberReader(new FileReader(TEMP_RESULTS_FILE_NAME));
                     fit = 0; int j; int bound;
                     while ((line = lnr.readLine()) != null)
                      {
                          if (line.startsWith("#"))
                              continue;
                          
                          st = new StringTokenizer(line, " \t");
                          j=0; token = "";

                          //System.out.println("================> RESULT line = "+line);
                          /***
                              How the output line looks like:
                              ---------------------------------------------------------------------------------------------
                              # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start,
                              q. end, s. start, s. end, evalue, bit score
                           ***/
                          bound = 3; 
                          while (j<bound)
                               {
                                 if (st.hasMoreTokens())
                                     token = st.nextToken();
                                 j++;
                               }

                          /*** find alignment ***/
                          alignment = Double.parseDouble(token);
                          j++;

                          /*** find evalue ***/
                          bound = 12; 
                          while (j<bound)
                               {
                                 if (st.hasMoreTokens())
                                     token = st.nextToken();
                                 j++;
                               }
                          evalue = Double.parseDouble(token);
                          if (DEBUG_LEVEL>=3)
                             System.out.println("++++++++++++> Found alignment="+alignment+" evalue= "+evalue);

                          if (evalue<= highestEvalue && alignment >= (double)alignmentBound)
                              fit++;
                          
                      }

                     lnr.close();
                     //fl.close();

                     if (DEBUG_LEVEL>=3 && fit!=0)
                      System.out.println("++++++++++++> Good matches "+fit); // these are all matches
                     // save the fragment and its count in the hashmap
                     if (fit >= SAVE_FRAGMENT_BOUND)
                     {
                         storeFragment(fragment, fit);
                     }                     
                }
                     // store results: redundancy measure of a segment and of each residue in a segment w.r.t.
                     // the given sequence
                     // update the redundancy count of each residue in range i..i+fragmentLength-1
                     effectiveLength = Math.min(i+fragmentLength-1, len-1);
                     for (int j=i; j<effectiveLength; j++)
                       {
                           // weight of the RC depends reset a location of j within i..i+fragmentLength-1
                           if (fit != 0)
                           {
                              weight = Math.min(j-i, effectiveLength-j); // could be log or linear
                              logweight = Math.log(weight);
                              if (weight != 0)
                              {
                                weights[j] += weight;   
                                weightedRC[j] += weight*fit;
                                if (DEBUG_LEVEL>=3)
                                   System.out.println("++++++++++++> weight : "+weight);
                              }
                           }
                           // the final score will be weightedRC[j]/weights[j]
                        }
                     
                }
                catch (IOException e1)
                {
                   System.out.println("Exception "+e1.toString()
                                   +" while matching FASTA sequence ");
                   e1.printStackTrace();
                }
                catch (InterruptedException e2)
                {
                    System.out.println("Exception "+e2.toString());
                }
            }
            // now, all RCs have been computed
            // store them in a file in format residue Rcount|residue Rcount|...|residue Rcount\n
            storeResults(fasta, weightedRC, weights); 
            whereToStoreResults = STORE_IN_FILE;
            storeResults(fasta, weightedRC, weights);          
        }


    /***
        Store the result of redundancy count in a file
        Parameters: computed redundancy count, source string and optional output file (if null, ignored)
        Output: stored in a global matrix matrix[structure_number][residue instance number][measure]
    ***/
    public void storeResults(String fasta, int[] weightedRC, int[] weights)
        {
            double score_i, score_j;
            double W1[];
            double W2[];

            W1 = new double[weightedRC.length]; // weight and log weight
            
            if (DEBUG_LEVEL>=2)
                System.out.println("++++++++++++> Storing results for sequence of length "+fasta.length());

            for (int i=0; i<weightedRC.length; i++)
                 {
                   if (weightedRC[i] != 0)
                       W1[i] = (double)weights[i]/(double)weightedRC[i]; // average redundancy score
                   else
                       W1[i] = 1; // default - redundancy takes no effect
                   if (DEBUG_LEVEL>=2 && W1[i]!=1)
                       System.out.println("final rweight["+i+"]="+W1[i]);
                  }
                
                   storeResults(W1);
                
                if (DEBUG_LEVEL>=2)
                    System.out.println("++++++++++++> storing atom scores... ");
                
        }

/***
    Store current redundancy scores in memory, file or database table
 ***/
private void storeResults(double W1[])
{
    switch (whereToStoreResults)
    {
    case STORE_IN_MEMORY: 
        this.redunancyMatrix.add(W1);
        break;

    case STORE_IN_FILE:
        if (matrixFile != null)
        {
                  if (DEBUG_LEVEL>=1)
		      System.out.print("=====> storing results in file "+matrixFileName+" ");

                  StringBuffer result = new StringBuffer();
                  result.append(currentResultNumber+">");
                  
                  for (int j=0; j<W1.length; j++)
                       result.append(" "+W1[j]);
                   result.append("\n");
                   
                  currentResultNumber++;
                  try
                  {
                     matrixFile.write(result.toString().getBytes());
                     matrixFile.flush();
                  }
                  catch (IOException e)
                  {
                      e.printStackTrace();
                  } 
        }
        break;

    case STORE_IN_DATABASE:
        //mysql_updater.storeRedundancy(dr.getCurrentFileName(), W1);
        String parameters = "(filename, residue, redu)";
        String values = "('"+dr.getCurrentFileName()+"',";
        for (int i=0; i<W1.length; i++)
        {
            mysql_updater.insertIntoTable(this.redundancyTable, parameters, values+i+","+W1[i]+")");
        }
        break;
        
    default:
        // no place to store results is defined
        break;
    }
     
}

    /***
        Store redundancy weight in the following format:
        for each sequence, write a line:
        "line number 0.."> <redundancy weight> .. <redundancy weight> 
     ***/
    public void storeMatrix(String fileName)
        {
            StringBuffer result = null;
            double W1[];
            
            try {
                FileOutputStream fos = new FileOutputStream(fileName);
                for (int i=0; i<this.redunancyMatrix.size(); i++)
                {
                   result = new StringBuffer();
                   result.append(i+">");
                   W1 = (double[])this.redunancyMatrix.get(i);
                   for (int j=0; j<W1.length; j++)
                       result.append(" "+W1[j]);
                   result.append("\n");
                   fos.write(result.toString().getBytes());
                }
                fos.close();
            }
            catch(IOException e)
            {
                System.out.println("Error writing to file "+fileName);
                e.printStackTrace();
            }
        }
 
    
    /*******************************************************************************************/
    /*************************  THE MAIN PROGRAM TO BE CALLED   ********************************/
    /*******************************************************************************************/
    public static void main(String[] args) throws IOException{
        CommandList commands = Utils.init(args, 1, "RedundancyMatrixGeneration");
        
        //if (args==null || args.length<4)
        // {
        //    System.out.println("FORMAT: RedudnancyMatrixGeneration <database directory>"+
        //                       "<full path to blast> <fragment length int> <alignment lower bound 1..100> "+
        //                       "<evalue lower bound 0..1> <output file - optional>");
        //    System.exit(1);
        //}
        //else
        //{
            RedundancyMatrixGeneration redu = new RedundancyMatrixGeneration(commands);
            //}
            
    }
}
