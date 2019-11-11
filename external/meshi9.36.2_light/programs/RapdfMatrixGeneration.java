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

import java.io.*;
import java.util.*;

/***
    Data: 26/12/2010
    Creator: N. Vanetik
    This programs reads files in PDB format with known structure and builds the Samudrala-Moult matrix
    for the RAPDF score computation. This matrix is stored in a file.
    The program receives following parameters:
        <directory that holds PDB files> 
***/
public class RapdfMatrixGeneration {

    public static final int    DEBUG_LEVEL = 1; // debug printouts
    public static final String PDB_ATOM_RECORD = "ATOM";
    public static final String PDB_EXTENSION = ".pdb";
    public static final double DEFAULT_REDUNDANCY_WEIGHT = 1.0;
    public static final int    PRINT_LIMIT = 20; // limit of data items per line for printing matrices

    public static final int    READ_REDUNDANCY_FROM_NOWHERE   = -1;
    public static final int    READ_REDUNDANCY_FROM_FILE      = 0;
    public static final int    READ_REDUNDANCY_FROM_DATABASE  = 1;

    // username, password and DB name for mysql access
    private static final String redundancyAttribute = "redu";
    private static final String residueAttribute = "residue";

    private int whereToGetRedundancyFrom = READ_REDUNDANCY_FROM_NOWHERE;
    
    /*** matrix data variables ***/
    private double[][][] matrix = null;
    private double PDABC[][][] = null;
    private double PDAB[] = null;
    private int numberOfAtomTypes = -1;
    private double NDAB_TOTAL = 0;
    private double[] NDAB_D = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];

    /*** redundancy data ***/
    private ArrayList redundancy = null;

    /*** Access to data storage ***/
    private CommandList  commands = null;
    private String       queryDirectory = null;
    private String       outputFile = null;
    private String       outputTable = null;
    private String       redundancyFile = null;
    private String       redundancyTable = null;
    private String       queryPath = null;

    MysqlHandler mysql_updater = null; // connection to a mysql database

    int yce_pair_count = 0; // for debugging

    /*** constructors & initializers ***/
    public RapdfMatrixGeneration(CommandList commands)
        {
            this.commands = commands;

            this.queryDirectory = commands.firstWord("rapdffilesrepository").secondWord();
            System.out.println("Found query directory="+this.queryDirectory);
            this.outputFile = commands.firstWord("rapdfoutputfile").secondWord();
            System.out.println("Found output file="+this.outputFile);
            this.outputTable = commands.firstWord("rapdfoutputtable").secondWord();
            System.out.println("Found output table="+this.outputTable);
            this.redundancyFile = commands.firstWord("rapdfredundancyfile").secondWord();
            System.out.println("Found redundancyFile ="+this.redundancyFile);
            this.redundancyTable = commands.firstWord("rapdfredundancytable").secondWord();
            System.out.println("Found redundancyTable ="+this.redundancyTable);

            if (commands.keyExists("rapdfquerypath"))
            {
                this.queryPath = commands.firstWord("rapdfquerypath").secondWord();
                System.out.println("Found query path ="+this.queryPath);
            }
            
            this.mysql_updater = new MysqlHandler(commands);
            
            whereToGetRedundancyFrom = findRedundancy();
            System.out.println("whereToGetRedundancyFrom="+whereToGetRedundancyFrom);
            
            initData(this.queryDirectory);
            compute();
            storeMatrix();
            mysql_updater.close();
            computeScores();
        }

    /***
        Decide where to get redundancy from.
        If a database connection exists, and the appropriate table exists, then redundancy should be
        read from the database - it is time-optimal.
        Otherwise, check that the file exists and is of appropriate format. If so, read redundancy from file.
        If no db table and no file was found, work without redundancy.
     ***/
    private int findRedundancy()
        {
            if (this.redundancyTable != null && mysql_updater!=null &&
                mysql_updater.tableExists(this.redundancyTable))
            {
                System.out.println("Reading redundancy from database table "+this.redundancyTable);
                return READ_REDUNDANCY_FROM_DATABASE; 
            }
            this.redundancy = new ArrayList();
            if (this.redundancyFile != null && loadRedundancy(this.redundancy, this.redundancyFile))
            {
                System.out.println("Reading redundancy from file "+this.redundancyFile);
                return READ_REDUNDANCY_FROM_FILE;
            }
            return READ_REDUNDANCY_FROM_NOWHERE;
        }

    /***
        Initialization: allocate space for the matrix and arrays to speed-up the computation
     ***/
    private void initData(String queryDirectory)
        {
            this.queryDirectory = queryDirectory;
            if (queryDirectory == null)
            {
                System.out.println("Directory for ramdf computation is not defined, stopping...");
                System.exit(0);
            }
            
            this.numberOfAtomTypes = AtomType.values().length;
            // init the matrix
            matrix = new double[numberOfAtomTypes][numberOfAtomTypes][];
            // init distance bins, and make the matrix diagonal
            for (int i=0; i<numberOfAtomTypes; i++)
                for (int j=i; j<numberOfAtomTypes; j++)
                {
                    matrix[i][j] = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
                    Arrays.fill(matrix[i][j], 1);
                }
        }


    /************************************************************************/
    /*********************   GETTERS & SETTERS         **********************/
    /************************************************************************/
    public double[][][] getMatrix() { return matrix; }

    /************************************************************************/
    /*********************  MATRIX COMPUTATION METHODS **********************/
    /************************************************************************/
    public void compute()
        {
            PDABC = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            PDAB = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            
            fillMatrix(queryDirectory);
            NDAB_TOTAL = RapdfScore.computePDAB(PDAB, NDAB_D, matrix);
            RapdfScore.computePDABC(matrix, PDABC);

            if (DEBUG_LEVEL >=3)
            {
               RapdfScore.printPDAB(PDAB);
               RapdfScore.printPDABC(PDABC);
            }
        }

    /*** Update the entire matrix:
         (1) go over all *.pdb files in a given directory
         (2) in each file, go over all ATOM records
         (2.5) determine atom type with the help of AtomType enum
         (3) store atom coordinates in local array
         (4) keep track of residues
         (5) for each pair of getAtoms NOT IN THE SAME RESIDUE, toggle the counters by calling updateDistance ***/
    public void fillMatrix(String pdbDirectory)
        {
            File[] files;
            LineNumberReader pdbFileReader;
            String previousResidue = "", line;
            int residueNumber = 0;
            String token, atomToken;
            StringTokenizer st;
            ArrayList atoms = new ArrayList();
            ArrayList coordinatesList = new ArrayList();
            ArrayList residues = new ArrayList();
            ArrayList Xcoordinate = new ArrayList();
            double distance;
            double[] coor1, coor2;
            AtomType atom1, atom2;
            int place;

            if (DEBUG_LEVEL>=1)
              System.out.println("=========> Reading directory "+pdbDirectory);
            
            try {
                 // open the Query directory - it should contain .pdb files
                 File path = new File(pdbDirectory);
                 files = path.listFiles();
                 Arrays.sort(files);
                 for (int i=0; i<files.length; i++)
                    {
                        
                        if (!files[i].getName().endsWith(PDB_EXTENSION))
                            continue;
                        if (DEBUG_LEVEL>=1)
                          System.out.println("=========> Reading PDB file "+files[i]);
                        
                        pdbFileReader = new LineNumberReader(new FileReader(files[i]));
                        while ((line = pdbFileReader.readLine()) != null)
                        {
                            // if this is an atom record
                            if (line.trim().startsWith(PDB_ATOM_RECORD))
                            {
                                //System.out.println("=========>  found ATOM record "+line);
                                double[] coordinates = new double[RapdfScore.NUMBER_OF_COORDINATES];
                                // read the ATOM record
                                
                                st = new StringTokenizer(line);

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                // 3rd element is the atom name
                                if (st.hasMoreTokens())
                                    atomToken = st.nextToken();
                                else continue;
                                
                                
                                if (st.hasMoreTokens())
                                    token = st.nextToken();
                                else continue;
                                
                                // 4th element is the residue name
                                if (!previousResidue.equals(token))
                                {
                                    residueNumber++;
                                    previousResidue = token; 
                                }

                                if (token.length()<3 || token.length()>3)
                                    continue;

                                //System.out.println("=========> line="+line);
                                //System.out.println("=========>  found ATOM of type "+atomToken+" residue "+token);
                                AtomType atom = AtomType.type(token, atomToken);
                                
                                /*** exclude hydrogen - Samudrala-Moult optimization ***/
                                if (atom.isHydrogen())
                                    continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                // 8th,9th and 10th elements are X,Y and Z coordinates
                                try {
                                    coordinates[RapdfScore.X] = Double.parseDouble(st.nextToken());
                                    coordinates[RapdfScore.Y] = Double.parseDouble(st.nextToken());
                                    coordinates[RapdfScore.Z] = Double.parseDouble(st.nextToken());

                                    //System.out.println("read atom "+atom+" of residue "+
                                    //                   previousResidue+" with coordinates x="+coordinates[X]
                                    //                   +" y="+coordinates[Y]+" z="+coordinates[Z]);
                                    Double Xcoor = new Double(coordinates[RapdfScore.X]);
                                    //place = - insertion point - 1
                                    
                                    place = Arrays.binarySearch(Xcoordinate.toArray(), Xcoor);
                                    //System.out.println("place = "+place);
                                    // now - insertion point
                                    if (place<0)
                                        place = (place +1 )*(-1);
                                    else
                                        place = place+1;
                                    //System.out.println("now place = "+place);
                                    
                                    // insert sorted
                                    atoms.add(place, atom);
                                    residues.add(place, new Integer(residueNumber));
                                    coordinatesList.add(place, coordinates);
                                    Xcoordinate.add(place, Xcoor);

                                    // for a matrix entry                                
                                    /*getAtoms.add(atom);
                                    residues.add(new Integer(residueNumber));
                                    coordinatesList.add(coordinates);*/
                                }
                                catch (NumberFormatException e2)
                                {
                                }
                                
                            } // if
                        } // while

                        pdbFileReader.close();

                        ArrayList currentRedundancy = null;
                        if (whereToGetRedundancyFrom == READ_REDUNDANCY_FROM_DATABASE)
                            currentRedundancy = getAllWeightsFromDatabase(files[i].getName());


                        for (int j=0; j<atoms.size(); j++)
                        {
                            coor1 = (double[])coordinatesList.get(j);
                            atom1 = (AtomType)atoms.get(j);
                            
                             for (int k=j+1; k<atoms.size(); k++)    
                             {
                                 atom2 = (AtomType)atoms.get(k);
                                 coor2 = (double[])coordinatesList.get(k);

                                 
                                 if (Math.abs(coor1[RapdfScore.X]-coor2[RapdfScore.X])>RapdfScore.MAX_DISTANCE_BIN+1)
                                 {
                                     k = atoms.size();
                                     continue;
                                 }

                                 // update matrix entries
                                 if ((((Integer)residues.get(j)).intValue() != ((Integer)residues.get(k)).intValue())
                                     && Math.abs(coor1[RapdfScore.Y]-coor2[RapdfScore.Y])<=RapdfScore.MAX_DISTANCE_BIN+1
                                     && Math.abs(coor1[RapdfScore.Z]-coor2[RapdfScore.Z])<=RapdfScore.MAX_DISTANCE_BIN+1)
                                 {
                                     //coor2 = (double[])coordinatesList.get(k);
                                    distance = Math.sqrt((coor1[RapdfScore.X]-coor2[RapdfScore.X])*(coor1[RapdfScore.X]-coor2[RapdfScore.X])+
                                        (coor1[RapdfScore.Y]-coor2[RapdfScore.Y])*(coor1[RapdfScore.Y]-coor2[RapdfScore.Y])+
                                        (coor1[RapdfScore.Z]-coor2[RapdfScore.Z])*(coor1[RapdfScore.Z]-coor2[RapdfScore.Z]));

                                    updateDistance(atom1, atom2, distance,
                                                   i, currentRedundancy,
                                                   ((Integer)residues.get(j)).intValue(),
                                                   ((Integer)residues.get(k)).intValue());
                                 }
                             }
                        }

                        // clearRows lists
                        Xcoordinate.clear();
                        atoms.clear();
                        residues.clear();
                        coordinatesList.clear();
                        if (currentRedundancy != null)
                            currentRedundancy.clear();

                        if (DEBUG_LEVEL>=3)
                            System.out.println("yce_pair_count="+yce_pair_count);
                        
                    } //for
            } // try
           catch (IOException e)
            {
                System.out.println("Exception "+e.toString()
                                   +" while reading file "+pdbDirectory);
            }
            
        }


     /*** Update the matrix pointwise, given a pair of getAtoms and a distance ***/
    public void updateDistance(AtomType a, AtomType b, double distance, int sequence, ArrayList dbRedundancy, int res1, int res2)
        {
          if (a.compareTo(b) > 0)
            {
                AtomType tmp = a;
                a = b;
                b = tmp;
            }
          // check that neither atom is a hydrogen
          if (a.isHydrogen() || b.isHydrogen())
              return;

          //if (a == AtomType.XXX && b == AtomType.XXX)
          //   return;

          // only count distances between min distance and max distance
          int place = (int)Math.floor(distance);
          double w;
          
          if (distance>=RapdfScore.MIN_DISTANCE_BIN && distance<RapdfScore.MAX_DISTANCE_BIN+1)
            {
                if (DEBUG_LEVEL>=3 && a == AtomType.YCE && b == AtomType.YCE)
                   yce_pair_count++;
                
                // if redundancy data is present, update the count according to redundancy
                w = getRedundancyWeight(res1, res2, sequence, dbRedundancy);

                if (DEBUG_LEVEL>=3 && w!=1)
                    System.out.println("-----------> w="+w);
                
                // all counts are < 1 right now
                matrix[a.ordinal()][b.ordinal()][place-RapdfScore.MIN_DISTANCE_BIN] += w;
            }
        }

    /***
        Compute redundancy measure (weight) for two specific pair of residues in a specific sequence (file)
     ***/
    private double getRedundancyWeight(int res1, int res2, int sequence, ArrayList dbRedundancy)
        {
            switch (whereToGetRedundancyFrom)
            {
            case READ_REDUNDANCY_FROM_FILE:
                return getRedundancyWeightFromMemory(res1, res2, sequence);
            case READ_REDUNDANCY_FROM_DATABASE:
                return getRedundancyWeightFromDatabase(res1, res2, dbRedundancy);
            }
            return DEFAULT_REDUNDANCY_WEIGHT;
        }

    /***
        Compute redundancy measure (weight) for two specific pair of residues in a specific sequence (file)
     ***/
    private double getRedundancyWeightFromMemory(int res1, int res2, int sequence)
        {
            if (redundancy==null || sequence<0 || sequence>=redundancy.size())
                return DEFAULT_REDUNDANCY_WEIGHT;
            
            double w[] = (double[])redundancy.get(sequence);
            if (res1>=w.length || res2>=w.length)
                return DEFAULT_REDUNDANCY_WEIGHT;

            if (DEBUG_LEVEL>=2)
                System.out.println("++++++> redundancy weight ="+w[res1]*w[res2]);

            return w[res1]*w[res2];
            //return (w[res1]+w[res2])/2;
        }

    /***
        Compute redundancy measure (weight) for two specific pair of residues in a specific sequence (file)
     ***/
    private double getRedundancyWeightFromDatabase(int res1, int res2, ArrayList dbRedundancy)
        {
            //if (DEBUG_LEVEL>=1)
            //   System.out.println("++++++> dbRedundancy.size()="+dbRedundancy.size()+" res1="+res1+" res2="+res2);
            
            if (dbRedundancy==null || res1>=dbRedundancy.size() || res2>=dbRedundancy.size())
                return DEFAULT_REDUNDANCY_WEIGHT;

            double w1 = ((Double)dbRedundancy.get(res1)).doubleValue();
            double w2 = ((Double)dbRedundancy.get(res2)).doubleValue();

            if (DEBUG_LEVEL>=3)
               System.out.println("++++++> redundancy weight ="+w1*w2);
                
            if (w1>0 && w2>0)
            {
                return w1*w2;  
            }            
            return DEFAULT_REDUNDANCY_WEIGHT;  
        }

    /***
       Compute redundancy measure (weight) for all residues in a specific file
    ***/
    private ArrayList getAllWeightsFromDatabase(String fileName)
        {
            String whereCondition = "filename='"+fileName+"' ORDER BY "+residueAttribute;
            String selectCondition = redundancyAttribute+","+residueAttribute;
            
            return  mysql_updater.selectDoubleList(this.redundancyTable,
                                                         whereCondition,
                                                         selectCondition,
                                                         redundancyAttribute);
        }


    /***************************************************************************************/
    /*****                               LOAD REDUNDANCY                               *****/
    /***************************************************************************************/
    public static boolean loadRedundancy(ArrayList redundancy, String filename)
        {
            try
            {
               // open file
                if (!(new File(filename)).exists())
                    return false;

                // open existing file
               LineNumberReader lnr = new LineNumberReader(new FileReader(filename));
               String line;
               StringTokenizer st;
               int num;

               
               if (DEBUG_LEVEL >=1)
                   System.out.println("====> In loadRedundancy");
                   
               // read line by line
               while ((line = lnr.readLine()) != null)
               {                   
                   // every line is redundancy of file #<current line>
                   st = new StringTokenizer(line);
                   // skip the first token - it is a line number with >
                   if (st.hasMoreTokens())
                   {
                       st.nextToken();
                       num = st.countTokens();
                       double W[] = new double[num];
                       try
                       {
                         for (int i=0; i<num; i++)
                             W[i] = Double.parseDouble(st.nextToken());
                         redundancy.add(W);
                       }
                      catch (NumberFormatException e1)
                      {
                      }
                   }
               }
               if (DEBUG_LEVEL >=1)
                   System.out.println("====> In loadRedundancy found "+redundancy.size()+" lines");
               return true;
               
            }
            catch (IOException e)
            {
                e.printStackTrace();
                redundancy = null;
                return false;
            }
        }


    /***************************************************************************************/
    /*****                               STORAGE METHODS                               *****/
    /***************************************************************************************/
    public void storeMatrix()
        {
            if (this.outputFile != null)
                storeMatrixInFile(this.outputFile);
            if (this.outputTable != null)
                storeMatrixInTable(this.outputTable);
        }

    /***
        store in database 
     ***/
    public void storeMatrixInTable(String outTable)
        {
            
            boolean rc = mysql_updater.createTable(outTable, "(i varchar(10), j varchar(10), k integer, value double)");
            String parameters = "(i, j, k, value)";
            String values;
            if (rc)
            {
              for (int i=0; i<numberOfAtomTypes; i++)
               {
                for (int j=i; j<numberOfAtomTypes; j++)
                {
                    for (int k=0; k<matrix[i][j].length; k++)
                    {
                        values = "('"+AtomType.type(i)+"','"+AtomType.type(j)+"',"+k+","+matrix[i][j][k]+")";
                        mysql_updater.insertIntoTable(outTable, parameters, values);
                    }
                }
              }
            }              
        }

    /***
        store in file 
     ***/
    public void storeMatrixInFile(String outFile)
        {
            StringBuffer sb;

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 
              for (int i=0; i<numberOfAtomTypes; i++)
               {
                for (int j=i; j<numberOfAtomTypes; j++)
                {
                    sb = new StringBuffer();
                    sb.append(AtomType.type(i)+" "+AtomType.type(j));
                    for (int k=0; k<matrix[i][j].length; k++)
                        sb.append(" "+matrix[i][j][k]);
                    sb.append("\n");
                    fos.write(sb.toString().getBytes());
                }
              }
              fos.close();
            }
            catch (IOException e)
            {
                System.out.println("Error while storing RAPDF matrix in file "+outFile);
                e.printStackTrace();
            }
        }

    /***************************************************************************************/
    /*****             COMPUTE SCORE FROM THE PREPARED MATRIX                          *****/
    /***************************************************************************************/
    public void computeScores()
        {
            if (this.queryPath==null)
                return;

            File qf = new File(this.queryPath);
            if (!qf.exists())
                return;

            computeSingleScore(this.queryPath);
        }

    /*** Read PDB structure from a file and compute its score ***/
    public void computeSingleScore(String pdbFileName)
        {
           LineNumberReader pdbFileReader;
           File  pdbFile, resultsFile;
           File  files[];
           StringTokenizer st;
           String previousResidue = "", line;
           ArrayList atoms = new ArrayList();
           ArrayList coordinatesList = new ArrayList();
           ArrayList residues = new ArrayList();
           int residueNumber = 0;
           String token, atomToken;
           
           try {
              
              pdbFile = new File(pdbFileName);
              if (pdbFile.isDirectory())
              {
                  files = pdbFile.listFiles(); 
              }
              else
              {
                  files = new File[1];
                  files[0] = pdbFile;
              }

              for (int i=0; i<files.length; i++)
              {
                  if (!files[i].getName().endsWith(PDB_EXTENSION))
                  {
                      if (DEBUG_LEVEL >=2)
                        System.out.println("=========>  File "+files[i].getName()+" is not a .pdb file, skipped");
                      continue;
                  }
                  
               pdbFileReader = new LineNumberReader(new FileReader(files[i]));
               while ((line = pdbFileReader.readLine()) != null)
                {
                  // if this is an atom record
                  if (line.trim().startsWith(PDB_ATOM_RECORD))
                   {
                       //System.out.println("=========>  found ATOM record "+line);
                                double[] coordinates = new double[RapdfScore.NUMBER_OF_COORDINATES];
                                // read the ATOM record
                                
                                st = new StringTokenizer(line);

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                // 3rd element is the atom name
                                if (st.hasMoreTokens())
                                    atomToken = st.nextToken();
                                else continue;
                                
                                
                                if (st.hasMoreTokens())
                                    token = st.nextToken();
                                else continue;
                                
                                // 4th element is the residue name
                                if (!previousResidue.equals(token))
                                {
                                    residueNumber++;
                                    previousResidue = token; 
                                }

                                AtomType atom = AtomType.type(token, atomToken);
                                //System.out.println("=========>  found ATOM of type "+atomToken+" setResidue to type "+atom);

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                if (st.hasMoreTokens())
                                    st.nextToken();
                                else continue;

                                // 8th,9th and 10th elements are X,Y and Z coordinates
                                try {
                                    coordinates[RapdfScore.X] = Double.parseDouble(st.nextToken());
                                    coordinates[RapdfScore.Y] = Double.parseDouble(st.nextToken());
                                    coordinates[RapdfScore.Z] = Double.parseDouble(st.nextToken());

                                    //System.out.println("read atom "+atom+" of residue "+
                                    //                   previousResidue+" with coordinates x="+coordinates[X]
                                    //                   +" y="+coordinates[Y]+" z="+coordinates[Z]);
                                }
                                catch (NumberFormatException e2)
                                {
                                }

                                // for a matrix entry
                                atoms.add(atom);
                                residues.add(new Integer(residueNumber));
                                coordinatesList.add(coordinates);
                   } // if - this is an atom record
               } // while - reading the file

               
               
               double totalScore = RapdfScore.totalScore(atoms, residues, coordinatesList, PDABC, PDAB);
               if (DEBUG_LEVEL >=1)
                 System.out.println("=========> PDB file = "+files[i].getName()+" score = "+totalScore);
               
               pdbFileReader.close();
               // clearRows the structures for the next computation
               atoms.clear();
               residues.clear();
               coordinatesList.clear();
               
              } // for - files in a directory
              
           } // try
            catch (IOException e)
            {
                System.out.println("Exception "+e.toString()
                                   +" while reading file "+pdbFileName);
            }

    }


    /*******************************************************************************************/
    /*************************  THE MAIN PROGRAM TO BE CALLED   ********************************/
    /*******************************************************************************************/
    public static void main(String[] args) throws IOException{
        
        if (args==null || args.length<1)
        {
            System.out.println("FORMAT: RapdfMatrixGeneration <commands file>");
            System.exit(1);
        }
        else
        {
            CommandList commands = Utils.init(args, 1, "RedundancyMatrixGeneration");
            RapdfMatrixGeneration rapdf = new RapdfMatrixGeneration(commands);
        }
            
    }
}
