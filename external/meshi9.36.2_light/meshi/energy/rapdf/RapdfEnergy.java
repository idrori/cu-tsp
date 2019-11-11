package meshi.energy.rapdf;
import meshi.energy.simpleEnergyTerms.*;

/**
 * Created by N. Vanetik
 */

import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.energy.*;
import meshi.util.CommandList;
import meshi.util.Utils;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;
import meshi.util.mysql.*;
import meshi.parameters.*;

import java.util.*;
import java.io.*;

public class RapdfEnergy extends SimpleEnergyTerm{
    /**
     * This constructor initializes RAPDF energy matrix and
     * computes the RAPDF score of Samudrala-Moult for the given structure
     **/
    
    /*** matrix data variables ***/
    private double[][][] rapdfMatrix = null;
    private double initialLocationsMatrix[][] = null;
    public static final int    DEBUG_LEVEL = 1; // debug printouts

    public static final int X = 0;
    public static final int Y = 1;
    public static final int Z = 2;

    private double PDABC[][][] = null;
    private double PDAB[] = null;
    private int numberOfAtomTypes = -1;
    private double NDAB_TOTAL = 0;
    private double[] NDAB_D = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];

    public static final int RAPDF_TABLE_ATTRIBUTE_NUMBER = 4;

    private AtomList atomList = null;
    private String tableName, rapdfFileNameKey;

    /***
        Constructors
     ***/
    public RapdfEnergy() {}
    
    public RapdfEnergy(AtomList atomList, double[][] inputMatrix, EnergyInfoElement info, String comment,String rapdfFileNameKey) {
        super(toArray(),null, info,EnergyType.NON_DIFFERENTIAL);
        initialLocationsMatrix = inputMatrix;
        this.atomList = atomList;

        elementsList = new ArrayList();
        this.weight = weight;
        for (int i = 0; i < atomList.size(); i++) {
            //for (Iterator baseElements = atomList.iterator(); baseElements.hasNext();) {
            Atom atom = atomList.get(i);
            if (! atom.nowhere())  {
                EnergyElement newElement = createElement(atom );
                if (! newElement.frozen()) {
                    elementsList.add(newElement);
                }
            }
        }
        this.comment = comment;
        this.tableName = "rapdfoutputtable";
        this.rapdfFileNameKey = rapdfFileNameKey;
    }

    /***
        Compute Samudrala-Moult score of this protein
    ***/
    public EnergyInfoElement evaluate() {

     double energy = RapdfScore.totalScore(this.atomList, null, null, PDABC, PDAB);
     if (DEBUG_LEVEL>=3)
         System.out.println("+++++++++++++++> Samudrala-Moult score reset "+this.atomList.size()+" getAtoms is "+energy);
     
     if (info == null) throw new RuntimeException("null info in "+this);
     info.setValue(energy);   
	 return info;
    }


    public EnergyElement createElement(Object baseElement) {
        Atom atom = (Atom) baseElement;
        EnergyElement out = new RapdfEnergyElement(atom, initialLocationsMatrix[0][atom.number()],
                initialLocationsMatrix[1][atom.number()] ,initialLocationsMatrix[2][atom.number()], weight);
        return out;
    }

       public boolean isRapdfEnergy() {
        return true;
      }

   /**
     * Testing of one atom in all energy elements
     *
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom        an criminal <code>Atom</code> value
     * At this point, do nothing!
     * Derivability is not important here !!!
     */
    public void test(TotalEnergy totalEnergy, Atom atom) {
        System.out.println("----------------------------------------> IN RAPDF TEST!!!");
        if (!on) System.out.println("" + this + " is off");
        return;
    }
    
    /***
     Read input parameters, determine where to get the RAPDF matrix from and store it in memory.
     Compute auxiliary matrices from the one read from DB or a file.
     ***/
    public void set(CommandList commands) {
        
        if (commands.keyExists(tableName)) {
                String rapdfTable = commands.firstWord(tableName).secondWord();
                /*** establish mysql conection ***/
                MysqlHandler mysql_handler = new MysqlHandler(commands);
                
                /*** read rapdf data into memory
                     and use it for further computation ***/
                readRapdfMatrixFromDb(mysql_handler, rapdfTable);
                if (rapdfMatrix != null)
                    computeAuxiliaryMatrices();
                
                /*** close mysql connection ***/
                mysql_handler.close();

        }

        if (rapdfMatrix==null && commands.keyExists(rapdfFileNameKey))  {
            /*** open file ***/
            String rapdfFileName = commands.firstWord(rapdfFileNameKey).secondWord();
            
            /*** read rapdf dat into memory
             and use it for further computation ***/
            readRapdfMatrixFromFile(rapdfFileName);
            computeAuxiliaryMatrices();
        }

        if (rapdfMatrix==null)  {
            throw new  RuntimeException("No RAPDF data table or file specified in input, cannot compute Samudrala score");
        }
    }


    /*** compute additional data based reset matrix ***/
    private void computeAuxiliaryMatrices()
        {
           PDABC = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
           PDAB = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
           NDAB_TOTAL = RapdfScore.computePDAB(PDAB, NDAB_D, rapdfMatrix);
           RapdfScore.computePDABC(rapdfMatrix, PDABC);

           if (DEBUG_LEVEL >=3)
                 {
                   RapdfScore.printPDAB(PDAB);
                   RapdfScore.printPDABC(PDABC);
                 }
        }

    /********************************************************************************************/
    /***********************             ENERGY METHODS  ***************************************/
    /********************************************************************************************/
    public void update(){}
    public void update(int i ){}

    public EnergyElement createElement(Object baseElement,Parameters p){return null;}

    public void reset() {
        for (EnergyElement e : elementsList) {
            RapdfEnergyElement element = (RapdfEnergyElement) e;
            element.reset();
        }
    }

    public void remove(Atom atom) {
        for (EnergyElement e : elementsList) {
            RapdfEnergyElement element = (RapdfEnergyElement) e;
            if (element.atom == atom) {
                element.remove();
            }
        }
    }

    /************************************************************************/
    /*********************    MATRIX READING METHODS   **********************/
    /************************************************************************/
    private void readRapdfMatrixFromFile(String filename)
      {
          if (DEBUG_LEVEL>=1)
              Utils.println("+++++++++++++++> Reading rapdf matrix from file="+filename);
          try
          {
             LineNumberReader lnr = new LineNumberReader(new FileReader(filename));
             String line;
             int size, a, b;
             StringTokenizer st;
             AtomType type1, type2;
             
             this.numberOfAtomTypes = AtomType.values().length;
             // init the matrix
             rapdfMatrix = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
             
             while ((line = lnr.readLine()) != null)
             {
                 st = new StringTokenizer(line);
                 size = st.countTokens();
                 if (size>=RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1+2)
                 {
                   /*** First two tokens are atom types ***/
                     type1 = AtomType.type(st.nextToken());
                     type2 = AtomType.type(st.nextToken());
                     a = type1.ordinal();
                     b = type2.ordinal();
                   /*** following tokens are distance histogram ***/
                     Arrays.fill(rapdfMatrix[a][b], 1);

                     for (int i=RapdfScore.MIN_DISTANCE_BIN; i<=RapdfScore.MAX_DISTANCE_BIN; i++)
                     {
                         try {
                             rapdfMatrix[a][b][i-RapdfScore.MIN_DISTANCE_BIN] = Double.parseDouble(st.nextToken()); }
                         catch (NumberFormatException e1) {}
                     }
                 }
             }
          }
          catch (IOException e)
          {
              e.printStackTrace();
              on = false;
          }
      }

    private void readRapdfMatrixFromDb(MysqlHandler mysql_handler, String tablename)
      {
             if (DEBUG_LEVEL>=1)
              Utils.println("+++++++++++++++> Reading rapdf matrix from db table="+tablename);
             
             int size, a, b, k;
             AtomType type1, type2;
             double value;

             if (!mysql_handler.connected() || !mysql_handler.tableExists(tablename))
             {
                 Utils.println("+++++++++++++++> Db table="+tablename+" was not found");
                 rapdfMatrix = null;
                 return;
             }
             
             this.numberOfAtomTypes = AtomType.values().length;
             // init the matrix
             rapdfMatrix = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];

             ArrayList tableContents = mysql_handler.selectAll(tablename, RAPDF_TABLE_ATTRIBUTE_NUMBER);
             String row[] = null;

             for (int i=0; i<numberOfAtomTypes; i++)
                 for (int j=0; j<numberOfAtomTypes; j++)
                     Arrays.fill(rapdfMatrix[i][j], 1);
             
             for (int i=0; i<tableContents.size(); i++)
             {
                 row = (String[])tableContents.get(i);
                 /*** First two tokens are atom types ***/
                 type1 = AtomType.type(row[0]);
                 type2 = AtomType.type(row[1]);
                 a = type1.ordinal();
                 b = type2.ordinal();

                 /*** following tokens are distance histogram ***/
                 try {
                         k = Integer.parseInt(row[2]);
                         value = Double.parseDouble(row[3]);
                         rapdfMatrix[a][b][k] = value;
                    }
                         catch (NumberFormatException e1) {}                 
             }
          }

    public boolean evaluatesResidues() {
        return false;
    }

}


