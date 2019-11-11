package meshi.energy.rapdf;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by N. Vanetik
 */

public class RapdfScore {
    /**
     * This class contains static methods for Samudrala-Moult RAPDF score computation,
     * within an optimization process or outside of it
     **/
    /*** matrix constants ***/    
    public static final int MIN_DISTANCE_BIN = 3; // in angstroms
    public static final int MAX_DISTANCE_BIN = 20; // in angstroms

    public static final int    DEBUG_LEVEL = 1; // debug printouts

    public static final int X = 0;
    public static final int Y = 1;
    public static final int Z = 2;
    public static final int NUMBER_OF_COORDINATES = 3;

    /************************************************************************/
    /*********************  SCORE COMPUTATION METHODS  **********************/
    /************************************************************************/
    /*** compute the score for a pair of getAtoms:
         total score is
         S({a,b}^ij) = - sum_{ij} ln (P_d_ab^ij|c)/P_d_ab^ij - ln P (C | d_ab^ij)
         i and j are the getAtoms of types a and b;

         a single pair score is:
         S(a,b) = - ln (P_d_ab|C)/P_d_ab - P (C | d_ab)  ****/
    public static double score(AtomType a, AtomType b, double distance,
                               double[][][] PDABC, double[] PDAB) {
          // check that neither atom is a hydrogen - was checked before
          // now, compute the score
          double score = 0;
          int place = (int)Math.floor(distance)-MIN_DISTANCE_BIN;

          //System.out.println("score: found distance place "+place);
          if (place<0 || place>=MAX_DISTANCE_BIN-MIN_DISTANCE_BIN+1)
              return 0;

          if (a.compareTo(b) > 0) {
                AtomType tmp = a;
                a = b;
                b = tmp;
            }
          score = (double)PDABC[a.ordinal()][b.ordinal()][place]/(double)PDAB[place];

          //if (DEBUG_LEVEL>=3 && a.ordinal()==88 && b.ordinal()==139)
          //    System.out.println("+++++++++++++++> score["+a.ordinal()+"]["+b.ordinal()+"]["+place+"]="+score);
              
          return score;
        }

        /**** compute the total score for a structure;
          input: array of atom types, atom coordinates and serial residue number
                 for each atom ****/
    public static double totalScore(ArrayList atoms, ArrayList residues,
                                    ArrayList coordinates,
                                    double[][][] PDABC, double[] PDAB)
        {
           double score = 0;
           AtomType a, b;
           Object obj_a, obj_b;
           double distance;
           double rscore = 1;
           double current;
           int residueA, residueB;
           double [] coorA, coorB;
           int atomCount = 0, pairCount = 0;
           Atom atom_a, atom_b;

           coorA = new double[NUMBER_OF_COORDINATES];
           coorB = new double[NUMBER_OF_COORDINATES];

           for (int i=0; i<atoms.size(); i++) {
               obj_a = atoms.get(i);
               if (obj_a instanceof AtomType) /* we have list of atom types */
               {
                  a = (AtomType)obj_a;
                  residueA = ((Integer)residues.get(i)).intValue();
                  coorA = (double[])coordinates.get(i);
               }
               else /*** should be of type Atom ***/
               {
                   atom_a = (Atom)obj_a;
                   if (atom_a.nowhere()) continue;
                   a = atom_a.samudralaType();
                   residueA = atom_a.residue().position();
                   coorA[X] = atom_a.x();
                   coorA[Y] = atom_a.y();
                   coorA[X] = atom_a.z();
               }
               
               if (a.isHydrogen()) continue;
               atomCount++;
               
               for (int j=i+1; j<atoms.size(); j++) 
                 {                    
                    obj_b = atoms.get(j);
                    if (obj_b instanceof AtomType) /* we have list of atom types */
                    {
                       b = (AtomType)obj_b;
                       residueB = ((Integer)residues.get(j)).intValue();
                       coorB = (double[])coordinates.get(j);
                    }
                    else /*** should be of type Atom ***/
                    {
                        atom_b = (Atom)obj_b;
                        if (atom_b.nowhere()) continue;
                        b = atom_b.samudralaType();
                        residueB = atom_b.residue().position();
                        coorB[X] = atom_b.x();
                        coorB[Y] = atom_b.y();
                        coorB[X] = atom_b.z();
                    }
                    
                    if (b.isHydrogen()) continue;
                    
                    if (residueA!=residueB)
                    {
                      pairCount++;  
                      distance = Math.sqrt((coorA[X]-coorB[X])*(coorA[X]-coorB[X])
                                           +(coorA[Y]-coorB[Y])*(coorA[Y]-coorB[Y])+
                                            (coorA[Z]-coorB[Z])*(coorA[Z]-coorB[Z]));

                      current = score(a, b, distance, PDABC, PDAB);
                      int place = (int)Math.floor(distance)-MIN_DISTANCE_BIN;

                      if (DEBUG_LEVEL>=3 && i==0 && j<20 && current!=0)
                        System.out.println("+++++++++++++++> distance ("+i+","+j+")="
                                           +distance+" current="+current
                                           +" types ("+a+","+b+")"
                                           +" place="+place+"\n                 score "+score
                                           +" => "+(score-Math.log(current)));
                      
                      if (current!=0)
                      {
                          //System.out.println(" score="+current+" log score = "+(-1)*Math.log1p(current-1)); 
                         score -= Math.log(current);
                      }
                         rscore *= current;
                         //System.out.println(" current score="+score+" current rscore="+rscore);                      
                    }
                 }
            }

           // now, compute P(C|d_ab) for given node pairs
           if (rscore == 0) rscore = 1;
           score = score - Math.log(rscore);

           if (DEBUG_LEVEL>=3)
           {
               //System.out.println("++++++> rscore="+(-1)*Math.log(rscore)+" final score = "+score);
               System.out.println("++++++> non-hydrogen atom count ="+atomCount);
               System.out.println("++++++> atom count ="+atoms.size());
               System.out.println("++++++> pair count ="+pairCount);
           }
           
           return score;
        }

    
    /***************************************************************************************/
    /*****                  COMPUTING AUXILIARY MATRICES                              *****/
    /***************************************************************************************/
        /*** Compute conditional probability of observing atom types a and b at a distance d:
         P(d_ab|C)=N(d_ab)/sum_{d}N(d_ab)
         N(d_ab) -- the number of data of getAtoms types a and b at distance d
         sum_{d}N(d_ab) -- total number of data of getAtoms types a and b at any distance ***/
    public static void computePDABC(double[][][] matrix, double[][][] PDABC)
        {
           double distanceCount;
            
           for (int i=0; i<matrix.length; i++)
                for (int j=i; j<matrix[i].length; j++)
                {
                    distanceCount = 0; // do not divide by zero 
                    for (int k=0; k<matrix[i][j].length; k++)
                        distanceCount += matrix[i][j][k];
                    
                    /** P(d_ab|C)=N(d_ab)/sum_{d}N(d_ab) **/
                    for (int k=0; k<PDABC[i][j].length; k++)
                        PDABC[i][j][k] = matrix[i][j][k]/distanceCount; 
                }
        }

    /*** Compute unconditional probability of observing atom types a and b at a distance d:
         P(d_ab)=sum_{ab}N(d_ab)/sum_{d}sum_{ab}N(d_ab)
         sum_{ab}N(d_ab) - total number of data of getAtoms of any type at distance d
         sum_{d}sum_{ab}N(d_ab) - total number of data (any atom types, any distance)
    ***/
    public static double computePDAB(double[] PDAB, double[] NDAB_D, double[][][] matrix)
        {
           Arrays.fill(NDAB_D, 0);
           double NDAB_TOTAL = 0;
            
           for (int i=0; i<matrix.length; i++)
              for (int j=i; j<matrix[i].length; j++)
                  for (int k=0; k<matrix[i][j].length; k++)
                  {
                      NDAB_D[k] += matrix[i][j][k];
                      NDAB_TOTAL += matrix[i][j][k];
                  }

           for (int k=0; k<PDAB.length; k++)
              PDAB[k] = NDAB_D[k]/NDAB_TOTAL;

           return NDAB_TOTAL;
        }


    /***************************************************************************************/
    /*****                               PRINTING METHODS                              *****/
    /***************************************************************************************/
    public static void printPDABC(double[][][] PDABC)
        {
           System.out.println("====> PDABC");
           for (int i=0; i<PDABC.length; i++) 
               for (int j=i; j<PDABC[i].length; j++)
               {
                   for (int k=0; k<PDABC[i][j].length; k++)
                       System.out.println(PDABC[i][j][k]);
               }
        }

    /*** Print a single entry in PDABC ***/
    public static void printPDABC(double[][][] matrix, double[][][] PDABC, int i, int j, int place)
        {
           System.out.println("====> PDABC["+i+"]["+j+"]["+place+"]="+PDABC[i][j][place]);
           int distanceCount = 0;
           for (int k=0; k<matrix[i][j].length; k++)
               distanceCount += matrix[i][j][k];
           System.out.println("====> matrix["+i+"]["+j+"]["+place+"]="+matrix[i][j][place]+
               " distance["+i+"]["+j+"]="+distanceCount);
           
        }

    public static void printNDAB_TOTAL(double NDAB_TOTAL)
        {
             System.out.println("====> NDAB_TOTAL = "+NDAB_TOTAL);
        }

    public static void printPDAB(double[] PDAB)
        {
           System.out.println("====> PDAB of length "+PDAB.length);
           
           for (int k=0; k<MAX_DISTANCE_BIN-MIN_DISTANCE_BIN+1; k++)
              System.out.print(PDAB[k]+" ");
           System.out.println();
        }

    public static void printPDAB(double[] PDAB, double[] NDAB_D, double NDAB_TOTAL, int place)
        {
           System.out.println("====> PDAB["+place+"]="+PDAB[place]);
           System.out.println("====> NDAB_D["+place+"]="+NDAB_D[place]+" NDAB_TOTAL="+NDAB_TOTAL);
        }

     public static void printMatrix(double[][][] matrix)
        {
            int count;
            
            System.out.println("====> distance bin matrix: ");
            
            for (int i=0; i<matrix.length; i++)
            {
                for (int j=i; j<matrix[i].length; j++)
                {
                    System.out.print("m["+i+"]["+j+"]=");
                    for (int k=0; k<matrix[i][j].length; k++)
                        System.out.print(" "+matrix[i][j][k]);
                    System.out.println();
                }
            } 
        }

    public static void printMatrix(double[][][] matrix, int i, int j)
        {
            int count;
            
            System.out.println("====> distance bin matrix ["+i+"]["+j+"]");            
            if (i>=0 && j>=i && i<matrix.length && j<matrix[i].length)
               {
                  System.out.print("m["+i+"]["+j+"]=");
                  for (int k=0; k<matrix[i][j].length; k++)
                        System.out.print(" "+matrix[i][j][k]);
                  System.out.println();
               }
        }

    public static void printNDAB_D(double[] NDAB_D)
        {
            System.out.println("====> NDAB_D: ");
            for (int k=0; k<NDAB_D.length; k++)
                System.out.print(NDAB_D[k]+" ");
            System.out.println(); 
        }

}
