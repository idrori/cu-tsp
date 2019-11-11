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
import java.text.*;

/***
    Data: 17/02/2011
    Creator: N. Vanetik
    This programs reads two files containing Samudrala-Moult matrices
    and computes several measures reset them.
    Following measures are computed for each pair of getAtoms separately for each matrix:
       - median for distance counts
       - L1 distance of distance counts
       - L2 distance of distance counts 
    The program receives following parameters:
        <first Samudrala-Moult matrix file> <second Samudrala-Moult matrix file> <output file> 
***/
public class RapdfMatrixComparison{

    private static final int DEBUG_LEVEL = 1;
    private double rapdf1[][][]; // matrix #1
    private double rapdf2[][][]; // matrix #2

    private int numberOfAtomTypes = -1;
    
    /*** constructors & initializers ***/
    public RapdfMatrixComparison(String matrix1, String matrix2, String outputFileName)
        {
            this.numberOfAtomTypes = AtomType.values().length;
            rapdf1 = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            rapdf2 = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];

            double pdabc1[][][] = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            double pdab1[] = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            double[] ndab_d1 = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];

            double pdabc2[][][] = new double[numberOfAtomTypes][numberOfAtomTypes][RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            double pdab2[] = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            double[] ndab_d2 = new double[RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1];
            if (!readRapdfMatrixFromFile(rapdf1, matrix1))
		{
                    System.out.println("Error while reading file"+matrix1);
		    System.exit(0);
		}

            if (!readRapdfMatrixFromFile(rapdf2, matrix2))
		{
                    System.out.println("Error while reading file"+matrix2);
		    System.exit(0);
		}

            double ndab_total_1 = computeAuxiliaryMatrices(rapdf1, pdabc1, pdab1, ndab_d1);
            double ndab_total_2 = computeAuxiliaryMatrices(rapdf2, pdabc2, pdab2, ndab_d2);


            if (DEBUG_LEVEL>=1)
               System.out.println("Comparing matrix "+matrix1+" with matrix "+matrix2);
            
            /*if (!computeObservedVsExpected(rapdf1, rapdf2, outputFileName, pdabc1, 
                         pdabc2, pdab1, pdab2, ndab_d1, ndab_d2,
                         ndab_total_1, ndab_total_2))*/
            if (!computeResampling(rapdf1, rapdf2, outputFileName, false))
               {
                    System.out.println("Error while comparing files");
		    System.exit(0);
		}
        }

    /************************************************************************/
    /*********************    MATRIX READING METHOD    **********************/
    /************************************************************************/
    private boolean readRapdfMatrixFromFile(double[][][] rapdfMatrix, String filename)
      {
          if (DEBUG_LEVEL>=1)
              System.out.println("+++++++++++++++> Reading rapdf matrix from file="+filename); 
          try
          {
             LineNumberReader lnr = new LineNumberReader(new FileReader(filename));
             String line;
             int size, a, b;
             StringTokenizer st;
             AtomType type1, type2;
             
             this.numberOfAtomTypes = AtomType.values().length;
             // fill the matrix
             
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
                                rapdfMatrix[a][b][i-RapdfScore.MIN_DISTANCE_BIN] = Double.parseDouble(st.nextToken()); 
                                rapdfMatrix[b][a][i-RapdfScore.MIN_DISTANCE_BIN] = rapdfMatrix[a][b][i-RapdfScore.MIN_DISTANCE_BIN];
                             }
                         catch (NumberFormatException e1) {}
                     }
                 }
             }
             return true;
          }
          catch (IOException e)
          {
              e.printStackTrace();
              return false;
          }
      }

      /*** compute additional data based reset matrix ***/
    private double computeAuxiliaryMatrices(double[][][] rapdfMatrix, double[][][] PDABC, double[] PDAB, double[] NDAB_D)
        {
           double NDAB_TOTAL = RapdfScore.computePDAB(PDAB, NDAB_D, rapdfMatrix);
           RapdfScore.computePDABC(rapdfMatrix, PDABC);

           if (DEBUG_LEVEL >=3)
                 {
                   RapdfScore.printPDAB(PDAB);
                   RapdfScore.printPDABC(PDABC);
                 }
           return NDAB_TOTAL;
        }

      /************************************************************************/
      /*********************        WRITING METHOD       **********************/
     /************************************************************************/
     /***
      Compute distances for each pair of getAtoms: difference of means, L1, L2 and Linfinity
      ***/
    private boolean computeDistances(double[][][] matrix1, double[][][] matrix2, String outFile)
    {
           StringBuffer sb;
           double median1, median2, L1, L2, Linfinity;
           double medianDifference = 0, L1difference = 0, L2difference = 0, Linfinity_difference=0;
           int count1, count2;
	   int total = 0;
           
            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <pair count 1> <pair count 2> <median abs difference> <L1> <L2> <Linfty>\n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               {
                   
		   for (int j=0; j<numberOfAtomTypes; j++) // was j=i
                {
                    total++;
                    // here, compute median
                    // compute L1 distance
                    // compute L2 distance
                    L1 = 0; L2 = 0; median1 = 0; median2 = 0; Linfinity = 0; 
		    count1 = 0; count2 = 0; 
                    double[] arr1 = new double[matrix1[i][j].length];
                    double[] arr2 = new double[matrix2[i][j].length];
                    arr1 = Arrays.copyOf(matrix1[i][j], matrix1[i][j].length);
                    arr2 = Arrays.copyOf(matrix2[i][j], matrix1[i][j].length);
                    Arrays.sort(arr1);
                    Arrays.sort(arr2);

                    if (arr1.length%2==0)
			{ // find the middle pair and compute the average
			    median1 = (arr1[arr1.length/2-1]+arr1[arr1.length/2])/2;
                            median2 = (arr2[arr2.length/2-1]+arr2[arr2.length/2])/2;
			}
		    else
			{
			    // use the middle element
			    median1 = arr1[arr1.length/2];
                            median2 = arr2[arr2.length/2]; 
			}

                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    L1 += Math.abs(matrix1[i][j][k]-matrix2[i][j][k]);
                            L2 += (matrix1[i][j][k]-matrix2[i][j][k])*(matrix1[i][j][k]-matrix2[i][j][k]);
                            count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                            if (Linfinity < Math.abs(matrix1[i][j][k]-matrix2[i][j][k]))
				Linfinity = Math.abs(matrix1[i][j][k]-matrix2[i][j][k]);
                        }

                    L2 = Math.sqrt(L2);
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                           sb = new StringBuffer(); 
                           sb.append(AtomType.type(i)+"\t"+AtomType.type(j)+
                                     "\t"+(median1-median2)
                                     +"\t"+L1+"\t"+L2+"\t"+Linfinity+"\n");
		           fos.write(sb.toString().getBytes());
			}

                    medianDifference += Math.abs(median1-median2);
                    L1difference += L1;
                    L2difference += L2;
		    Linfinity_difference += Linfinity;
                }		    
	       }
              fos.close();
              if (DEBUG_LEVEL>=1)
		  System.out.println("Average difference is median="+(medianDifference/total)
                                     +" L1="+(L1difference/total)+" L2="+(L2difference/total)+
				     " Linfinity="+(Linfinity_difference/total));
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }


    /***
     Compute Kullback-Lieber divergencies for distance distributions of every pair of atom types
     ***/
    public boolean computeKLdivergency(double[][][] matrix1, double[][][] matrix2, String outFile)
    {
           StringBuffer sb;
           int count1, count2;
	   double KL1, KL2;

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <KL 1->2 divergency> <KL 2->1 divergency>\n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {
                    
                    count1 = 0; count2 = 0;
                    double[] arr1 = new double[matrix1[i][j].length];
                    double[] arr2 = new double[matrix2[i][j].length];
                    arr1 = Arrays.copyOf(matrix1[i][j], matrix1[i][j].length);
                    arr2 = Arrays.copyOf(matrix2[i][j], matrix1[i][j].length);
                    Arrays.sort(arr1);
                    Arrays.sort(arr2);
                    
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                            KL1 = 0; KL2 = 0;
                    
			    for (int k=0; k<matrix1[i][j].length; k++)
		        	{
				    // compute Kullback-Lieber divergence in both directions

                                    //\sum P(i)*Math.log(P(i)/Q(i));
				    KL1 += (matrix1[i][j][k]/count1)*Math.log((matrix1[i][j][k]/count1)/(matrix2[i][j][k]/count2));
                                    //\sum Q(i)*Math.log(Q(i)/P(i));
				    KL2 += (matrix2[i][j][k]/count2)*Math.log((matrix2[i][j][k]/count2)/(matrix1[i][j][k]/count1));                           
			        }
                            
                            if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
				{
                                   sb = new StringBuffer(); 
                                   sb.append(AtomType.type(i)+"\t"+AtomType.type(j)
				             +"\t"+KL1+"\t"+
                                             KL2+"\n");                            
                                   fos.write(sb.toString().getBytes());
				}
              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }

    /***
     Compute mean and variance of distance distributions for each pair of atom types
     ***/
     public boolean computeMeanAndVariance(double[][][] matrix1, double[][][] matrix2, String outFile)
    {
           StringBuffer sb;
           int count1, count2;
           double mean1 = 0; double mean2 = 0;
           double variance1 = 0; double variance2 = 0;

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <mean1> <mean2> <variance1> <variance2>\n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {                    
                    count1 = 0; count2 = 0;
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                    // compute mean and variance when count1 and count2 is known
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                            
                            mean1 = 0; mean2 = 0; 
                    
			    for (int k=0; k<matrix1[i][j].length; k++)
		        	{
				    // compute means
                                    mean1 += (matrix1[i][j][k]/count1)*matrix1[i][j][k];
                                    mean2 += (matrix2[i][j][k]*matrix2[i][j][k])/count2;
			        }

			    
                             for (int k=0; k<matrix1[i][j].length; k++)
			        {
			           variance1 += (matrix1[i][j][k]/count1)*(matrix1[i][j][k]-mean1)*(matrix1[i][j][k]-mean1);
                                   variance2 += (matrix2[i][j][k]/count2)*(matrix2[i][j][k]-mean2)*(matrix2[i][j][k]-mean2);
			        }
                             
                             sb = new StringBuffer(); 
                             sb.append(AtomType.type(i)+"\t"+AtomType.type(j)+
                                       "\t"+mean1+"\t"+mean2
                                       +"\t"+variance1+"\t"+variance2+"\n");
			     fos.write(sb.toString().getBytes());
			}

              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }


    /***
     Compute resampling for each pair of atom types:
         0. Compute absolute difference |M1-M2| of means in each case.
         1. Re-shuffle samples - 2^18 options
         2. Compute means X1, X2 for each reshuffling.
         3. Compute the ratio of reshuffled samples for which |X1-X2|>=|M1-M2|
     ***/
    public boolean computeResampling(double[][][] matrix1, double[][][] matrix2, String outFile, boolean absolute)
    {
           StringBuffer sb;
           int count1, count2;
           double mean1 = 0; double mean2 = 0;
           double mean_diff;
           double rmean1, rmean2, rmean_diff;
           double variance1 = 0; double variance2 = 0;
           int larger_count = 0, totalNumber = 0;
           
           totalNumber = 1<<(RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1);

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <absolute permutation ratio>\n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {                    
                    count1 = 0; count2 = 0;
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                    // compute mean and variance when count1 and count2 is known
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                            
                            mean1 = 0; mean2 = 0;

                            if (AtomType.type(i).toString().startsWith("ICG1") 
                                && AtomType.type(j).toString().startsWith("ICG1")) 
                                {
				    printDistances(matrix1[i][j], "First matrix ");
				    printDistances(matrix2[i][j], "Second matrix ");
				}
                    
			    for (int k=0; k<matrix1[i][j].length; k++)
		        	{
				    // compute original means
                                    mean1 += (matrix1[i][j][k]*matrix1[i][j][k])/count1;
                                    mean2 += (matrix2[i][j][k]*matrix2[i][j][k])/count2;
			        }
                            if (absolute)
                               mean_diff = Math.abs(mean1-mean2);
			    else
                               mean_diff = mean1-mean2;

                            larger_count = 0;

			    // go over all combinations
			     for (int current=0; current<totalNumber; current++)
				 {
				     /*** with the help of a mask, compute reshuffled means***/
                                     rmean1 = 0; rmean2 = 0; 
                                     for (int k=0, mask=0x1; k<RapdfScore.MAX_DISTANCE_BIN-RapdfScore.MIN_DISTANCE_BIN+1; k++, mask = mask <<1)
					 {
                                             int x = current & mask;
                                             //System.out.println("current="+current+" mask="+mask);
					     if (x!=0) // swap matrix1[i][j] and matrix2[i][j]  data
						 {
                                                    // compute the new means
						    rmean2 += (matrix1[i][j][k]*matrix1[i][j][k])/count1;
                                                    rmean1 += (matrix2[i][j][k]*matrix2[i][j][k])/count2;
						 }
					     else // leave the data as is
						 {
                                                    rmean1 += (matrix1[i][j][k]/count1)*matrix1[i][j][k];
                                                    rmean2 += (matrix2[i][j][k]*matrix2[i][j][k])/count2;
						 }
					 } 
                                     if (AtomType.type(i).toString().startsWith("ICG1") 
                                         && AtomType.type(j).toString().startsWith("ICG1") && current <10) 
					    printReshuffledDistances(current, matrix1[i][j], matrix2[i][j], "Reshuffling #"+current);

				     /*** compute new mean difference ***/
                                     
                                     if (absolute)
                                       rmean_diff = Math.abs(rmean1-rmean2);
			             else
                                       rmean_diff = rmean1-rmean2;

                                     if (rmean_diff>=mean_diff)
					 larger_count++;
				 }
			    
                             sb = new StringBuffer(); 
                             sb.append(AtomType.type(i)+"\t"+AtomType.type(j)+
                                       "\t"+((double)larger_count/(double)totalNumber)+"\n");
			     fos.write(sb.toString().getBytes());
                             if (DEBUG_LEVEL>=3)
                               System.out.println("Atom pair ("+i+","+j+") is done, total="+totalNumber+" larger="+larger_count);
			}

              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }


    /***
     Compute Kolmogorov-Smirnov test for each distance distrivution of a pair of atom types
     ***/
     public boolean computeKolmogorovSmirnovTest(double[][][] matrix1, double[][][] matrix2, String outFile)
    {
           StringBuffer sb;
           int count1, count2;
	   double Fn1[], Fn2[]; // for Kolmogorov-Smirnov test
	   double KS1, KS2;
            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <KS1> <KS2> \n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {
                    
                    count1 = 0; count2 = 0;
                    double[] arr1 = new double[matrix1[i][j].length];
                    double[] arr2 = new double[matrix2[i][j].length];
                    arr1 = Arrays.copyOf(matrix1[i][j], matrix1[i][j].length);
                    arr2 = Arrays.copyOf(matrix2[i][j], matrix1[i][j].length);
                    Arrays.sort(arr1);
                    Arrays.sort(arr2);
                    
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                    // compute mean and variance when count1 and count2 is known
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                            // Kolmogorov-Smirnov function is to be computed for every value of arr1[] 
                            Fn1 = new double[arr1.length];
                            Fn2 = new double[arr1.length];

                            // it is easier to compute KS indication reset a sorted array
                            // Fn[x] is constant
		            for (int k=0; k<arr1.length; k++)
		         	{
                                   Fn1[k] = (double)k/(double)arr1.length;
                                   Fn2[k] = (double)k/(double)arr2.length;
                                   
		              	}
		            // now we need to find sup_{x}|Fn[x]-F[x]| where F[x] is the distribution
		            KS1 = 0; KS2 = 0;
                            for (int k=0; k<arr1.length; k++)
				{
                                  if (KS1 < Math.abs(Fn1[k]-matrix1[i][j][k]/count1))
			              KS1 = Math.abs(Fn1[k]-matrix1[i][j][k]/count1);
                                  if (KS2 < Math.abs(Fn2[k]-matrix2[i][j][k]/count2))
			              KS2 = Math.abs(Fn2[k]-matrix2[i][j][k]/count2);
				}
                           
			     sb = new StringBuffer(); 
                             sb.append(AtomType.type(i)+"\t"+AtomType.type(j)
					      +"\t"+
                                              KS1+"\t"+KS2+"\n");
                             fos.write(sb.toString().getBytes());                            
			}

              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }

    /***
     Compute Mann-Whitney test for each pair of atom types
     ***/
    public boolean computeMannWhitneyTest(double[][][] matrix1, double[][][] matrix2, String outFile)
    {
           StringBuffer sb;
           int count1, count2;
	   double MW1, U1, MW2, U2; // for Mann-Whitney test

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <pair count 1> <pair count 2> <median abs difference> <L1> <L2> <Linfty>\n");
              //fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {
                    
                    count1 = 0; count2 = 0;
                    
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                    // compute mean and variance when count1 and count2 is known
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                            U1 = 0; U2 = 0;

                            for (int k=0; k<matrix1[i][j].length; k++)
		         	{
                                   for (int m=0; m<k-1; m++)
				       {
                                         U1 += matrix2[i][j][m]/count2;
                                         U2 += matrix1[i][j][m]/count1;
				       }
                                   U1+=matrix2[i][j][k]/(2*count2);
                                   U2+=matrix1[i][j][k]/(2*count2);
		              	}

			     sb = new StringBuffer(); 
                             sb.append(AtomType.type(i)+"\t"+AtomType.type(j)+U1+"\t"+U2+"\n");
                             fos.write(sb.toString().getBytes());
                             
			}

              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }

    /***
     Compute observed/expected ration for each distance bin and each pair of atom types
     ***/
    public boolean computeObservedVsExpected(double[][][] matrix1, double[][][] matrix2, String outFile,
                           double[][][] padbc1, double[][][] padbc2, double[] pdab1, double[] pdab2,
                           double[] ndab_d1, double[] ndab_d2, double ndab_total1, double ndab_total2)
    {
           StringBuffer sb;
           int count1, count2;

            try
            {
              FileOutputStream fos = new FileOutputStream(outFile); 

              String s = new String ("# Atom1 Atom2 <distance bin> <observed/expected 1> <observed/expected 2>\n");
              fos.write(s.getBytes());
              
              for (int i=0; i<numberOfAtomTypes; i++)
               for (int j=0; j<numberOfAtomTypes; j++) 
                {
                    
                    count1 = 0; count2 = 0;
                    
                    for (int k=0; k<matrix1[i][j].length; k++)
                        {
			    count1 += matrix1[i][j][k];
                            count2 += matrix2[i][j][k];
                        }
                    
                    // compute mean and variance when count1 and count2 is known
                    if (!(count1==18 && count2==18) && !AtomType.type(i).isHydrogen() && !AtomType.type(j).isHydrogen())
			{
                           
			    for (int k=0; k<matrix1[i][j].length; k++)
		        	{
				    // compute observed/expected for each distance value
                                    sb = new StringBuffer(); 
                                    double exp_obs1 = (padbc1[i][j][k]*ndab_total1)/ndab_d1[k];
                                    double exp_obs2 = (padbc2[i][j][k]*ndab_total1)/ndab_d1[k];
                                    sb.append(AtomType.type(i)+"\t"+AtomType.type(j)
					      +"\t"+k+"\t"+exp_obs1+"\t"+
                                              exp_obs2+"\n");
                                    fos.write(sb.toString().getBytes());
			        }
			}

              }
              fos.close();
	      return true;
            }
            catch (IOException e)
            {
                System.out.println("Error while storing results in file "+outFile);
                e.printStackTrace();
		return false;
            }
    }

    /******************************************************************/
    /*                          HELPERS                               */
    /******************************************************************/

    private double roundDouble(double d) {
        	DecimalFormat twoDForm = new DecimalFormat("#.##");
                //System.out.println("Rounding "+d+" to "+twoDForm.format(d));
		return Double.valueOf(twoDForm.format(d));
     }

    private void printDistances(double d[], String header)
    {
        System.out.print(header);
	for (int i=0; i<d.length; i++)
	    System.out.print(" "+d[i]);
	System.out.println();
    }


    private void printReshuffledDistances(int indication, double d1[], double d2[], String header)
    {
        System.out.print(header);
        int mask = 0x1;
	for (int i=0; i<d1.length; i++, mask=mask<<1)
	    if ((indication & mask)!=0)
              System.out.print(" "+d1[i]);
            else
              System.out.print(" "+d2[i]);
	System.out.println();
    }

    /******************************************************************/
    /*                          THE MAIN METHOD                       */
    /******************************************************************/
    public static void main(String[] args) throws IOException{
        
        if (args==null || args.length<3)
        {
            System.out.println("FORMAT: RapdfMatrixGenerationComparison <rapdf matrix file #1> <rapdf matrix file #2> <output file>");
            System.exit(1);
        } 
        RapdfMatrixComparison rmc = new RapdfMatrixComparison(args[0], args[1], args[2]);
    }
}