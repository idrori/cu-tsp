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
    This programs performs bootstrapping tests of distance distribution between pairs of atom types.
    The program receives following parameters:
        <atom type 1 (3 letters)> <atom type 2 (3 letters)> <distance (3-20)> <PDB database directory> <sample size> <number of samples>
***/
public class BootStrapping{
    public static final String PDB_EXTENSION = ".pdb";
    public static final String PDB_ATOM_RECORD = "ATOM";
    public static final int    DEBUG_LEVEL = 1; // debug printouts

    private String redundancyFile = null;
    private ArrayList redundancy = null;

    // constructor
    public BootStrapping(String atype1, String atype2, int distance, File dbDirectory, 
                         int sampleSize, int numOfSamples, String outFile, String redundancyFile, int sliceSize)
    {
        // initialize the sampling process
        File files[] = dbDirectory.listFiles();
        Arrays.sort(files);

        int sample[] = new int[sampleSize];
        double measurements[] = new double[numOfSamples];
	Arrays.fill(measurements,0);
	Random rand = new Random(System.currentTimeMillis());
        this.redundancyFile = redundancyFile;
        if (redundancyFile != null)
	    {
		System.out.println("Loading redundancy from file "+redundancyFile);
                this.redundancy = new ArrayList();
                boolean rc = RapdfMatrixGeneration.loadRedundancy(this.redundancy, redundancyFile);
                if (!rc) redundancy = null;
                else
                  System.out.println("Loaded redundancy of length "+redundancy.size());
	    }
        
  
	// do the test numOfSamples times
	for (int i=0; i<numOfSamples; i++)
	    {
		generateRandomSample(sample, sampleSize, files.length, rand);
                for (int j=0; j<sample.length; j++)
		    {
			if (DEBUG_LEVEL>=2)
                          System.out.println("Reading file["+sample[j]+"]="+files[sample[j]].getName());
			measurements[i] += measure(files[sample[j]], atype1, atype2, distance, sample[j]);
		    }
                if (DEBUG_LEVEL>=1)
                   System.out.println("sample="+i+" measurements=="+measurements[i]);
	    }
        if (sliceSize<=0)
         saveInFile(measurements, outFile);
        else
	    {
		// compute values distribution
                // find min and max
                double dmin = Integer.MAX_VALUE;
                double dmax = Integer.MIN_VALUE;
                for (int i=0; i<numOfSamples; i++)
		    {
			if (measurements[i]<dmin)
			    dmin = measurements[i];
                        if (measurements[i]>dmax)
			    dmax = measurements[i];
		    }
		int min = (int)Math.floor(dmin);
		int max = (int)Math.ceil(dmax);
                if (DEBUG_LEVEL>=1)
                   System.out.println("min="+min+" max="+max);
                // compute the number of cells
		int cells = (max-min)/sliceSize+2;
                if (DEBUG_LEVEL>=1)
                   System.out.println("cells="+cells);

                // init distribution array
                double distribution[] = new double[cells];
		Arrays.fill(distribution, 0);
                int mincell = min/sliceSize;
                // fill the distribution array
                for (int i=0; i<numOfSamples; i++)
		    {
                        if (DEBUG_LEVEL>=1)
			    System.out.println("value="+measurements[i]+" goes to cell "+(int)Math.floor(measurements[i]/sliceSize));
		        distribution[(int)Math.floor(measurements[i]/sliceSize)-mincell]++; // was += measurements[i];
		    }
                // save array and cell data in file
                saveInFile(min, sliceSize, distribution, outFile);
	    }
    }

    /***
     Save distribution in file
     ***/
    public void saveInFile(int min, int slice, double distribution[], String outFile)
    {
	try {
              FileOutputStream fos = new FileOutputStream(outFile); 
              StringBuffer sb;
              for (int i=0; i<distribution.length; i++, min=min+slice)
               {
                    sb = new StringBuffer();
                    sb.append(min+"\t"+distribution[i]+"\n");
                    fos.write(sb.toString().getBytes());
                }
              fos.close();
	}
         catch (IOException ioe)
		     {
			ioe.printStackTrace();
			System.exit(0);
		    }
    }

    /***
     Save results in file
     ***/
    public void saveInFile(double measurements[], String outFile)
    {
	try {
              FileOutputStream fos = new FileOutputStream(outFile); 
              StringBuffer sb;
              for (int i=0; i<measurements.length; i++)
               {
                    sb = new StringBuffer();
                    sb.append(measurements[i]+"\n");
                    fos.write(sb.toString().getBytes());
                }
              fos.close();
	}
         catch (IOException ioe)
		     {
			ioe.printStackTrace();
			System.exit(0);
		    }
    }

    /***
     Generate random numbers in range 0..max and place then in an array
     ***/
    public void generateRandomSample(int arr[], int size, int max, Random generator)
       {
	   for (int i=0; i<size; i++)
	       {
		   arr[i] = generator.nextInt(max);
	       }
           if (DEBUG_LEVEL>=2)
	       {
                 System.out.print("generated sample ");
                 for (int i=0; i<size; i++) System.out.print(arr[i]+" ");
	         System.out.println();
	       }
       }

    /***
     Read file in PDB format and compute the number of getAtoms of types atom1 and atom2 at distance d
     ***/
    public double measure(File pdbFile, String atom1, String atom2, int distance, int pdbFileNumber)
    {
	if (pdbFile.exists() && pdbFile.getName().endsWith(PDB_EXTENSION))
	    {
                try  { // read the PDB file
                 
                 LineNumberReader pdbFileReader = new LineNumberReader(new FileReader(pdbFile));
                 String line, token, atomToken, previousResidue = "";
                 StringTokenizer st;
		 int residueNumber = 0;
                 ArrayList atoms = new ArrayList();
                 ArrayList coordinatesList = new ArrayList();
                 ArrayList residues = new ArrayList();
                 ArrayList Xcoordinate = new ArrayList();
		
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
                                    
                                    int place = Arrays.binarySearch(Xcoordinate.toArray(), Xcoor);
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
                                }
                                catch (NumberFormatException e2)
                                {
                                }
			    } // and if ATOM_RECORD
			} // end while
                  pdbFileReader.close();
		  /*** now go over atom pairs and search for specified types and distance ***/
                  return computeNumberOfAppearances(atoms, residues, coordinatesList, atom1, atom2, distance, pdbFileNumber);
		}
	       catch (IOException ioe)
		     {
			ioe.printStackTrace();
			System.exit(0);
		    }
	    }
	return 0;
    }

    /***
        Compute redundancy measure (weight) for a residue in the specific sequence (file)
     ***/
    private double getRedundancyWeightFromMemory(int res, int sequence)
        {
            if (this.redundancy==null || sequence<0 || sequence>=this.redundancy.size())
                return 1;
            
            double w[] = (double[])this.redundancy.get(sequence);
            if (res<0 || res>=w.length)
                return 1;

            if (w[res]<0)
                return 1;
            
            if (DEBUG_LEVEL>=2)
                System.out.println("++++++> redundancy weight ="+w[res]);

            return w[res];
        }


    /*** Find in the atom list pairs of residues of required types and at the required distance
         and count their number while taking redundancy in consideration (if requested). ***/
    public double computeNumberOfAppearances(ArrayList atoms, ArrayList residues, ArrayList coordinatesList,
                                             String type1, String type2, int requiredDistance, int pdbFileNumber)
    {
	 boolean first;
         String atom1, atom2;
	 double finalCount=0;
         double[] coor1, coor2;
         double delta_count_a, delta_count_b;
         String selectAttribute = "redu";
         String whereCondition;

         for (int j=0; j<atoms.size(); j++)
             {
                 // find the first required atom of the pair
                 atom1 = ((AtomType)atoms.get(j)).toString();
		 if (!atom1.equals(type1) && !atom1.equals(type2))
		     continue;

                 // atom found
                 first = true;
                 if (!atom1.equals(type1))
		     first = false;

                 // get the redundancy weight of its residue
                 delta_count_a = getRedundancyWeightFromMemory(((Integer)residues.get(j)).intValue(), pdbFileNumber);
                 if (DEBUG_LEVEL>=2)
		    System.out.println("redudancy weight for file="+pdbFileNumber+" and residue="
                                       +((Integer)residues.get(j)).intValue()+" is "+delta_count_a);

                 coor1 = (double[])coordinatesList.get(j);                          
                            
                 for (int k=j+1; k<atoms.size(); k++)    
                     {
                        // find the second required atom of the pair
			atom2 = ((AtomType)atoms.get(k)).toString();
                        if (!((first && atom2.equals(type2)) || (!first && atom2.equals(type1))))
		          continue;
 
                        coor2 = (double[])coordinatesList.get(k);

                                 
                        if (Math.abs(coor1[RapdfScore.X]-coor2[RapdfScore.X])>RapdfScore.MAX_DISTANCE_BIN+1)
                           {
                             k = atoms.size();
                             continue;
                           }                      

                         // compute the distance if it is not too large
                         if ((((Integer)residues.get(j)).intValue() != ((Integer)residues.get(k)).intValue())
                                     && Math.abs(coor1[RapdfScore.Y]-coor2[RapdfScore.Y])<=RapdfScore.MAX_DISTANCE_BIN+1
                                     && Math.abs(coor1[RapdfScore.Z]-coor2[RapdfScore.Z])<=RapdfScore.MAX_DISTANCE_BIN+1)
                                 {
                                     double distance = Math.sqrt((coor1[RapdfScore.X]-coor2[RapdfScore.X])*(coor1[RapdfScore.X]-coor2[RapdfScore.X])+
                                        (coor1[RapdfScore.Y]-coor2[RapdfScore.Y])*(coor1[RapdfScore.Y]-coor2[RapdfScore.Y])+
                                        (coor1[RapdfScore.Z]-coor2[RapdfScore.Z])*(coor1[RapdfScore.Z]-coor2[RapdfScore.Z]));
                                     int roundDistance = (int)Math.floor(distance);
				     if (DEBUG_LEVEL>=3)
                                      System.out.println("Found atom pair "+atom1+"-"+atom2+" at distance "+roundDistance);
                                     if (roundDistance==requiredDistance)
					 {
                                             delta_count_b = getRedundancyWeightFromMemory(((Integer)residues.get(k)).intValue(), pdbFileNumber);
                                             
					     finalCount += delta_count_a*delta_count_b;
					 }
                                     
                                 }
                             }
                        }
                 // clearRows lists
                 atoms.clear();
                 residues.clear();
                 coordinatesList.clear();
		 return finalCount;
    }

    /******************************************************************/
    /*                          THE MAIN METHOD                       */
    /******************************************************************/
    public static void main(String[] args) throws IOException{
        
        if (args==null || args.length<7)
        {
            System.out.println("FORMAT: PdbDatabaseSeparation <source directory> <target directory> <number of residues> <out file>");
            System.exit(1);
        } 
        try {
            String atype1 = args[0].toUpperCase();
            String atype2 = args[1].toUpperCase();
            if (atype1.length()>4 || atype2.length()>4)
		{
                   System.out.println("Atom type should consist of three-four letters while current lengths are "+atype1.length()+" and "+atype2.length());
                   System.exit(1);
		}
            int distance = Integer.parseInt(args[2]);
            if (distance<3 || distance>20)
                {
                   System.out.println("Distance should lie in the interval 3..20");
                   System.exit(1);
		}
            File dbDirectory = new File(args[3]);
            if (!dbDirectory.exists() || !dbDirectory.isDirectory())
                 {
                   System.out.println("Database direcotry "+args[3]+" does not exist or is not a direcotry");
                   System.exit(1);
		}

            int sampleSize = Integer.parseInt(args[4]);
            if (sampleSize<0)
                {
                   System.out.println("Sample size should be greater than 0");
                   System.exit(1);
		}

            int numOfSamples = Integer.parseInt(args[5]);
            if (numOfSamples<0)
                {
                   System.out.println("Number of samples should be greater than 0");
                   System.exit(1);
		}
            
            /*** if redundancy table is specified, compute with redundancy ***/
            String redundancyTableName = null;
            if (args.length>7 && !args[7].equals("none"))
		redundancyTableName = args[7];

            int sliceSize = -1;
            if (args.length>8)
		sliceSize = Integer.parseInt(args[8]);

	    BootStrapping bs = new BootStrapping(atype1, atype2, distance, dbDirectory, 
						 sampleSize, numOfSamples, args[6], redundancyTableName, sliceSize);
	}
	catch (NumberFormatException e)
	    {
		e.printStackTrace();
	    }
    }
}