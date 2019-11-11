package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ExtendedAtomsProtein;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentException;
import meshi.util.Histogram;
import meshi.util.Rms;
import meshi.util.Stat;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.formats.Format;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 */
public class Consensus {
    /*
    * In the set of CASP9_FM native structures there are 1560 helix residues, 816 beta & 1763 coil
    * That is 0.74 is the ratio C/(H+E)
     */
    public static final double CHE = 0.75;
    public static final int SEPARATION = 2;
    public static int MAX_CONSENSUS_RESIDUES = 1000;
    public static double CONSENSUS_FRACTION = 0.5;
    public static double MINIMAL_NUMBER = 10;
    public static final int MAX_LENGTH = 600; // maximum length of protein domain
    private static Stat statAll = new Stat(0,1,0.1);
    private static double[] referenceStd = new double[MAX_LENGTH];

    public static void main(String[] argv) throws IOException,AlignmentException{
        ContactMap contactMap;
        ArrayList<ConsensusModel> models;
        if(argv.length != 2)
            throw new RuntimeException("Usage:Consensus prefix nativeExists/noNative");
        String prefix = argv[0];
        boolean nativeFlag = argv[1].equals("nativeExists");
        Protein nativeStructure;
        getReferenceStd(prefix);
        ArrayList<File> targets  = Utils.getDirectories(prefix);
        for (File target : targets) {
            System.out.println("Processing target "+target.getName());
            if (nativeFlag) nativeStructure = getNativeStructure(target);
            else nativeStructure = null;
            contactMap = getContactMap(target);
            models = sortModels(target, nativeStructure, contactMap, CONSENSUS_FRACTION);
            contactMap = getContactMap(target, models);
            printContactMap(target, contactMap);
            cluster(target,models);
        }
        statAll.getHistogram().print();
        System.out.println("Max & average fraction of best model "+statAll.getMax()+" "+statAll.getMean());
    }

    public static void cluster(File target, ArrayList<ConsensusModel> models) throws IOException{
        MeshiWriter clusterWriter = new MeshiWriter(target.getName()+".clusters");
        ArrayList<ConsensusModel> usedClusters;

    }
    public static void getReferenceStd(String prefix) throws IOException{
        Reference[] references;
        references = new Reference[MAX_LENGTH];

        File referenceFile = new File("referenceStd.dat");
        String line;
        if (referenceFile.exists()) {
            MeshiLineReader reader = new MeshiLineReader("referenceStd.dat");
            while ((line = reader.readLine()) != null) {
                String[] words = line.split(",");
                if (words.length == 2) {
                    int index = Integer.valueOf(words[0].trim());
                    double std = Double.valueOf(words[1].trim());
                    System.out.println("read reference: "+index+" "+std);
                    referenceStd[index]= std;
                }
            }
            for (int i = 0; i < referenceStd.length; i++) {
                if (referenceStd[i]<= 0) referenceStd[i] = 1;
            }
        }
        else {
            ArrayList<File> targets  = Utils.getDirectories(prefix);
            for (File target : targets) {
                System.out.println("reference from target "+target.getName());
                Protein model;
                ArrayList <File> modelFiles = getModelFiles(target);
                for (int iModel = 0; iModel < modelFiles.size(); iModel++) {
                   File modelFile = modelFiles.get(iModel);
                   System.out.println("reference from "+target.getName()+" "+modelFile.getName());
                   try {
                     model = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(modelFile);
                   }
                   catch (Exception e) { System.out.println("Failed"); continue;}
                   AtomList caAtoms = model.atoms().CAFilter();
                   for (int i = 0; i <caAtoms.size(); i++)
                       for(int j = i + SEPARATION ; j < caAtoms.size(); j++) {
                           Atom iAtom =  caAtoms.get(i);
                           Atom jAtom =  caAtoms.get(j);
                           int iResidueNumber = iAtom.residueNumber();
                           int jResidueNumber = jAtom.residueNumber();
                           double distance = iAtom.distanceFrom(jAtom) ;
                           int separation = jResidueNumber-iResidueNumber;
                           if (separation < 0) throw new RuntimeException("This is weird");
                           if (references[separation] == null) references[separation] = new Reference();
                           references[separation].add(distance);
                       }
                }
            }
            MeshiWriter writer = new MeshiWriter(referenceFile);
            for (int i = 0; i < references.length; i++){
                Reference reference = references[i];
                if (reference != null)
                    writer.println(String.format("%-10d , %-10.4f",i,reference.getStd()));

            }
            writer.close();
            getReferenceStd(prefix);
        }
    }
    public static Protein getNativeStructure(File target) {
        ArrayList<File> files = getFiles(target,true);
        if (files.size() != 1) throw new RuntimeException("Failed to find native structure for "+target.getName());
        File nativeFile = files.get(0);
        Protein nativeStructure =  ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(nativeFile);
        System.out.println("Native structure for "+target.getName()+": "+nativeFile.getName());
        return nativeStructure;
    }

    /*public static ArrayList<File> getModelFiles(File target) {
        return getModelFiles(target,null);
    } */
    public static ArrayList<File> getModelFiles(File target) {
          ArrayList<File> out =  getFiles(target,false);
         //for (File file : out) System.out.println(file.getName());
         return out;
    }

    public static ArrayList<File> getFiles(File target,boolean nativeFlag) {
        ArrayList out = new ArrayList<File>();
        File[] files = target.listFiles();
        for (File file : files) {
            if (file.getName().endsWith("pdb")) {
                if (nativeFlag)  {
                    if (file.getName().indexOf(".N.pdb") >= 0) {
                        out.add(file);
                        break;
                    }
                }
                else {
                    if ((file.getName().indexOf(".N.pdb") < 0) &&
                            (file.getName().indexOf("native") < 0)) {
                        out.add(file);
                    }
                }
            }
        }
        return out;
    }

    public static ContactMap getContactMap(File target) {
        return getContactMap(target,null);
    }
    public static ContactMap getContactMap(File target, ArrayList<ConsensusModel> models) {
           Protein model;
        ArrayList <File> modelFiles = getModelFiles(target);
        ContactMap contactMap = new ContactMap();
           for (int iModel = 0; iModel < modelFiles.size(); iModel++) {
                File modelFile = modelFiles.get(iModel);
                System.out.println(target.getName()+" "+modelFile.getName());
                try {
                    model = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(modelFile);
                   // Utils.AssignDSSP(model, modelFile.getAbsolutePath());
                }
                catch (Exception e) {
                    System.out.println("Failed");
                    //e.printStackTrace();
                    //throw new RuntimeException(e);
                    continue;
                }
                for (Residue residue : model.residues()) {
                    ConsensusResidue consensusResidue = contactMap.getResidue(residue.number());
                    consensusResidue.addResidue(residue.type().nameOneLetter(), residue.getSecondaryStructure());
                }
                AtomList caAtoms = model.atoms().CAFilter();
                for (int i = 0; i <caAtoms.size(); i++)
                    for(int j = i + SEPARATION ; j < caAtoms.size(); j++) {
                        Atom iAtom =  caAtoms.get(i);
                        Atom jAtom =  caAtoms.get(j);
                        contactMap.add(iAtom.residueNumber(),jAtom.residueNumber(),iAtom.distanceFrom(jAtom)) ;
                    }
            }
        for (int i = contactMap.minResidue(); i<= contactMap.maxResidue();i++) {
            ConsensusResidue iResidue = contactMap.getResidue(i);
              for (int j = i+2; j <= contactMap.maxResidue(); j++)   {
                  ContactMapCell cell = contactMap.getCell(i,j);
                  if (cell.numberOfObservations() >= MINIMAL_NUMBER) {
                      int separation = j-i;
                      double score = -Math.log(cell.getStd()/referenceStd[separation]);
                      if (Utils.isAbnormal(score)) {
                          throw new RuntimeException("Abnormal score "+score+" "+cell.getStd()+" "+
                                                      separation+" "+referenceStd[separation]+" "+i+" "+j);
                      }
                      cell.setScore(score);
                      contactMap.cells.add(cell);
                      ConsensusResidue jResidue = contactMap.getResidue(j);
                      iResidue.add(score);
                      jResidue.add(score);
                  }
              }
        }
        return contactMap;
    }

    public static void printContactMap(File target, ContactMap contactMap) throws IOException{
        MeshiWriter cmWriter = new MeshiWriter(target.getName()+".CM");
        MeshiWriter consensusWriter = new MeshiWriter(target.getName()+".consensus");
        MeshiWriter secondaryStructureWriter = new MeshiWriter(target.getName()+".ss");
        for (int i = contactMap.minResidue(); i<= contactMap.maxResidue();i++) {
            ConsensusResidue iResidue = contactMap.getResidue(i);
              for (int j = i+2; j <= contactMap.maxResidue(); j++)   {
                  ContactMapCell cell = contactMap.getCell(i,j);
                  if (cell.numberOfObservations() >= MINIMAL_NUMBER) {
                      int separation = j-i;
                      double score = -Math.log(cell.getStd()/referenceStd[separation]);
                      cell.setScore(score);
                      if (Utils.isAbnormal(score)) {
                                      throw new RuntimeException("Abnormal score "+score+" "+cell.getStd()+" "+
                                                                  separation+" "+referenceStd[separation]+" "+i+" "+j);
                                  }
                      contactMap.cells.add(cell);
                      ConsensusResidue jResidue = contactMap.getResidue(j);
                      iResidue.add(score);
                      jResidue.add(score);
                  }
              }
        }
        Object[] cellsArray = contactMap.cells.toArray();

        Arrays.sort(cellsArray);
        for (int i = 0; i < cellsArray.length; i++) {
            ContactMapCell cell = (ContactMapCell) cellsArray[i];
            cmWriter.printf("%-10d %-10d  %12.6f %12.6f %12.6f %-10d\n",cell.iResidue,cell.jResidue,cell.getScore(),
                             cell.getMean(),cell.getStd(),cell.numberOfObservations());
        }
        cmWriter.close();
        ConsensusResidue[] visited = visitedResidues(contactMap.residues);

        double[] weights = getNonSsWeight(visited);
        double ssWeight = weights[0];
        double nonSsWeight = weights[1];
        for (int i = 0; i < visited.length;i++) {
            secondaryStructureWriter.println(visited[i].secondaryStructureToString(ssWeight,nonSsWeight));
        }
        secondaryStructureWriter.close();

        Arrays.sort(visited);
        for (int i = 0; i < visited.length;i++) {
            consensusWriter.println(visited[i]);
        }
        consensusWriter.close();
    }

    public static double[] getNonSsWeight(ConsensusResidue[] residues){
        int nSS = 0;
        int nNonSS = 0;
        for (ConsensusResidue residue:residues){
            nSS += residue.secondaryStructure.getBin(0)+residue.secondaryStructure.getBin(1);
            nNonSS += residue.secondaryStructure.getBin(2);
        }
        double fSS = 1.0*nSS/(nSS+nNonSS);
        double fNonSS = 1.0*nNonSS/(nSS+nNonSS);
        double ssWeight = 1/(fSS*(CHE+1));
        double nonSsWeight = CHE/(fNonSS*(CHE+1));
        double out[] = {ssWeight,nonSsWeight};
        return out;
    }

    public static ConsensusResidue[] visitedResidues(ConsensusResidue[] residues) {
        int i = 0;
        for (ConsensusResidue residue : residues) if (residue.numberOfObservations()>= MINIMAL_NUMBER)i++;
        ConsensusResidue[] out = new ConsensusResidue[i];
        i = 0;
        for (ConsensusResidue residue : residues) if (residue.numberOfObservations()>= MINIMAL_NUMBER){
            out[i] = residue;
            i++;
        }
        return out;
    }
    public static ArrayList<ConsensusModel> sortModels(File target,
                                                       Protein nativeStructure,
                                                       ContactMap contactMap,
                                                       double fraction) throws IOException,AlignmentException{
        Protein model;
        double Q;
        Atom iAtom,jAtom;
        Stat stat = new Stat(0,1,0.01);
        MeshiWriter modelsWriter = new MeshiWriter(target.getName()+".models");
        ArrayList<ConsensusModel> models = new ArrayList<ConsensusModel>();
        ArrayList <File> modelFiles = getModelFiles(target);
        double[] gdt = {-1,-1,-1,-1,-1};
        Stat realStat = new Stat();
        for (int iModel = 0; iModel < modelFiles.size(); iModel++) {
                Q = 0;
                File modelFile = modelFiles.get(iModel);
                System.out.println(target.getName()+" testing "+modelFile.getName());
                try {
                    model = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(modelFile);
                }
                catch (Exception e) {
                    System.out.println("Failed");
                    continue;
                }
                double sum = 0;
                AtomList caAtoms = model.atoms().CAFilter();
                if (nativeStructure != null)
                  gdt = Rms.gdt(nativeStructure, model );
                if (gdt[0] >= 1 ) System.out.println("************************ GDT 1 ***************************************************");
                realStat.add(gdt[0]);
                for (int i = 0; i <caAtoms.size(); i++) {
                    iAtom =  caAtoms.get(i);
                    for(int j = i + SEPARATION ; j < caAtoms.size(); j++) {
                        jAtom =  caAtoms.get(j);
                        ContactMapCell cell =  contactMap.getCell(iAtom.residueNumber(), jAtom.residueNumber());
                        double weight = Math.exp(10*cell.getScore());
                        Q += weight;
                        double z = weight * cell.zScore(iAtom.distanceFrom(jAtom)) ;
                        sum += Math.abs(z);
                    }
                }
             models.add(new ConsensusModel(gdt[0],sum/Q,modelFile.getName()));
            }
        Object[] modelsArray = models.toArray();
        Arrays.sort(modelsArray);
        models = new ArrayList<ConsensusModel>();
        for (int i = 0; i < modelsArray.length; i++)
            models.add((ConsensusModel) modelsArray[i]);


        double modelNumber = 1;
        System.out.println("# columns: "+ConsensusModel.header()+" cumulativeMeanGdt cumulativeMaxGdt meanGdt maxGdt fraction");
        for (ConsensusModel consensusModel : models) {
            stat.add(consensusModel.gdt);
            if (consensusModel.gdt == realStat.getMax()) statAll.add(modelNumber/models.size());
             modelsWriter.printf("%s %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n",
                                   consensusModel.toString(),
                                   stat.getMean(),stat.getMax(),
                                   realStat.getMean(),realStat.getMax(),
                                   modelNumber/models.size());
            modelNumber += 1;
        }
        modelsWriter.close();
        return models;
    }
       //-------------------------------------------- private classes  ---------------------------------------------------------------------
    private static class ConsensusResidue extends Stat {
        private int number;
        String name;
        private Histogram secondaryStructure;

        public ConsensusResidue(int number) {
            this.number = number;
            if (number < 0) add(10000);
            secondaryStructure = new Histogram(0,2,1);
        }

        public String toString() {return String.format("%-10d %12.6f %12.6f %-20d",number,getMean(),getStd(),numberOfObservations());  }

        public int compareTo(Object o) { // Reverses the order. Want high scoring residues first.
            return -super.compareTo(o);
        }

        public void add(double addMe) {
            super.add(addMe);
            double mean = getMean();
            if (Utils.isAbnormal(mean)) {
                            throw new RuntimeException("Abnormal score "+mean+" "+addMe);
                        }

            setScore(getMean());
        }

           public void addResidue(String name, SecondaryStructure secondaryStructure) {
               this.secondaryStructure.add(secondaryStructure.ordinal());
               this.name = name;
           }

           public String secondaryStructureToString(double ssWeight, double nonSsWeight){
               double[] probabilities = secondaryStructure.fractionalHistogram();
               double max = -1;
               int maxIndex = -1;
               SecondaryStructure ss;
               double weightH = ssWeight*probabilities[0];
               double weightE = ssWeight*probabilities[1];
               double weightC = nonSsWeight*probabilities[2];
               double sum = weightC+weightE+weightH;
               weightC = weightC/sum;
               weightE = weightE/sum;
               weightH = weightH/sum;
               if (weightH > weightE) {
                   if(weightH>weightC) { ss = SecondaryStructure.HELIX;}
                   else                { ss = SecondaryStructure.COIL;}
               }
               else {
                   if(weightE>weightC) { ss = SecondaryStructure.SHEET;}
                   else                { ss = SecondaryStructure.COIL;}
               }

               String out = name+" "+
                            Format.fintL(number,10)+" "+ss+" "+
                            Format.fdoubleL(weightH,3,10)+" "+
                            Format.fdoubleL(weightE,3,10)+" "+
                            Format.fdoubleL(weightC,3,10);
               return out;
           }
    }
    private static class ContactMapCell extends Stat {
        private int iResidue,jResidue;
        public ContactMapCell(int iResidue, int jResidue) {
            super(3,10,0.5);
            this.iResidue = iResidue;
            this.jResidue = jResidue;
        }

        public int compareTo(Object o){ // Reverses the order. We want high score contacts first.
            return -super.compareTo(o);
        }
        public void add(double distance) {
            if (distance < 0) throw new RuntimeException("Negative distance");
            super.add(distance);
        }
        public double getMean() {
            if (numberOfObservations() < MINIMAL_NUMBER) return -1;
            return super.getMean();
        }
        public double getStd() {
            if (numberOfObservations() < MINIMAL_NUMBER) return -1;
            return super.getStd();
        }
    }

    private static class ContactMap {
        private ContactMapCell[][] matrix;
        private ConsensusResidue[]  residues;
        private ArrayList<ContactMapCell>cells;
        private static final ConsensusResidue dummyResidue = new ConsensusResidue(-1);
        private int minResidue, maxResidue;
        public int minResidue() {return minResidue;}
        public int maxResidue() {return maxResidue; }

        public ContactMap() {
            minResidue = 1000;
            maxResidue = -1000;
            matrix = new ContactMapCell[1][1];
            residues = new ConsensusResidue[MAX_CONSENSUS_RESIDUES];
            for (int i = 0 ; i < MAX_CONSENSUS_RESIDUES; i++) residues[i] = dummyResidue;
            cells = new ArrayList<ContactMapCell>();
        }
        public int size() {return matrix.length;}

        public ConsensusResidue getResidue(int i) {
            if (i > MAX_CONSENSUS_RESIDUES) throw new RuntimeException("weird i");
            if (residues[i] == dummyResidue) residues[i] = new ConsensusResidue(i);
            return residues[i];
        }
        private ContactMapCell[][] copyMatrix(int newSize, int newMin) {
            ContactMapCell[][] newMatrix = new ContactMapCell[newSize][newSize];
            for (int i = minResidue; i <= maxResidue; i++)
                for (int j = minResidue; j <= maxResidue; j++)  {
                    newMatrix[i-newMin][j-newMin] = matrix[i- minResidue][j- minResidue];
                }
            return newMatrix;
        }


        public ContactMapCell getCell(int i, int j) {
            if (matrix[i- minResidue][j- minResidue] == null)
                matrix[i- minResidue][j- minResidue] = new ContactMapCell(i,j);
            return matrix[i- minResidue][j- minResidue];
        }
        public void add(int iRow, int jColumn, double distance) {
            int newMin;
            int newMax;
            boolean change = false;
            if (iRow < minResidue) {
                newMin = iRow;
                change = true;
            }
            else newMin = minResidue;
            if (jColumn < newMin) {
                change = true;
                newMin = jColumn;
            }
            if (iRow > maxResidue) {
                newMax = iRow;
                change = true;
            } else newMax = maxResidue;
            if(jColumn > newMax) {
                change = true;
                newMax = jColumn;
            }
            if (change) {
                int newSize = newMax-newMin+1;
                matrix = copyMatrix(newSize,newMin);
                minResidue = newMin;
                        maxResidue = newMax;
            }
            ContactMapCell cell = getCell(iRow,jColumn);
            cell.add(distance);
        }
    }
    private static class Reference extends Stat{
        public Reference() {
            super();
        }
    }

    public static class ConsensusModel implements Comparable {
        public final double score;
        public final String name;
        public final double gdt;

        public ConsensusModel(double gdt, double score, String name) {
            this.gdt = gdt;
            this.name = name;
            this.score = score;
        }
         public int compareTo(Object o) {
             ConsensusModel other = (ConsensusModel) o;
             if (score > other.score) return 1;
             if (score < other.score) return -1;
             return 0;
         }
        public String toString() {
            return String.format("%-35s %-10.4f %-10.4f", name,score,gdt);
        }
        public static String header() {return "name score gdt";}
    }

}