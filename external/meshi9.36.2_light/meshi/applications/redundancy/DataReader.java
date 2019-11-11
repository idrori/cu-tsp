package meshi.applications.redundancy;
import java.io.*;
import java.util.*;
import meshi.parameters.ResidueType;
/***
    This class reads the data in PDB or FASTA format serially upon request, and return the current
    FASTA sequence.
 ***/
public class DataReader {

    public static final int    DEBUG_LEVEL = 1; // debug printouts
    public static final String FASTA_FORMAT = "fasta";
    public static final String PDB_FORMAT = "pdb";
    public static final int    PRINT_LINE_LENGTH = 80;
    public static final String PDB_EXTENSION = ".pdb";
    public static final int MAX_LINE = 80;
    public static final  String PDB_SEQRES_RECORD = "SEQRES";
    
    private String  dataPath;   // where to get the data
    private String  dataFormat; // format of the data
    private File[]  filesToRead = null;
    private int     currentFileNumber = -1;
    private LineNumberReader currentLineReader = null;
    
    /*** constructor ***/
    public DataReader (String dataPath, String dataFormat)
        {
            this.dataPath = dataPath;
            this.dataFormat = dataFormat;

            if (DEBUG_LEVEL>=1)
              System.out.println("=========> data path="+dataPath+" data format = "+dataFormat);
        }

    /*** current file being read ***/
    public String getCurrentFileName()
        {
            if (filesToRead!=null && currentFileNumber-1>=0 && currentFileNumber-1<filesToRead.length)
                return filesToRead[currentFileNumber-1].getName();
            return null;
        }

    /*** sequential reader ***/
    public String getNext()
        {
            if (filesToRead == null)
                reset();

            if (DEBUG_LEVEL>=2 && currentFileNumber<filesToRead.length)
                System.out.println("=========> reading "+filesToRead[currentFileNumber].getName());

            LineNumberReader prevLineReader = currentLineReader;
            currentLineReader = getReader();

            try {
            if (currentLineReader != prevLineReader && prevLineReader!=null)
                prevLineReader.close(); }
            catch (IOException e) {}
            
            if (currentLineReader != null)
            {
                if (dataFormat.equals(PDB_FORMAT))
                    return readPdb();
                if (dataFormat.equals(FASTA_FORMAT))
                    return readFasta();
            }

            return null;
        }

    /***
        Read FASTA sequence from a fasta file
    ***/
    private String readFasta()
        {
           try
            {
                String line, fastaSequence = "";
                int mode = 1;

                if (DEBUG_LEVEL>=2)
                    System.out.println("=========> reading file "+filesToRead[currentFileNumber].getName());
                
                while ((line = currentLineReader.readLine()) != null)
                    {
                        System.out.println("read line "+line+" mode = "+mode);
                        if (DEBUG_LEVEL>=1)
                           System.out.println("TESTED => isFasta="+Fasta.isFasta(line));
                        if (Fasta.isFasta(line))
                        {
                            fastaSequence += line;
                            mode = 0;
                        }
                        else
                        {
                            if (mode == 0) // fasta sequence has ended
                            {
                                // query the db with sequences from the query file and blast executable
                                fastaSequence.trim();
                                break;                          
                            }
                            mode = 1;
                         }
                        if (line == null)
                        {
                            currentLineReader.close();
                            currentLineReader = null;
                            currentFileNumber++;
                        }
                }
                if (DEBUG_LEVEL>=2)
                    System.out.println("=========> read fasta sequence "+fastaSequence);
                return fastaSequence;
            }
            catch (IOException e1)
                {
                   System.out.println("Exception "+e1.toString()); 
                }
            return null;    
        }

    /***
        Read FASTA sequence from a pdb file
    ***/
    private String readPdb()
        {
            try
            {
                int lineCounter = 0, j;
                String line, token, oneLetterCode;
                StringTokenizer st;
                StringBuffer fasta = new StringBuffer();
                int total = 0;
                
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
                                    if (total % MAX_LINE == 0)
                                        fasta.append("\n");
                                }
                            }
                            lineCounter++;
                        }

                        currentLineReader.close();
                        currentLineReader = null;
                        currentFileNumber++;
                        fasta.append("\n");
                        return fasta.toString();
            }
            catch (IOException e1)
                {
                   System.out.println("Exception "+e1.toString()); 
                }
            return null;
        }
    

    /*** get current reader ***/
    private LineNumberReader getReader()
        {
            if (currentFileNumber >= filesToRead.length)
                return null;
            if (currentLineReader != null)
                return currentLineReader;
            try {                
                    if (DEBUG_LEVEL>=2)
                       System.out.println("=========> format="+dataFormat+" files num="+filesToRead.length);
                    
                    while (currentFileNumber<filesToRead.length &&
                           (filesToRead[currentFileNumber].isDirectory() ||
                            (dataFormat.equals(PDB_FORMAT)
                             && !filesToRead[currentFileNumber].getName().endsWith(PDB_EXTENSION))))
                        currentFileNumber++;
                    
                    if (DEBUG_LEVEL>=2)
                       System.out.println("=========> currentFileNumber="+currentFileNumber);
                    if (currentFileNumber >= filesToRead.length) return null;
                

                if (DEBUG_LEVEL>=2)
                       System.out.println("=========> got reader");
                return new LineNumberReader(new FileReader(filesToRead[currentFileNumber])); 
            }
            catch (IOException e1)
                {
                   System.out.println("Exception "+e1.toString()); 
                }
            return null;
        }

    /*** start over ****/
    public void reset()
        {
          discoverFiles(dataPath);  
        }

    /*** store files to read in an internal structure ***/
    private void discoverFiles(String path)
        {
                if (DEBUG_LEVEL>=1)
                   System.out.println("=========> Building data reader reset "+path);
                File inp = new File(path);
                if (inp.isDirectory())
                {
                    filesToRead = inp.listFiles();
                    Arrays.sort(filesToRead);
                }
                else
                {
                    filesToRead = new File[1];
                    filesToRead[0] = inp; 
                }
                currentFileNumber = 0;        
        }
}
