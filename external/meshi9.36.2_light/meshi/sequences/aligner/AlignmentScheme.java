package meshi.sequences.aligner;

import java.io.IOException;

import meshi.util.file.MeshiLineReader;

public class AlignmentScheme extends SubstitutionMatrix{

    
    public AlignmentScheme(String matrixFilePath, double gapStartPenalty, double inGapPenalty) throws IOException {
        super(gapStartPenalty,inGapPenalty);
		MeshiLineReader matrixFileReader;
		String matrixLine;
	
	    matrixFileReader=new MeshiLineReader(matrixFilePath);
	    
	    matrixLine=matrixFileReader.readLine();//reads order	    
		for (int matLineInd=0,orderInd=0;matLineInd<matrixLine.length();matLineInd++){//collects order
			if (matrixLine.charAt(matLineInd)==' '){
				continue;
			}
			else{
				order[orderInd]=matrixLine.charAt(matLineInd);
				orderInd++;
			}
		}//end of order collection
		
		//make a new order, revOrder, to enable quick access to the contents of the matrix by letters
        setLetterToIndex(order);

		matrixLine=matrixFileReader.readLine();//now for the rest of the matrix
		for (int lineInd=0;lineInd<numAminoAcids+numWildChars;lineInd++){
			if(matrixLine.length()<2){
				matrixLine=matrixFileReader.readLine();
				lineInd--;
				continue;	
			}
			
			matrixLineSplit=matrixLine.split(" ");
			for (int splitInd=0,matrixColInd=0;splitInd<matrixLineSplit.length;splitInd++){	
				if (matrixLineSplit[splitInd].compareTo("")==0||
						(matrixLineSplit[splitInd].compareTo("A")>=0&&matrixLineSplit[splitInd].compareTo("A")<numEnglishLetters)){
					continue;
				}
				substitutionMatrix[lineInd][matrixColInd]=Integer.parseInt(matrixLineSplit[splitInd]);
				matrixColInd++;
				
			}
			matrixLine=matrixFileReader.readLine();
			
		}
		matrixFileReader.close();
		//End of matrix production	
    }

	public double minScore() {
		return 0;
	}
    
}
