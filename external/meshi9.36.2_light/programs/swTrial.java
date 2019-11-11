package programs;

import java.io.IOException;

import meshi.molecularElements.atoms.*;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.*;
import meshi.sequences.AlignmentException;
import meshi.sequences.aligner.AlignmentScheme;
import meshi.util.Rms;
import meshi.util.file.MeshiLineReader;


public class swTrial {

	public static void main(String[] args) throws IOException, AlignmentException {
		AlignmentScheme alSchem=new AlignmentScheme("Blosum62.txt",-10,-0.5);//Added by Tommer 23.9.14
		
		//alSchem.print();	
		
		//routine to test alignment
		
		Protein[] examples=swTrial.CreateProteins();
		//long startTime = System.currentTimeMillis();	
		double[] gdt = Rms.gdt(examples[0], examples[1]);//prints sequences
 		//long estimatedTime = System.currentTimeMillis() - startTime;
		PrintMismatchesAndGaps();
		//System.out.println("\nestimated time is"+estimatedTime);
		//End of routine to test alignment
		
	}
	
	private static Protein[] CreateProteins(){
		//AtomList atomList=new AtomList("..\\work\\QA1\\T0759\\T0759-D1\\T0759-D1.server02_TS1.pdb");
		AtomList atomList=new AtomList("3L22.pdb");
		
		AtomList atomListN=new AtomList("4Q69.pdb");
		//AtomList atomListN=new AtomList("..\\work\\nativeStructures.26.8.14\\T0770_4q69.pdb");
		Protein protein=new Protein(atomList, new ResidueExtendedAtomsCreator());
		Protein proteinN=new Protein(atomListN, new ResidueExtendedAtomsCreator());
		Chain chainN=proteinN.chains().get(0);
		Chain chain=protein.chains().get(0);
		proteinN = new Protein(chainN.atoms(),new ResidueExtendedAtomsCreator());
		protein = new Protein(chain.atoms(),new ResidueExtendedAtomsCreator());
		Protein[] out= {protein, proteinN};
		return out;
	}
	
	private static void PrintMismatchesAndGaps() throws IOException{
		int countMismatches=0, countGapStarts=0, countAllGap=0, inGap=0, inGapTemp;
		double scoreBackTrack=0;
		MeshiLineReader sequences;
		AlignmentScheme matrix=new AlignmentScheme("Blosum62.txt",-10,-0.5);
		String chains[]=new String[2];
		sequences= new MeshiLineReader("..\\work\\workOut\\sequences.txt");
		sequences.readLine();
		chains[0]=sequences.readLine();
		chains[1]=sequences.readLine();
		sequences.close();
		for(int i=0;i<chains[0].length()&&chains[0].charAt(i)!='*';i++){
			//double tempScore=scoreBackTrack;
			if(chains[0].charAt(i)!='-'&&chains[1].charAt(i)!='-'
					&&chains[0].charAt(i)!=chains[1].charAt(i)){
				countMismatches++;
				inGap=0;
				scoreBackTrack+=matrix.substitutionMatrix()
						[matrix.letterToIndex()[chains[0].charAt(i)-65]][matrix.letterToIndex()[chains[1].charAt(i)-65]];//for score count
				//System.out.print("X");
			}
			else{
				if(chains[0].charAt(i)=='-'||chains[1].charAt(i)=='-'){
					
					scoreBackTrack-=0.5;
					countAllGap++;
					inGapTemp=inGap;
					if(chains[0].charAt(i)=='-'){
						inGap=1;//gap in first sequence
					}
					else{
						inGap=2;//gap in second sequence
					}
					if(inGap!=inGapTemp){
						scoreBackTrack-=9.5;
						countGapStarts++;
					}
				}
				else {
					inGap = 0;
					scoreBackTrack += matrix.substitutionMatrix()[matrix.letterToIndex()[chains[0]
							.charAt(i)-65]][matrix.letterToIndex()[chains[1].charAt(i)-65]];// for score count
				}
			}
			//System.out.print("Score addition: "+(scoreBackTrack-tempScore)+" ");
			
		}
		System.out.print("\n"+countMismatches+" mismatches, "+countGapStarts+
				" gap starts, "+countAllGap+" total gaps\n Backtrack score: "+scoreBackTrack);
	}
	
	


}
