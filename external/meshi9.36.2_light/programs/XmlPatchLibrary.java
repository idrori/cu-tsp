package programs;

import meshi.util.file.MeshiWriter;
import meshi.util.info.MeshiInfoElementList;
import meshi.util.info.MeshiInfoXMLwriter;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * @author Avia Ohana
 *
 */
public class XmlPatchLibrary {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		try {
		    FileInputStream fstream = new FileInputStream(args[0]);
		    DataInputStream in = new DataInputStream(fstream);
		    BufferedReader br = new BufferedReader(new InputStreamReader(in));
            MeshiWriter writer = new MeshiWriter(args[0]+".xml");

            writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>");
            writer.println("<PatchLibrary size=\"350\">");
            String strLine;
		
			// reads until we get to the pivots
			while ( ! (strLine = br.readLine()).startsWith("1)") )  {
                System.out.println("1 " + strLine);
            }
			
			for(int i=1; i<=350; i++) {
				System.out.println("2 "+strLine);
				String[] splited = strLine.split(" ");
				String pivotName = splited[3];
				br.readLine();            //read the line "patch getAtoms"
                PatchLibraryElement patchLibraryElement = new PatchLibraryElement(pivotName);
				int atomsCounter = 1;
				
				// read all pivot's getAtoms
				while ( (strLine = br.readLine()).startsWith(""+i+ "") )  {
				     System.out.println("3 "+strLine);
					// strLine: 1.1) Residue Name = HIS	 Residue Number = 62	 
					//          0    1       2    3 4    5       6      7 8
					//Atom Name = NE2	Coordinates (14.942)(-28.83)(-11.252)
					//9    10  11 12                 13
					
					String[] atomArgs = strLine.split(" ");
					atomArgs[12] = atomArgs[12].replace(')', ' ').replace('(',' ');
					String[] atomCoors = atomArgs[12].split(" ");
                    String x = atomCoors[1];
                    String y = atomCoors[3];
                    String z = atomCoors[5];

					String residueName = (atomArgs[4].split("\t"))[0];
                    String residueNumber = (atomArgs[7].split("\t"))[0];
                    String type = atomArgs[10];
                    patchLibraryElement.atoms.add(new PatchAtom(residueName, residueNumber, type, x,y,z));
					atomsCounter++;
				
				} //while
                writer.println(patchLibraryElement.toXml());
			
			} //for

            writer.println("</PatchLibrary>");
			in.close();
            writer.close();
		}
		catch (IOException e) {
			System.out.println("error!\n" + e.toString());
		}

	}

    private static class PatchLibraryElement{
        String pivot;
        ArrayList<PatchAtom> atoms = new ArrayList<PatchAtom>();

        public PatchLibraryElement(String pivot) {
            this.pivot = pivot;
        }
        public String toXml() {
            String out = "\t<PatchLibraryElement pivot=\""+pivot+"\">\n";
             for (PatchAtom patchAtom : atoms)
                 out+=patchAtom.toXML()+"\n";
            out +=  "\t</PatchLibraryElement>";
            return out;
        }
    }

    private static class PatchAtom {
        String residueName;
        String residueNumber;
        String x,y, z;
        String atomName;
        public PatchAtom(String residueName, String residueNumber, String atomName, String x, String y, String z){
            this.residueName = residueName;
            this.residueNumber = residueNumber;
            this.x = x;
            this.y = y;
            this.z = z;
            this.atomName = atomName;
        }
        public String toXML() {
            return "\t\t<PatchAtom residueName=\""+residueName+"\" residueNumber=\""+residueNumber+"\"  atomName=\""+atomName+"\" x=\""+x+"\""+" y=\""+y+"\""+" z=\""+z+"\"/>";
        }
    }
}
