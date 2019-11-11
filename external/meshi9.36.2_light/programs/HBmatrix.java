package programs;

import com.sun.org.apache.regexp.internal.RE;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomBuilder;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class HBmatrix {
    public static void main(String[] args) throws IOException{
        double[][] hArray, oArray;
        double[] nowhere = {-9999, -9999, -9999};
        File file = new File(args[0]);
        Protein protein = new Protein(new AtomList(file.getAbsolutePath()), ResidueExtendedAtoms.creator);
        Chain chain = protein.chain();
        JSONObject jsonObject = new JSONObject();
        int first = chain.firstNonDummyResidue().number();
        hArray = new double[chain.size() - first][];
        oArray = new double[chain.size() - first][];

        jsonObject.put("file", file.getAbsolutePath());
        jsonObject.put("firstResidue", chain.firstNonDummyResidue().number());
        JSONArray hBonds = new JSONArray();
        jsonObject.put("hBonds", hBonds);

        hArray[0] = nowhere; // Cannot locate without previous residue
        Residue prevResidue = chain.residueAt(first);
        Atom atom = prevResidue.carbonylO();
        double[] firstO = {atom.x(), atom.y(), atom.z()};
        oArray[0] = firstO;
        for (int i = 0; i < oArray.length; i++) {
            Residue residue = chain.residueAt(i+first);
            if (residue.dummy()) {
                hArray[i] = nowhere;
                oArray[i] = nowhere;
            } else {
                Atom oAtom = residue.carbonylO();
                if (oAtom.nowhere())
                    oArray[i] = nowhere;
                else {
                    double[] o = {oAtom.x(), oAtom.y(), oAtom.z()};
                    oArray[i] = o;
                }
                Atom hAtom = residue.amideH();
                Atom cAtom = prevResidue.carbonylC();
                Atom nAtom = residue.amideN();
                Atom caAtom = residue.ca();

                if ((hAtom != null) && (!caAtom.nowhere()) && (!nAtom.nowhere()) && (!cAtom.nowhere())) {
                    //$H1=gen_atom_LAT(PI,$H_N_CA_ang,$H_N_len,$N,$CA,$pC);
                    AtomBuilder.buildByTorsion(Math.PI, AtomBuilder.H_N_CA_ang, 0.95, nAtom, caAtom, cAtom, hAtom);
                    double[] h = {hAtom.x(), hAtom.y(), hAtom.z()};
                    hArray[i] = h;
                } else {
                    hArray[i] = nowhere;
                }
            }
            prevResidue = residue;
        }

        double[] h, o;
        for (int i = 0; i < hArray.length; i++) {
            if (hArray[i][0] > -9000) {
                h = hArray[i];
                for (int j = 0; j < oArray.length; j++) {
                    if (Math.abs(i-j)<3) continue;
                    if (oArray[i][0] > -9000) {
                        o = oArray[j];
                        double dx = h[0] - o[0];
                        double dy = h[1] - o[1];
                        double dz = h[2] - o[2];
                        double r2 = dx * dx + dy * dy + dz * dz;
                        double hb = 1 / (1 + Math.exp(2 * (r2 - 5)));
                        if (hb > 0.01) {
                            Number[] triplet = {new Integer(i+4), new Integer(j+4), new Double(hb)};
                            hBonds.add(triplet);
                        }
                    }
                }
            }
        }
        MeshiWriter writer = new MeshiWriter(file.getAbsolutePath().substring(0,file.getAbsolutePath().length()-4)+".hBonds.json");
        writer.print(jsonObject);
        writer.close();
    }
}
