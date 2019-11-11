package meshi.util;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.util.info.ChainsInfo;

/**
 * Created by chen on 24/08/2017.
 */
public class ResidueData  {
    double[][] chainsData;
    public ResidueData(ChainsInfo chainsInfo) {
        chainsData = new double[chainsInfo.size()][];
        for (int i = 0; i < chainsData.length; i++) {
            chainsData[i] = new double[chainsInfo.get(i).size()];
        }
    }

    public void add(Atom atom, double value) {
        int residueNumber = atom.residueNumber();
        int chainNumber = atom.residue().getChainNumber();
        chainsData[chainNumber][residueNumber] += value;
    }

    public double get(int chainNumber, int residueNumber) {
        return chainsData[chainNumber][residueNumber];
    }
}