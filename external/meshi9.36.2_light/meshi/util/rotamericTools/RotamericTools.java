/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.rotamericTools;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.TorsionList;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.ResidueType;
import meshi.util.KeyWords;
import meshi.util.overlap.Overlap;


public class RotamericTools implements KeyWords, Rot1Arrays {


    // This method puts all the residues of prot, that have at least ONE NON-FROZEN atom, into their rot1 position.
    // It also returns an Nx4 double array where row n is the phi,psi,type,{rotamer 1 prob} of residue n.
    // Note that the distance matrix is discarded after this method is run, as it is the best way to ensure,
    // that it is updated properly after this sidechain change.
    // In version meshi.6.22 (7.3.2010) chen added the boolean flag "conservative",which if setResidue "true" indicate that the side chain should not be moved.
    public static double[][] putIntoRot1(Protein prot, DistanceMatrix dm, DunbrackLib lib) {
        return putIntoRot1(prot, dm, lib);
    }

    public static double[][] putIntoRot1(Protein protein, DistanceMatrix dm, DunbrackLib lib, boolean conservative) {

        double[][] pp = phipsi(protein, dm);
        double[] tmprot;

        int maxRes = -1;
        Chain chain = protein.chain();
        for (int kk = 0; kk < protein.atoms().size(); kk++)
            if (protein.atoms().atomAt(kk).residueNumber() > maxRes)
                maxRes = protein.atoms().atomAt(kk).residueNumber();
        double[][] largepp = new double[maxRes + 1][];
        for (int kk = 0; kk < largepp.length; kk++)
            if ((chain.residueAt(kk) != null) && (chain.residueAt(kk).type().ordinal() > -1))
                if ((kk < pp.length) && (pp[kk] != null))
                    largepp[kk] = pp[kk];
                else {
                    largepp[kk] = new double[3];
                    largepp[kk][0] = -60 * Math.PI / 180.0;
                    largepp[kk][1] = -40 * Math.PI / 180.0;
                    largepp[kk][2] = chain.residueAt(kk).type().ordinal();
                }
        pp = largepp;

        // putting Rot1 into all instances of this residue type
        if (!conservative) {
            for (int kk = 0; kk < pp.length; kk++) {
                if ((chain.residueAt(kk) != null) && (!chain.residueAt(kk).dummy())) {
                    if ((pp[kk] != null) && (pp[kk][2] > 0) && (pp[kk][2] != 5) &&
                            (chain.residueAt(kk).getAtoms().defrostedAtoms().size() > 0)) {
                        tmprot = lib.getRotamer((int) pp[kk][2], pp[kk][0], pp[kk][1], 0);
                        ResidueBuilder.build(chain.residueAt(kk),
                                chain.residueAt(kk).type(),
                                tmprot);
                    }
                }
            }
        }

        double[][] tmppp = new double[pp.length][];
        for (int iResidue = 0; iResidue < pp.length; iResidue++) {
            if ((chain.residueAt(iResidue) != null) && (!chain.residueAt(iResidue).dummy())) {
                if ((pp[iResidue] != null) && (pp[iResidue][2] > 0) && (pp[iResidue][2] != 5) &&
                        (chain.residueAt(iResidue).getAtoms().defrostedAtoms().size() > 0)) {
                    tmppp[iResidue] = new double[4];
                    tmppp[iResidue][0] = pp[iResidue][0];
                    tmppp[iResidue][1] = pp[iResidue][1];
                    tmppp[iResidue][2] = pp[iResidue][2];
                    tmppp[iResidue][3] = lib.getRotamerProb((int) pp[iResidue][2], pp[iResidue][0], pp[iResidue][1], 0);
                }
            }
        }
        // Discarding the distance matrix.
        dm.terminator.kill("The distance matrix is not updated correctly in PutIntoRot1. It has to be created anew.");
        return tmppp;
    }  // of putIntoRot1


    // Returns the an Nx3 double array where row n is the phi,psi,type of residue n.
    public static double[][] phipsi(Protein protein, DistanceMatrix dm) {
        int maxResNum = -999;
        int resNum;
        double[][] pp;

        TorsionList phiList =  TorsionList.createTorsionList(protein,
                dm).filter(new TorsionList.FilterPhi());
        TorsionList psiList =  TorsionList.createTorsionList(protein,
                dm).filter(new TorsionList.FilterPsi());

        for (int i = 0; i < phiList.size(); i++)
            if (phiList.get(i).getTorsionResNum() > maxResNum)
                maxResNum = phiList.get(i).getTorsionResNum();

        for (int i = 0; i < psiList.size(); i++)
            if (psiList.get(i).getTorsionResNum() > maxResNum)
                maxResNum = psiList.get(i).getTorsionResNum();

        pp = new double[maxResNum + 1][];

        for (int i = 0; i < phiList.size(); i++) {
            resNum = phiList.get(i).getTorsionResNum();
            if (pp[resNum] == null) {
                pp[resNum] = new double[3];
                pp[resNum][0] = -60 * Math.PI / 180.0;
                pp[resNum][1] = -40 * Math.PI / 180.0;
                pp[resNum][2] = -1;
            }
            pp[resNum][0] = phiList.get(i).torsion();
            pp[resNum][2] = ResidueType.type(phiList.get(i).getTorsionResName()).ordinal();
        }
        for (int i = 0; i < psiList.size(); i++) {
            resNum = psiList.get(i).getTorsionResNum();
            if (pp[resNum] == null) {
                pp[resNum] = new double[3];
                pp[resNum][0] = -60 * Math.PI / 180.0;
                pp[resNum][1] = -40 * Math.PI / 180.0;
                pp[resNum][2] = -1;
            }
            pp[resNum][1] = psiList.get(i).torsion();
            pp[resNum][2] = ResidueType.type(psiList.get(i).getTorsionResName()).ordinal();
        }

        return pp;
    } // of phipsi


    // This method finds the most similar (in term of rms) rotamer of residue 'res'.
    // If it returns an array than a similar rotamer was found in the Dunbrack library, and it is returned. 
    // If it return null than an error occured, probably some atoms are missing in the side chain. 
    /*public static double[] nearestRot(DunbrackLib lib, Residue res, double phi, double psi) {
        AtomList al1 = res.atoms();
        AtomList al2 = Utils.duplicateInAnewMolecularSystem(al1);
        ResidueType restype = res.type();
        int numOfRotamers = lib.getRotamerNum(restype.ordinal(), phi, psi);
        double[] bestrot = new double[lib.getChiMax(restype.ordinal())];
        ;
        double[] tmprot;
        double rms, minrms = 999999999.9;
        boolean cont = true;
        for (int ind = 0; cont && (ind < numOfRotamers); ind++) {
            tmprot = lib.getRotamer(restype.ordinal(), phi, psi, ind);
            try {
                ResidueBuilder.build(res, restype, tmprot);
            }
            catch (Exception e) {
                cont = false;
            }
            rms = calcRMS(al1, al2);
            if (rms < 0.0)
                cont = false;
            else if (rms < minrms) {
                minrms = rms;
                for (int tmp = 0; tmp < tmprot.length; tmp++)
                    bestrot[tmp] = tmprot[tmp];
            }
            //             if ((restype==2) || (restype==3) || (restype==4) || (restype==14) || (restype==19)) {
            //                 if ((restype==2) || (restype==4) || (restype==19))
            //                     tmprot[1] += Math.PI;
            //                 if (restype==3)
            //                     tmprot[2] += Math.PI;
            //                 if (restype==14)
            //                     tmprot[3] += Math.PI;
            if ((restype == ResidueType.ASP) || (restype == ResidueType.GLU) || (restype == ResidueType.PHE) || (restype == ResidueType.ARG) || (restype == ResidueType.TYR)) {
                if ((restype == ResidueType.ASP) || (restype == ResidueType.GLU) || (restype == ResidueType.TYR))
                    tmprot[1] += Math.PI;
                if (restype == ResidueType.GLU)
                    tmprot[2] += Math.PI;
                if (restype == ResidueType.ARG)
                    tmprot[3] += Math.PI;
                try {
                    ResidueBuilder.build(res, restype, tmprot);
                }
                catch (Exception e) {
                    cont = false;
                }
                rms = calcRMS(al1, al2);
                if (rms < 0.0)
                    cont = false;
                else if (rms < minrms) {
                    minrms = rms;
                    for (int tmp = 0; tmp < tmprot.length; tmp++)
                        bestrot[tmp] = tmprot[tmp];
                }
            }
        }
        // returning the positions in residue to the original
        for (int c = 0; c < al2.size(); c++) {
            if (!al1.atomAt(c).name.equals(al2.atomAt(c).name)) {
                throw new RuntimeException("The lists are not ordered the same. This should not have happened!!");
            }
            al1.atomAt(c).setXYZ(al2.atomAt(c).x(), al2.atomAt(c).y(), al2.atomAt(c).z());
        }
        if (cont)
            return bestrot;
        else
            return null;
    }
*/

    // al1,al2 - are two atom lists from two instances a certain residue. 
    // This method overlap the backbone atoms of the two residues and returns the RMS of their
    // (non-hydrogen) sidechain atoms (including the CB atom).
    // If the two lists do not contain exactly the same atoms (in term of names not instances) then -1 is returned.
    // The order in the lists is irrelevent.
    public static double calcRMS(AtomList al1, AtomList al2) {
        al1 = al1.filter(new AtomList.NonHydrogen()).noOXTFilter();
        al2 = al2.filter(new AtomList.NonHydrogen()).noOXTFilter();
        if (al1.size() != al2.size())
            return -1;
        boolean found;
        double rms = 0.0;
        int backboneCount = 0;
        int anum = 0;
        int arCount = 5;
        double[][] co = new double[3][al1.size()];
        double[][] co2 = new double[3][al1.size()];
        int[] aux = {0, 1, 2, 3, 4};
        for (int c1 = 0; c1 < al1.size(); c1++) {
            found = false;
            for (int c2 = 0; c2 < al2.size(); c2++)
                if (al1.atomAt(c1).name.equals(al2.atomAt(c2).name)) {
                    found = true;
                    if (al1.atomAt(c1).name.equals("N")) {
                        co[0][0] = al1.atomAt(c1).x();
                        co[1][0] = al1.atomAt(c1).y();
                        co[2][0] = al1.atomAt(c1).z();
                        co2[0][0] = al2.atomAt(c2).x();
                        co2[1][0] = al2.atomAt(c2).y();
                        co2[2][0] = al2.atomAt(c2).z();
                        backboneCount++;
                    } else if (al1.atomAt(c1).name.equals("CA")) {
                        co[0][1] = al1.atomAt(c1).x();
                        co[1][1] = al1.atomAt(c1).y();
                        co[2][1] = al1.atomAt(c1).z();
                        co2[0][1] = al2.atomAt(c2).x();
                        co2[1][1] = al2.atomAt(c2).y();
                        co2[2][1] = al2.atomAt(c2).z();
                        backboneCount++;
                    } else if (al1.atomAt(c1).name.equals("C")) {
                        co[0][2] = al1.atomAt(c1).x();
                        co[1][2] = al1.atomAt(c1).y();
                        co[2][2] = al1.atomAt(c1).z();
                        co2[0][2] = al2.atomAt(c2).x();
                        co2[1][2] = al2.atomAt(c2).y();
                        co2[2][2] = al2.atomAt(c2).z();
                        backboneCount++;
                    } else if (al1.atomAt(c1).name.equals("O")) {
                        co[0][3] = al1.atomAt(c1).x();
                        co[1][3] = al1.atomAt(c1).y();
                        co[2][3] = al1.atomAt(c1).z();
                        co2[0][3] = al2.atomAt(c2).x();
                        co2[1][3] = al2.atomAt(c2).y();
                        co2[2][3] = al2.atomAt(c2).z();
                        backboneCount++;
                    } else if (al1.atomAt(c1).name.equals("CB")) {
                        anum++;
                        co[0][4] = al1.atomAt(c1).x();
                        co[1][4] = al1.atomAt(c1).y();
                        co[2][4] = al1.atomAt(c1).z();
                        co2[0][4] = al2.atomAt(c2).x();
                        co2[1][4] = al2.atomAt(c2).y();
                        co2[2][4] = al2.atomAt(c2).z();
                        backboneCount++;
                    } else {
                        if (arCount >= co[0].length)
                            return -1;
                        anum++;
                        co[0][arCount] = al1.atomAt(c1).x();
                        co[1][arCount] = al1.atomAt(c1).y();
                        co[2][arCount] = al1.atomAt(c1).z();
                        co2[0][arCount] = al2.atomAt(c2).x();
                        co2[1][arCount] = al2.atomAt(c2).y();
                        co2[2][arCount] = al2.atomAt(c2).z();
                        arCount++;
                    }
                }
            if (!found)
                return -1;
        }
        if (backboneCount < 5) {
            return -1;
        }

        Overlap.rmsPartial(co, co2, aux);
        for (int c1 = 4; c1 < co[0].length; c1++)
            rms = rms + (co[0][c1] - co2[0][c1]) * (co[0][c1] - co2[0][c1]) +
                    (co[1][c1] - co2[1][c1]) * (co[1][c1] - co2[1][c1]) +
                    (co[2][c1] - co2[2][c1]) * (co[2][c1] - co2[2][c1]);

        rms = rms / anum;
        return rms;
    }

    // This methods assign every sidechain torsion a random angle.
    // All the atoms must be present in the protein. 
    public static void jumble(Protein prot) {
        double[] tmprot = {-999, -999, -999, -999};
        for (int kk = 0; kk < prot.residues().size(); kk++) {
            tmprot[0] = 2 * Math.PI * (0.5 - Math.random());
            tmprot[1] = 2 * Math.PI * (0.5 - Math.random());
            tmprot[2] = 2 * Math.PI * (0.5 - Math.random());
            tmprot[3] = 2 * Math.PI * (0.5 - Math.random());
            if (prot.residue(kk).type() != ResidueType.DMY)
                ResidueBuilder.build(prot.residue(kk), prot.residue(kk).type(), tmprot);
        }
    }


    public static double[] getMean(int type, double prob, int what) {
        double[] result = new double[2];
        int best = -999, secondbest = -999;
        for (int c = 0; c < 10; c++)
            if (mean[type][c][what] > -0.999) {
                if (Math.abs(c / 10.0 + 0.05 - prob) <= Math.abs(best / 10.0 + 0.05 - prob)) {
                    secondbest = best;
                    best = c;
                } else if (Math.abs(c / 10.0 + 0.05 - prob) <= Math.abs(secondbest / 10.0 + 0.05 - prob)) {
                    secondbest = c;
                }
            }
        if (Math.abs(best / 10.0 + 0.05 - prob) > 0.162)
            System.out.println("\n\nWarning!!! A prob of: " + prob + " that is too far from the closest:" +
                    (best / 10.0 + 0.05) + " type: " + type + "\n\n");
        result[0] = mean[type][best][what] +
                (mean[type][secondbest][what] - mean[type][best][what]) / (secondbest / 10.0 - best / 10.0) * (prob - (best / 10.0 + 0.05));
        result[1] = std[type][best][what] +
                (std[type][secondbest][what] - std[type][best][what]) / (secondbest / 10.0 - best / 10.0) * (prob - (best / 10.0 + 0.05));

        return result;
    }

} // of class

