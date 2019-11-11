/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;
import meshi.sequences.AlignmentException;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.info.*;

import java.util.ArrayList;

/**

 */
public class ModelAnalyzer {
    public final Protein nativeStructure, model, originalModel;
    private TotalEnergy energy;
    private ResidueAlignmentMethod residueAlignmentMethod;
    private enum ContactsMode {CA,CB}

    public void setEnergy(TotalEnergy energy) {
        this.energy = energy;
    }

    public ModelAnalyzer(Protein model,
                         Protein nativeStructure,
                         Protein originalModel,
                         TotalEnergy energy,
                         ResidueAlignmentMethod residueAlignmentMethod) {
        this.nativeStructure = nativeStructure;
        this.model = model;
        this.originalModel = originalModel;
        this.energy = energy;
        this.residueAlignmentMethod = residueAlignmentMethod;
    }

    public double rms() throws Exception{
        return Rms.rms(nativeStructure, model, residueAlignmentMethod);

    }
    public double change() throws Exception{
        return Rms.rms(originalModel, model, residueAlignmentMethod);

    }

    public double rmsHeavy() throws Exception{
        return Rms.rmsHeavy(nativeStructure, model, residueAlignmentMethod);

    }


    public ProteinInfoOLd analyze(String comment) throws AlignmentException {
        return analyze(comment, null);
    }

    public ProteinInfoOLd analyze(String comment, ChainsInfo chainsInfo) throws AlignmentException{

            ProteinInfoOLd out = new ProteinInfoOLd(comment, model.sourceFile(), model.name(), model);
        if (nativeStructure != null) {
            double[] gdt_ts;
            gdt_ts = Rms.gdt(nativeStructure, model);
            out.add(new DoubleInfoElement(InfoType.GDT_TS, "GDT_TS form native structure", gdt_ts[0]));
            out.add(new DoubleInfoElement(InfoType.GDT_TS1, "GDT_TS form native structure", gdt_ts[1]));
            out.add(new DoubleInfoElement(InfoType.GDT_TS2, "GDT_TS form native structure", gdt_ts[2]));
            out.add(new DoubleInfoElement(InfoType.GDT_TS4, "GDT_TS form native structure", gdt_ts[3]));
            out.add(new DoubleInfoElement(InfoType.GDT_TS8, "GDT_TS form native structure", gdt_ts[4]));
            double[] originalGdt = Rms.gdt(nativeStructure, originalModel);
            out.add(new DoubleInfoElement(InfoType.DELTA_GDT_TS, "delta GDT_TS with respect to the original unrefined", gdt_ts[0] - originalGdt[0]));
            double[] gdt_ha = Rms.gdt(nativeStructure, model,GDTcalculator.Type.HA);
            out.add(new DoubleInfoElement(InfoType.GDT_HA, "GDT_HA form native structure", gdt_ha[0]));
            double tmScore = (new TMscore(model, nativeStructure)).findMaxScore();
            out.add(new DoubleInfoElement(InfoType.TM_SCORE, "TM-score", tmScore));
            double [] contacts = getContacts(nativeStructure, model, ContactsMode.CA, chainsInfo);
            out.add(new DoubleInfoElement(InfoType.CONTACTS6_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[0]));
            out.add(new DoubleInfoElement(InfoType.CONTACTS8_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[1]));
            out.add(new DoubleInfoElement(InfoType.CONTACTS10_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[2]));
            contacts = getContacts(nativeStructure, model, ContactsMode.CB);
            out.add(new DoubleInfoElement(InfoType.CB_CONTACTS6_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[0]));
            out.add(new DoubleInfoElement(InfoType.CB_CONTACTS8_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[1]));
            out.add(new DoubleInfoElement(InfoType.CB_CONTACTS10_MCC, "Matthews correlation coefficient of the predicted and native contact maps ", contacts[2]));
            try {
                double[] originalGdtHa = Rms.gdt(nativeStructure, originalModel,GDTcalculator.Type.HA);
                out.add(new DoubleInfoElement(InfoType.DELTA_GDT_HA, "delta GDT_HA with respect to the original unrefined", gdt_ha[0] - originalGdtHa[0]));
                double rms = Rms.rms(nativeStructure, model, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.RMS, "RMS form native structure", rms));
                double originalRms = Rms.rms(nativeStructure, originalModel, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.DELTA_RMS, "delta RMS with respect to the original unrefined model", rms - originalRms));
                double rmsHeavy = Rms.rmsHeavy(nativeStructure, model, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.RMS_HEAVY, "Heavy getAtoms RMS form native structure", rmsHeavy));
                double originalRmsHeavy = Rms.rmsHeavy(nativeStructure, originalModel, residueAlignmentMethod);
                out.add(new DoubleInfoElement(InfoType.DELTA_RMS_HEAVY, "delta RMS of heavy getAtoms with respect to the original unrefined model", rmsHeavy - originalRmsHeavy));
            }catch (Exception ex ) {
                ex.printStackTrace();
                System.out.println();
                throw new RuntimeException("__________________________________\n"+ex.getMessage());
            }
            double rmsFromOriginal;
            try {
                rmsFromOriginal = Rms.rms(originalModel, model, residueAlignmentMethod);
        out.add(new DoubleInfoElement(InfoType.CHANGE, "structural change (in RMS) due to refinement", rmsFromOriginal));
            } catch (Exception ex ) {throw new RuntimeException(ex.getMessage());}
        }

        if (energy != null) {
            double e = energy.evaluateAll(chainsInfo);
            out.add(new DoubleInfoElement(InfoType.ENERGY, "energy", e));
            for (AbstractEnergy energyElement : energy.energyTerms()) {
                out.add(new DoubleInfoElement(energyElement.info()));
                for (MeshiInfo element : energyElement.info().getChildren()) {
                    if ((element != null) && (element.getValue()!= null))
                        out.add(new DoubleInfoElement(element));
                }
            }
        }

        if (chainsInfo != null)
            for (ChainInfo chainInfo : chainsInfo) {
                for (ResidueInfo residueInfo : chainInfo)
                    out.add(residueInfo);
            }
        return out;
    }

    private static double[] getContacts(Protein nativeStructure, Protein model, ContactsMode contactsMode) {
        return getContacts(nativeStructure,model,contactsMode, null);
    }
    private static double[] getContacts(Protein nativeStructure, Protein model, ContactsMode contactsMode, ChainsInfo chainsInfo) {
        double[] out = {getContacts(nativeStructure, model, 6, contactsMode), getContacts(nativeStructure, model, 8, contactsMode, chainsInfo) ,
                        getContacts(nativeStructure, model, 10, contactsMode)};
        return out;
    }

    private static Atom getAtom(Residue residue, ContactsMode contactsMode) {
        if ((residue.type == ResidueType.GLY) | (contactsMode == ContactsMode.CA)) return residue.ca();
        else return residue.cb();
    }

    public static double getContacts(Protein nativeStructure, Protein model, double threshold, ContactsMode contactsMode) {
        return getContacts(nativeStructure, model, threshold, contactsMode, null);
    }
    public static double getContacts(Protein nativeStructure, Protein model, double threshold, ContactsMode contactsMode, ChainsInfo chainsInfo) {
        Atom atomInative, atomImodel, atomJnative, atomJmodel;
        ResidueAlignment residueAlignment = null;
        try {
            ChainList chains0 = nativeStructure.chains();
            ChainList chains1 = model.chains();
            String name0 = nativeStructure.name();
            String name1 = model.name();
            residueAlignment = new ResidueAlignment(chains0, name0, chains1, name1, ResidueAlignmentMethod.IDENTITY);//modified
        } catch (AlignmentException ex) {
            Utils.throwException("Static function Rms.gdt", ex, "Failed to align " + nativeStructure + " and " + model);
        }

        ChainList chains = model.chains();
        ArrayList<int[]> tp = new ArrayList();
        ArrayList<int[]> tn = new ArrayList();
        ArrayList<int[]> fp = new ArrayList();
        ArrayList<int[]> fn = new ArrayList();
        ArrayList<int[]> nativeContacts = new ArrayList();
        ArrayList<int[]> modelContacts  = new ArrayList();
        for (Chain chain : chains) {
            int[] array1 = new int[chain.size()];
            int[] array2 = new int[chain.size()];
            int[] array3 = new int[chain.size()];
            int[] array4 = new int[chain.size()];
            int[] array5 = new int[chain.size()];
            int[] array6 = new int[chain.size()];
            for (int i = 0; i < chain.size(); i++)
                array1[i] = array2[i] = array3[i] = array4[i] = -1;  //that is not visited
            tp.add(array1);
            tn.add(array2);
            fp.add(array3);
            fn.add(array4);
            nativeContacts.add(array5);
            modelContacts.add(array6);

        }
        for (int iResidue = 0; iResidue < residueAlignment.size(); iResidue++) {
            ResidueAlignmentColumn columnI = residueAlignment.get(iResidue);
            Residue modelResidueI  = columnI.residue1();
            int chainNumber = modelResidueI.getChainNumber();
            int number = modelResidueI.number();
            tp.get(chainNumber)[number] = fp.get(chainNumber)[number]  = tn.get(chainNumber)[number]  = fn.get(chainNumber)[number] = 0;
        }
        for (int iResidue = 0; iResidue < residueAlignment.size(); iResidue++) {
            ResidueAlignmentColumn columnI = residueAlignment.get(iResidue);
            Residue nativeResidueI = columnI.residue0();
            Residue modelResidueI  = columnI.residue1();
            if (nativeResidueI.type != modelResidueI.type())
                throw new RuntimeException("Weird ResidueAlignmentColumn # "+ iResidue +" "+ columnI+"\n"+residueAlignment);
            atomInative = getAtom(nativeResidueI,contactsMode);
            atomImodel  = getAtom(modelResidueI,contactsMode);
            int modelChainI = modelResidueI.getChainNumber();
            int modelResidueInunber = modelResidueI.number();
            for (int jResidue = iResidue + 1; jResidue < residueAlignment.size(); jResidue++) {
                ResidueAlignmentColumn columnJ = residueAlignment.get(jResidue);
                Residue nativeResidueJ = columnJ.residue0();
                Residue modelResidueJ  = columnJ.residue1();
                if (nativeResidueJ.type != modelResidueJ.type())
                    throw new RuntimeException("Weird ResidueAlignmentColumn ## "+ jResidue +" "+ columnJ+"\n"+residueAlignment);
                if (nativeResidueJ.number() - nativeResidueI.number() < 5) continue;
                atomJnative = getAtom(nativeResidueJ,contactsMode);
                atomJmodel  = getAtom(modelResidueJ,contactsMode);
                int modelChainJ = modelResidueJ.getChainNumber();
                int modelResidueJnunber = modelResidueJ.number();
                double nativeDistance = atomInative.distanceFrom(atomJnative);
                double modelDistance  = atomImodel.distanceFrom(atomJmodel);
                boolean nativeContact = nativeDistance <= threshold;
                boolean modelContact  = modelDistance <= threshold;
                if (nativeContact) {
                    if (modelContact) {tp.get(modelChainI)[modelResidueInunber]++;   tp.get(modelChainJ)[modelResidueJnunber]++; }
                    else {fn.get(modelChainI)[modelResidueInunber]++;   fn.get(modelChainJ)[modelResidueJnunber]++;};
                } else {
                    if (modelContact) {fp.get(modelChainI)[modelResidueInunber]++;   fp.get(modelChainJ)[modelResidueJnunber]++; }
                    else {tn.get(modelChainI)[modelResidueInunber]++;   tn.get(modelChainJ)[modelResidueJnunber]++; };
                }
            }
        }

        int nativeContactsS = getContacts(nativeStructure,contactsMode, threshold, null);
        int modelContactsS  = getContacts(model,contactsMode, threshold, modelContacts);


        int tpS = 0, tnS = 0, fpS = 0, fnS = 0;
        for (int i = 0; i < tp.size(); i++) {
            int length = tp.get(i).length;
            ChainInfo chainInfo = null;
            if (chainsInfo != null) chainInfo = chainsInfo.get(i);
            for (int j = 0; j < length; j++) {
                if ((tp.get(i)[j] == -1) & (chainsInfo != null)) {
                    chainInfo.get(j).add(new DoubleInfoElement(InfoType.RESIDUE_CONTACTS8_MCC,
                            "Per-residue MCC of native and model contatacts", -2.0)); // A meaningless number for correlation
                    chainInfo.get(j).add(new DoubleInfoElement(InfoType.CONTACTS8, "contacts8", modelContacts.get(i)[j]));
                }
                else {
                    if (!chains.get(i).get(j).dummy()) {
                        tpS += tp.get(i)[j];
                        tnS += tn.get(i)[j];
                        fpS += fp.get(i)[j];
                        fnS += fn.get(i)[j];
                        if (chainInfo != null) {
                            double mcc;
                            mcc = matthewsCorrelationCoefficient(tp.get(i)[j], fn.get(i)[j], fp.get(i)[j], tn.get(i)[j]);
                            chainInfo.get(j).add(new DoubleInfoElement(InfoType.RESIDUE_CONTACTS8_MCC,
                                    "Per-residue MCC of native and model contatacts", mcc));
                            chainInfo.get(j).add(new DoubleInfoElement(InfoType.CONTACTS8, "contacts8", tp.get(i)[j] + fp.get(i)[j]));
                        }
                    }
                }
            }
        }
        // A bit redundant I know but clearer
        int nativeContactsOutOfTheAlignment = nativeContactsS - tpS - fnS;
        fnS = fnS + nativeContactsOutOfTheAlignment;
        int modelContactsOutOfTheAlignment = modelContactsS -tpS - fpS;
        fpS = fpS + modelContactsOutOfTheAlignment;
        return matthewsCorrelationCoefficient(tpS, fnS, fpS, tnS);
    }

    public static int getContacts(Protein protein, ContactsMode contactsMode, double threshold, ArrayList<int[]> contacts) {
        int contactsSum = 0;
        for (int iResidue = 0; iResidue < protein.residues().size(); iResidue++) {
            Residue residueI = protein.residues().get(iResidue);
            if (!residueI.dummy()) {
                Atom atomI = getAtom(residueI, contactsMode);
                for (int jResidue = iResidue + 5; jResidue < protein.residues().size(); jResidue++) {
                    Residue residueJ = protein.residues().get(jResidue);
                    if (!residueJ.dummy()) {
                        Atom atomJ = getAtom(residueJ, contactsMode);
                        double distance = atomI.distanceFrom(atomJ);
                        boolean contact = distance <= threshold;
                        if (contact) {
                            contactsSum += 2;
                            if (contacts != null) {
                                contacts.get(atomI.residue().getChainNumber())[atomI.residueNumber()]++;
                                contacts.get(atomJ.residue().getChainNumber())[atomJ.residueNumber()]++;
                            }
                        }
                    }
                }
            }
        }
        return contactsSum;
    }

    public static double matthewsCorrelationCoefficient(double tp, double fn, double fp, double tn)  {
        double sum = tp+fn+fp+tn;
        if (sum == 0) return 1;
        if (((tp + fp) == 0) & (fn > 0))  return 0;
        if ((tp+fp)*(tp+fn)*(tn+fn)*(tn + fn) == 0) return 0;
        double mcc = (tp*tn - fp*fn)/Math.sqrt((tp+fp)*(tp+fn)*(tn+fn)*(tn + fn));
        if (mcc < 0) return 0;
        return mcc;
    }

}

