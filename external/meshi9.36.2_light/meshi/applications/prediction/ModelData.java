/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.*
        ;
import meshi.util.*;

import java.util.*;

public class ModelData extends ArrayList<Pair> {
    private static final int NAME_INDEX = 0;
    private static final int AVG_ENERGY_INDEX = 1;
    private static final int RMS_INDEX = 2;
    private static final int GDT_INDEX = 3;
    private static final int BASE_INDEX = 4;

    public ModelData() {
        super();
    }

    public ModelData(Protein model, Protein reference,
                     TotalEnergy totalEnergy) throws UpdateableException,EvaluationException,AlignmentException {
        this(model, reference, totalEnergy, model.name());
    }

    public ModelData(Protein model, Protein reference, TotalEnergy totalEnergy,
                     String fileName) throws UpdateableException, EvaluationException,AlignmentException {
        super();
        add(new Pair("name", fileName));

        double energy = totalEnergy.evaluate();
        int numberOfAtoms = Utils.numberOfAtomsWithCoordinates(model);
        add(new Pair("avgEnergy", energy / numberOfAtoms));

        double rms = -1, gdt = -1;
        if (reference != null) {
            ResidueAlignment modelRefeenceAlignment = new ResidueAlignment(model.chain(), model.name(), reference.chain(), reference.name(),
                    ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            rms = Rms.rms(modelRefeenceAlignment, Rms.RmsType.CA);
            gdt = (Rms.gdt(modelRefeenceAlignment, Utils.numberOfNonDummyResidues((reference))))[0];
        }
        add(new Pair("RMS", rms));
        add(new Pair("GDT_TS", gdt));
        add(new Pair("zRMS", rms));
        add(new Pair("zGDT", gdt));

        for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
            totalEnergy.off();
            AbstractEnergy term = (AbstractEnergy) terms.next();
            term.on();
            energy = totalEnergy.evaluate() / numberOfAtoms;
            totalEnergy.evaluateAtoms();
            double weightedEnergy = getWeightedEnergy(model, numberOfAtoms);
            double coolAtomsEnergy = getCoolAtomsEnergy(model, numberOfAtoms);
            add(new Pair(term.comment(), energy));
            add(new Pair("w" + term.comment(), weightedEnergy));
            add(new Pair("ca" + term.comment(), coolAtomsEnergy));
        }
    }

    public ModelData(ModelData source) {
        this();
        for (Pair pair : source) {
            if (pair.value() instanceof String) {
                add(new Pair(pair.key, pair.value()));
            } else {
                add(new Pair(pair.key, 0));
            }
        }
    }

    public static ModelData processProtein(Protein protein, Protein reference,
                                           CommandList commands,
                                           EnergyCreator[] energyCreators) throws UpdateableException , EvaluationException,AlignmentException{
        DistanceMatrix distanceMatrix = protein.atoms().molecularSystem().getDistanceMatrix();
        TotalEnergy totalEnergy = new TotalEnergy(protein, energyCreators, commands,"Generic total energy");
        ModelData out = new ModelData(protein, reference, totalEnergy);
        return out;
    }

    private static double getWeightedEnergy(Protein model, int numberOfAtoms) {
        double cmX = 0, cmY = 0, cmZ = 0, norm = 0;
        double dx, dy, dz, d2;
        double energy, weight, wEnergy;
        Atom atom;

        for (Iterator atoms = model.atoms().iterator(); atoms.hasNext();) {
            atom = (Atom) atoms.next();
            if (!atom.nowhere()) {
                energy = atom.energy() / numberOfAtoms;
                weight = Math.exp(-energy / 100);
                cmX += atom.x() * weight;
                cmY += atom.y() * weight;
                cmZ += atom.z() * weight;
                norm += weight;
            }
        }
        cmX /= norm;
        cmY /= norm;
        cmZ /= norm;

        wEnergy = 0;
        norm = 0;
        for (Iterator atoms = model.atoms().iterator(); atoms.hasNext();) {
            atom = (Atom) atoms.next();
            if (!atom.nowhere()) {
                energy = atom.energy();
                dx = atom.x() - cmX;
                dy = atom.y() - cmY;
                dz = atom.z() - cmZ;
                d2 = dx * dx + dy * dy + dz * dz;
                weight = Math.exp(-d2 / 100);
                wEnergy += energy * weight;
                norm += weight;
            }
        }
        return wEnergy / norm;
    }

    private static double getCoolAtomsEnergy(Protein model, int numberOfAtoms) {
        Atom atom;
        double energy;
        double[] energies = new double[model.atoms().size()];
        int size = 0;
        for (Iterator atoms = model.atoms().iterator(); atoms.hasNext();) {
            atom = (Atom) atoms.next();
            if (!atom.nowhere()) size++;
        }

        int i = 0;
        for (Iterator atoms = model.atoms().iterator(); atoms.hasNext();) {
            atom = (Atom) atoms.next();
            if (!atom.nowhere())
                energies[i] = atom.energy() / numberOfAtoms;
            else energies[i] = 10000000;
            i++;
        }

        Arrays.sort(energies);
        double totalE = 0;
        for (i = 0; i < numberOfAtoms * 9 / 10; i++)
            totalE += energies[i];
        return totalE;
    }


    public String toString(boolean headerFlag, boolean nameFlag, boolean zScoresFlag, String lable) {
        String out = "";
        if (headerFlag) {
            out += "H_" + lable + "\t";
            for (Pair pair : this)
                out += String.format("%-20s\t", pair.key);
            out += "\n";
        }
        if (nameFlag & zScoresFlag) {
            out += "V_" + lable + "\t";
            for (Pair pair : this)
                if (pair.value() instanceof String)
                    out += String.format("%-20s\t", pair.value());
                else
                    out += String.format("%-20.4f\t", pair.value());
        }
        if (nameFlag & (!zScoresFlag)) {
            out += lable + "\t";
            out += String.format("%-20s\t", get(0).value());
        }
        if ((!nameFlag) & zScoresFlag) {
            out += lable + "\t";
            for (Pair pair : this)
                if (pair.value() instanceof Double)
                    out += String.format("%-20.4f\t", pair.value());
        }

        return out;
    }

    public String toString() {
        return toString(true, true, true, "ModelData");
    }

    public static ModelData sumModel(List<ModelData> models) {
        ModelData sum = new ModelData(models.get(0));
        int size = sum.size();
        for (ModelData model : models) {
            for (int i = 1; i < size; i++) { // the 0th element is the protein name
                Pair pair = model.get(i);
                Pair pairSum = sum.get(i);
                double value = pair.doubleValue();
                double sumVal = pairSum.doubleValue();
                pairSum.setValue(sumVal + value);
            }
        }
        return sum;
    }

    public static ModelData sum2Model(List<ModelData> models) {
        ModelData sum2 = new ModelData(models.get(0));
        int size = sum2.size();
        for (ModelData model : models) {
            for (int i = 1; i < size; i++) { // the 0th element is the protein name
                Pair pair = model.get(i);
                Pair pairSum = sum2.get(i);
                double value = pair.doubleValue();
                double sumVal = pairSum.doubleValue();
                pairSum.setValue(sumVal + value * value);
            }
        }
        return sum2;
    }

    public static ModelData avgModel(List<ModelData> models, ModelData sum) {
        ModelData avg = new ModelData(models.get(0));
        int size = sum.size();
        int numberOfModels = models.size();
        for (int i = 1; i < size; i++) { // the 0th element is the protein name
            Pair pairSum = sum.get(i);
            Pair pairAvg = avg.get(i);
            double sumVal = pairSum.doubleValue();
            double avgVal = sumVal / numberOfModels;
            pairAvg.setValue(avgVal);
        }
        return avg;
    }

    public static ModelData stdModel(ModelData avg, ModelData avg2) {
        ModelData std = new ModelData(avg);
        int size = avg.size();
        for (int i = 1; i < size; i++) { // the 0th element is the protein name
            Pair pairAvg = avg.get(i);
            Pair pairAvg2 = avg2.get(i);
            Pair pairStd = std.get(i);
            double avgVal = pairAvg.doubleValue();
            double avg2Val = pairAvg2.doubleValue();
            pairStd.setValue(Math.sqrt(avg2Val - avgVal * avgVal));
        }
        return std;
    }


    public static List<ModelData> zScoreOfModels(List<ModelData> models) {
        ModelData sumModel = sumModel(models);
        ModelData sum2Model = sum2Model(models);
        ModelData avgModel = avgModel(models, sumModel);
        ModelData avg2Model = avgModel(models, sum2Model);
        ModelData stdModel = stdModel(avgModel, avg2Model);

        int size = sumModel.size();
        List<ModelData> out = new ArrayList<ModelData>();
        for (ModelData modelData : models) {
            ModelData zScores = new ModelData(modelData);
            for (int i = 0; i < BASE_INDEX; i++)
                zScores.get(i).setValue(modelData.get(i).value());
            for (int i = BASE_INDEX; i < size; i++) {
                if (stdModel.get(i).doubleValue() == 0) zScores.get(i).setValue(0);
                else {
                    double zScore = (modelData.get(i).doubleValue() - avgModel.get(i).doubleValue()) / stdModel.get(i).doubleValue();
                    zScores.get(i).setValue(zScore);
                }
            }
            out.add(zScores);
        }
        return out;
    }

    public static ModelData[] toArray(List<ModelData> list) {
        ModelData[] out = new ModelData[list.size()];
        int i = 0;
        for (ModelData md : list) {
            out[i] = md;
            i++;
        }
        return out;
    }

    public static Comparator<ModelData> termComparator(String key) {
        return new TermComparator(key);
    }

    private static class TermComparator implements Comparator<ModelData> {
        private String key;

        public TermComparator(String key) {
            this.key = key;
        }

        public int compare(ModelData md1, ModelData md2) {
            Pair p1 = null, p2 = null;
            for (Pair p : md1)
                if (p.key.equals(key)) p1 = p;
            for (Pair p : md2)
                if (p.key.equals(key)) p2 = p;
            if ((p1 == null) || (p2 == null)) {
                System.out.println("********************* " + key);
                for (Pair p : md1) System.out.println("md1 " + p);
                for (Pair p : md1) System.out.println("md1 " + p);
                throw new RuntimeException("This is weird ");
            }
            double v1 = p1.doubleValue();
            double v2 = p2.doubleValue();
            if (v1 > v2) return 1;
            if (v1 < v2) return -1;
            return 0;
        }
    }
}
