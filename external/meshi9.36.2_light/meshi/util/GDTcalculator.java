
/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.applications.hubbard.HubbardPosition;
import meshi.applications.prediction.homology.CaAtomFilter;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.sequences.AtomAlignment;
import meshi.sequences.ResidueAlignment;
import meshi.util.file.MeshiWriter;


public class GDTcalculator {
    public static double[][] C1;
    public static double[][] C2;
    public static HubbardPosition[] Popo;
    public static int len;
    public static final int SUB = 4; //the length of the first sub protein that we check.
    //According to Hubbard plot - suppose to be 3 but there is a problem
    //to the Hubbard overlap to handle it.
    public static final double maximal_distance = 4.9;//the maximal distance between 2 Ca atoms for joining them to the
    //group of Atoms we send at the next iteration.
    public static final int num_of_loops = 20;//the maximal number of iterations over the protein in finding
    //the best superposition.
    public enum Type {
        TS(1.0, 2.0, 4.0, 8.0),
        HA(0.5, 1.0, 2.0, 4.0);
        private final double[] thresholds = new double[4];
        private Type(double t0, double t1, double t2, double t3) {
              thresholds[0] = t0;
              thresholds[1] = t1;
              thresholds[2] = t2;
              thresholds[3] = t3;
        }
    }


    public static double[] gdt(ResidueAlignment residueAlignment, int refLength) {
        return gdt(residueAlignment, refLength, Type.TS);
    }
    public static double[] gdt(ResidueAlignment residueAlignment, int refLength, Type type) {
        return gdt(residueAlignment, refLength, null,type);
    }

    public static double[] gdt(ResidueAlignment residueAlignment, int refLength, MeshiWriter output) {
        return gdt(residueAlignment,refLength, output, Type.TS);
    }
    public static double[] gdt(ResidueAlignment residueAlignment, int refLength, MeshiWriter output, Type type) {
          AtomAlignment atomAlignment = new AtomAlignment(ResidueAlignment.toAlignmentArray(residueAlignment),new CaAtomFilter());
        AtomList list0 = atomAlignment.atomList(0);
        AtomList list1 = atomAlignment.atomList(1);
        len = 0;
        for (int i = 0; i < list1.size(); i++) {
            Atom atom0 = list0.get(i);
            Atom atom1 = list1.get(i);
            if (atom0.residueName().equals(atom1.residueName()))
                len++;
            else {
                if (output == null) System.out.println("Warrning:\n" + atom0 + "\n" + atom1);
                else output.println("Warrning:\n" + atom0 + "\n" + atom1);
            }
        }
        C1 = new double[3][len];
        C2 = new double[3][len];
        int iAtom = 0;
        for (int i = 0; i < list1.size(); i++) {
            Atom atom0 = list0.get(i);
            Atom atom1 = list1.get(i);
            if (atom0.residueName().equals(atom1.residueName())) {
                C1[0][iAtom] = atom0.x();
                C1[1][iAtom] = atom0.y();
                C1[2][iAtom] = atom0.z();
                C2[0][iAtom] = atom1.x();
                C2[1][iAtom] = atom1.y();
                C2[2][iAtom] = atom1.z();
                iAtom++;
            }
        }
        return gdt(list0, list1, refLength,type);
    }

    public static double[] gdt(AtomList reference, AtomList protein1, int refLength) {
        return gdt(reference, protein1,refLength,Type.TS);
    }
    public static double[] gdt(AtomList reference, AtomList protein1, int refLength, Type type) {
        //creating two arrays (one for each given structure)
        //according to the atoms coordinates writen at the files.
        init(reference, protein1);

        double[] out = new double[5];
        Popo = new HubbardPosition[len + 1 - SUB];//represents all the superspositions of the protein.
        initialize(Popo);//initialize Popo (for each subunit of the 
        //protein at the base subunit size.
        loop_over_Popo(Popo);//for each start point (according to the base subunit)
        //finding the best superposition of the protein.
        double[] best_distance = new double[len];
        build_best_distance(best_distance, Popo);//building the best distances array.

        out[1] = fracBelowRMS(best_distance, type.thresholds[0], refLength);
        out[2] = fracBelowRMS(best_distance, type.thresholds[1], refLength);
        out[3] = fracBelowRMS(best_distance, type.thresholds[2], refLength);
        out[4] = fracBelowRMS(best_distance, type.thresholds[3], refLength);
        out[0] = 0.25 * (out[1] + out[2] + out[3] + out[4]);

        return out;

    }//gdt


    //---------------------------------------------------------------------------------------

    public static void initialize(HubbardPosition[] Popo) {
        int place;
        for (place = 0; place < Popo.length; place++)   //initialize the basic states
        {
            Popo[place] = new HubbardPosition(SUB, place, maximal_distance, C1, C2);
        }
    }//End of initializing


    //--------------------------------------------------------------------------------------
    public static void loop_over_Popo(HubbardPosition[] Popo) {
        for (int i = 0; i < Popo.length; i++)//go over all our HubbardPositions.
        {
            Popo[i].find_best_conformation(num_of_loops);
        }//each position now has an array of the smallest sum of distances
    }


    //---------------------------------------------------------------------------------------

    public static void build_best_distance(double[] best_distance, HubbardPosition[] Popo) {
        int i, j;
        for (i = 0; i < len; i++) {// building the best_distance array
            double temp_best = 100000000;
            for (j = 0; j < Popo.length; j++) {
                if (Popo[j].numbers[i] < temp_best) {
                    temp_best = Popo[j].numbers[i];
                }
            }
            best_distance[i] = Math.sqrt(temp_best);
        }//the best distance's array is ready
    }//build_best_distance


    //-----------------------------------------------------------------------------------------------
    // find the fraction of the protein overlap (length of len) that is below a certain rms
    // isOn the Hubbard plot curve.
    public static double fracBelowRMS(double[] best_distance, double rms, int referLength) {
        for (int c = 0; c < len; c++)
            if (best_distance[c] > rms) {
                return (c* 1.0) / referLength;
            }
        return (len * 1.0) / referLength;
    }

    //-----------------------------------------------------------------------------------------------

    public static void init(AtomList al1, AtomList al2) {
        if (al1.size() != al2.size()) throw new RuntimeException("The atom lists need to have the same size");
        al1 = al1.CAFilter();
        al2 = al2.CAFilter();
        if (al1.size() != al2.size()) throw new RuntimeException("The two CA atom lists need to have the same size");
        for (int i = 0; i < al1.size(); i++) {
            if (!al1.get(i).name().equals(al2.get(i).name())) {
                throw new RuntimeException("Atom " + al1.get(i) + "  does not fit\n" + al2.get(i));
            }
/*
            if (al1.get(i).residueNumber() != al2.get(i).residueNumber()) {
                throw new RuntimeException("residue number of Atom " + al1.get(i) + "  does not fit\n" + al2.get(i));
            }
*/
        }
        len = al1.size();
        C1 = new double[3][al1.size()];
        C2 = new double[3][al1.size()];
        for (int i = 0; i < al1.size(); i++) {
            Atom atom1 = al1.get(i);
            Atom atom2 = al2.get(i);
            if (!atom1.type().backboneCA())
                throw new RuntimeException("GDT_TS is applicable only to C-alpha atoms 1: " + atom1);
            if (!atom2.type().backboneCA())
                throw new RuntimeException("GDT_TS is applicable only to C-alpha atoms 2: " + atom2);
            C1[0][i] = atom1.x();
            C1[1][i] = atom1.y();
            C1[2][i] = atom1.z();
            C2[0][i] = atom2.x();
            C2[1][i] = atom2.y();
            C2[2][i] = atom2.z();
        }
    }//init


}//class 
    
    
	



        	
   
