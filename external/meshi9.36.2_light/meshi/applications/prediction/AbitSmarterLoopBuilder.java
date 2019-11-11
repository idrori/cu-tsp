/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction;

import meshi.applications.prediction.beautify.BeautifyAttribute;
import meshi.geometry.Coordinates;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.MeshiProgram;

import java.util.*;

public class AbitSmarterLoopBuilder {
    private int maxNumberOfClashes;
    public static final double CA_CA_DISTANCE = 3.8;
    public static final double CA_CA_DISTANCE2 = 3.8 * 3.8;
    public static final int MAX_CLASHES = 100;
    private static Loop nTerm;
    private static Loop cTerm;
    private static final Random random = MeshiProgram.randomNumberGenerator();

    public static void addLoops(Protein protein, double clashDistance,
                                int optimalNumberOfClashes, int nTrys) {
        nTerm = new Loop(clashDistance, optimalNumberOfClashes, nTrys);
        cTerm = new Loop(clashDistance, optimalNumberOfClashes, nTrys);
        ArrayList<Loop> loops = getLoops(protein, clashDistance,
                optimalNumberOfClashes, nTrys);
        for (Loop loop : loops) {
            loop.build(protein);
        }
        nTerm.build(protein);
        cTerm.build(protein);
        AtomList temp = protein.atoms();
    }


    private static ArrayList<Atom> getNeighbors(Atom atom, Protein protein, double clashDistance) {
        ArrayList<Atom> out = new ArrayList<Atom>();
        for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
            Atom other = (Atom) atoms.next();
            if ((!atom.nowhere()) && (!BeautifyAttribute.isVisible(atom.residue())) && (atom != other))
                if (atom.distanceFrom(other) < clashDistance * 3)
                    out.add(other);
        }
        return out;
    }

    //------------------------------------------------------------- get loops ---------------------------------------------------------
    private static ArrayList<Loop> getLoops(Protein model, double clashDistance,
                                            int optimalNumberOfClashes, int nTrys) {
        ArrayList<Loop> out = new ArrayList<Loop>();
        Chain chain = model.chain();
        int chainSize = chain.size();
        int firstReliable, lastReliable;
        Residue residue, prevResidue = null;
        //
        // Find the first and last residues that may serve as relaible anchors for loop building.
        boolean done = false;
        nTerm.add(new Residue(new ResidueIdentifier(0)));
        for (firstReliable = 0; (firstReliable < chainSize) & (!done); firstReliable++) {
            residue = chain.get(firstReliable);
            if ((!residue.dummy()) && (!BeautifyAttribute.isProblematicOrLoop(residue))) {
                done = true;
            }
            if (!residue.dummy()) {
                System.out.println("Adding " + residue + " to nTerm");
                nTerm.add(residue);
            }
        }
        if (firstReliable >= chainSize) {
            System.out.println("No relaible residues in " + model);
            return null;
        }
        firstReliable--;
        //
        done = false;
        cTerm.add(new Residue(new ResidueIdentifier(chainSize)));
        for (lastReliable = chainSize - 1; (lastReliable > firstReliable) & (!done); lastReliable--) {
            residue = chain.get(lastReliable);
            if ((!residue.dummy()) & (!BeautifyAttribute.isProblematicOrLoop(residue))) {
                done = true;
            }
            System.out.println("Adding " + residue + " to cTerm");
            cTerm.add(residue);
        }
        if (lastReliable <= firstReliable) {
            System.out.println("Too few relaible residues in " + model);
            return null;
        }
        lastReliable++;
        System.out.println("First and last reliable residues" + chain.get(firstReliable) + " " + chain.get(lastReliable));
        //
        //
        Loop loop = null;
        for (int iResidue = firstReliable; iResidue <= lastReliable; iResidue++) {
            residue = chain.get(iResidue);
            if (prevResidue != null) { // no loop without anchor
                if (loop == null) {
                    if (!BeautifyAttribute.isProblematicOrLoop(prevResidue)) { // PrevResidue may serve as anchor
                        if (BeautifyAttribute.isProblematicOrLoop(residue)) {   // Loop is needed
                            System.out.println("new loop starting at " + prevResidue);
                            loop = new Loop(clashDistance, optimalNumberOfClashes, nTrys);
                            System.out.println("adding " + prevResidue + " to the loop");
                            loop.add(prevResidue);
                            System.out.println("adding " + residue + " to the loop");
                            loop.add(residue);
                        } else {
                        } // no need for loop
                    } else {
                    } // PrevResidue may not serve as anchor. Thus we cannot create loop even if needed.
                } else {
                    if (BeautifyAttribute.isProblematicOrLoop(residue)) {
                        System.out.println("adding " + residue + " to the loop.");
                        loop.add(residue);
                    } else {
                        System.out.println("Closing loop with " + residue);
                        loop.add(residue);
                        out.add(loop);
                        loop = null;
                    }
                }
            }
            prevResidue = residue;
        }
        Object[] tempArray = out.toArray();
        Arrays.sort(tempArray, new LoopComperator());
        out = new ArrayList<Loop>();
        for (Object o : tempArray)
            out.add((Loop) o);
        return out;
    }


    private static boolean hasNowhereAtoms(Chain chain) {
        for (Iterator residues = chain.iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) {
                if (residue.ca().nowhere()) return true;
            }
        }
        return false;
    }

    //------------------------------------------------------------ Loop -------------------------------------------------------------
    private static class LoopComperator implements Comparator {
        public int compare(Object a, Object b) {
            Loop l1 = (Loop) a;
            Loop l2 = (Loop) b;
            if (l1.rank() > l2.rank()) return 1;
            if (l1.rank() < l2.rank()) return -1;
            return 0;
        }

        public boolean equals(Object a, Object b) {
            return (compare(a, b) == 0);
        }
    }

    private static class Loop extends ResidueList {
        private enum End {
            HEAD, TAIL, TERM
        }

        private double clashDistance;
        private int optimalNumberOfClashes, nTrys;
        private int rank = -1;

        public Loop(double clashDistance,
                    int optimalNumberOfClashes, int nTrys) {
            super();
            this.clashDistance = clashDistance;
            this.optimalNumberOfClashes = optimalNumberOfClashes;
            this.nTrys = nTrys;
        }

        public void build(Protein protein) {
            int first = 0;
            int last = size() - 1;
            if (get(0).dummy()) {// terminal
                System.out.println("Building terminal " + last);
                build(protein, first, last, End.TERM);
            } else {
                if (random.nextDouble() > 0.5) {
                    System.out.println("Building loop " + first + " " + last + " head");
                    build(protein, first, last, End.HEAD);
                } else {
                    System.out.println("Building loop " + first + " " + last + " tail");
                    build(protein, first, last, End.TAIL);
                }
            }

        }

        public void build(Protein protein, int first, int last, End flag) {
            boolean assignOK;
            System.out.println("Building loop " + first + " " + last + " " + flag);
            if ((first == last) || (first + 1 == last)) return;
            if (flag == End.HEAD) {
                assignOK = assignCoordinates(protein, first, first + 1, last);
                if (assignOK) {
                    if (random.nextDouble() > 0.8)
                        build(protein, first + 1, last, End.HEAD);
                    else
                        build(protein, first + 1, last, End.TAIL);
                } else build(protein, first, last, End.TAIL);
            }
            if (flag == End.TAIL) {
                assignOK = assignCoordinates(protein, last, last - 1, first);
                if (assignOK) {
                    if (random.nextDouble() > 0.8)
                        build(protein, first, last - 1, End.TAIL);
                    else
                        build(protein, first, last - 1, End.HEAD);
                } else build(protein, first, last, End.HEAD);
            }
            if (flag == End.TERM) {
                assignCoordinates(protein, last, last - 1, first);
                build(protein, first, last - 1, End.TERM);
            }
        }

        public boolean assignCoordinates(Protein protein, int fromI, int assignMeI, int toI) {
            double distance1, distance2, grade;

            System.out.println("assignCoordinates " + fromI + " " + assignMeI + " " + toI);
            Residue assignMeResidue = get(assignMeI);
            if (assignMeResidue.dummy())
                throw new RuntimeException("Cannot assign coordinates to dummy " + assignMeResidue);
            Atom from = (get(fromI)).ca();
            if (from.nowhere()) throw new RuntimeException("This is weird 10");
            Atom to = null;
            if (!(get(toI)).dummy()) {
                to = (get(toI)).ca();
                if (to.nowhere()) throw new RuntimeException("This is weird 11");
            }
            Atom assignMe = assignMeResidue.ca();
            double targetDistance1;
            if ((assignMeI + 1 == toI) || (assignMeI - 1 == toI)) {
                targetDistance1 = 3.8;
            } else targetDistance1 = 4.5;
            Coordinates targetCoordinates1 = null;
            if (to != null)
                targetCoordinates1 = new Coordinates(to.x(), to.y(), to.z());
            double targetDistance2;
            Coordinates targetCoordinates2;
            if (BeautifyAttribute.isLoop(assignMeResidue)) {
                targetCoordinates2 = targetCoordinates1;
                targetDistance2 = targetDistance1;
            } else {
                targetCoordinates2 = new Coordinates(assignMe.x(), assignMe.y(), assignMe.z());
                targetDistance2 = 0;
            }

            Coordinates[] best = new Coordinates[MAX_CLASHES];
            double[] dbest = new double[MAX_CLASHES];
            for (int i = 0; i < MAX_CLASHES; i++) {
                dbest[i] = 1000;
                best[i] = null;
            }

            ArrayList<Atom> neighbors = getNeighbors(from, protein, clashDistance);

            for (int i = 0; i < nTrys; i++) {
                grade = 0;
                Coordinates temp = new Coordinates((new Coordinates(from)), CA_CA_DISTANCE);
                assignMe.setXYZ(temp);
                int nClashes = getClashes(assignMe, neighbors, clashDistance);
                if (nClashes >= MAX_CLASHES) throw new RuntimeException("This is a black hole.");
                if (to != null)
                    distance1 = (new Coordinates(assignMe)).distanceFrom(targetCoordinates1);
                else distance1 = targetDistance1;
                //if (distance1 < 0.5) grade -= 1000;
                distance2 = (new Coordinates(assignMe)).distanceFrom(targetCoordinates2);
                //if ((targetDistance2 != 0) & (distance2 < 0.5)) grade -= 1000;
                grade += Math.sqrt((distance1 - targetDistance1) * (distance1 - targetDistance1));
                grade += 10 * Math.sqrt((distance2 - targetDistance2) * (distance2 - targetDistance2));
                if (dbest[nClashes] > grade) {
                    dbest[nClashes] = grade;
                    best[nClashes] = temp;
                }
            }


            Coordinates chosen = null;
            int i;
            for (i = optimalNumberOfClashes; (i >= 0) & (chosen == null); i--) {
                if (best[i] != null) chosen = best[i];
            }
            if (chosen == null)
                for (i = optimalNumberOfClashes + 1; (i < MAX_CLASHES) & (chosen == null); i++) {
                    if (best[i] != null) chosen = best[i];
                }
            if (chosen == null) throw new RuntimeException("Very weird");
            assignMe.setXYZ(chosen);
            System.out.println("Coordinates were assigned to " + assignMeResidue + " " + i);
            BeautifyAttribute.setVisible(assignMeResidue);
            return true;
        }

        private static int getClashes(Atom atom, ArrayList<Atom> neighbors, double clashDistance) {
            int out = 0;
            for (Atom neighbor : neighbors) {
                if (atom.distanceFrom(neighbor) < clashDistance) out++;
                if (atom.distanceFrom(neighbor) < 2) out += 20;
            }
            return out;
        }

        public int rank() {
            if (rank == -1) {
                Residue first = get(0);
                Residue last = get(size() - 1);
                rank = last.number() - first.number() + random.nextInt(3);
            }
            return rank;
        }
    }
}
