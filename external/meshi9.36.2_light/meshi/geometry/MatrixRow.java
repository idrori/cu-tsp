/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;
import meshi.parameters.AtomType;

import java.util.ArrayList;
import java.util.Iterator;

public class MatrixRow extends DistanceList {
    final DistanceMatrix matrix;
    final protected double rMax2;
    final double rMaxPlusBuffer2;
    final DistanceList nonBondedRow;
    final DistanceList newNonBondedRow;
    protected static long debugCounter = 0;
    public final static boolean debug = false;


    public MatrixRow(AtomCore atomOne, DistanceMatrix dm, DistanceList nonBondedRow, DistanceList newNonBondedRow) {
        super(dm.molecularSystem.size(), atomOne);
        this.matrix = dm;
        this.rMax2 = dm.rMax2;
        rMaxPlusBuffer2 = dm.rMaxPlusBuffer2;
        this.nonBondedRow = nonBondedRow;
        this.newNonBondedRow = newNonBondedRow;

    }


    public void update() {
        double x1,y1,z1;
        double dx,dy, dz;
        double d,d2;
        DistanceMode mode;
        Distance distance;
        AtomCore atom2;
        x1 = atomOne.x();
        y1 = atomOne.y();
        z1 = atomOne.z();
        nonBondedRow.clear();
        newNonBondedRow.clear();

        for (Iterator<Distance> distanceIterator = iterator(); distanceIterator.hasNext();) { // I do not use the generic foreach since I'll need the iterator to remove elements
            distance = distanceIterator.next();
            if (distance == null) continue;
            mode = distance.mode;
//            if (mode.mirror || mode.free)
//                    throw new RuntimeException("Weird distance :\n"+
//                            distance.atom1()+"\n"+distance.atom1().core+"\n"+
//                            distance.atom2()+"\n"+distance.atom2().core+"\n"+
//                            distance+"\n"+
//                            "Distance mirrors and free distances have nothing to do here\n" + distance);
            if (mode.frozen) continue;
            atom2 = distance.atom2;
            dx = distance.dx = x1 - atom2.x();
            dy = distance.dy = y1 - atom2.y();
            dz = distance.dz = z1 - atom2.z();
            d2 = dx*dx + dy*dy + dz*dz;
            if ((d2 < rMax2) || mode.bonded) {
                d = distance.distance = Math.sqrt(d2);
                distance.invDistance = 1 / d;
                if (!mode.bonded) {
                    nonBondedRow.add(distance);
                    if (mode == DistanceMode.INFINITE) { // That is a distance that was in the buffer zone and now entered rMax
                        distance.mode = DistanceMode.NEW;
                        newNonBondedRow.add(distance);
                    } else {
                        distance.mode = DistanceMode.NORMAL;
                    }
                }
            } else if (d2 < rMaxPlusBuffer2) { //That is in the buffer zone
                        distance.mode = DistanceMode.INFINITE;
            } else { //Outside the buffer zone
                  distanceIterator.remove();
                  distance.mode = DistanceMode.DEAD; // just in case somebody is looking at this object
            }
        }
        if (debug & (debugCounter%100 == 1))
            test();
    }

    public Distance search(int atom2number) {
        for (Distance distance : this) {
            if ((distance != null) &&
                    (distance.atom2Number == atom2number)) return distance;
        }
        return null;
    }


    public String toString() {
        String out = "MatrixRow atom = " + atomOne + " number = " + atomOneNumber + "\n";
        for (Distance dis : this) {
            if (dis != null) out += " " + dis.toString() + " ; ";
            else out += " null ; ";
        }
        return out;
    }

    public void addCell(GridCell cell) {
        double x = atomOne.x();
        double y = atomOne.y();
        double z = atomOne.z();
        double dx, dy, dz, d2, dis;
        double rMaxPlusBuffer2 = DistanceMatrix.rMaxPlusBuffer2;
        for (AtomCore cellAtom : cell) {
            int cellAtomNumber = cellAtom.number;
            dx = x - cellAtom.x();
            dy = y - cellAtom.y();
            dz = z - cellAtom.z();
            d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < rMaxPlusBuffer2) {
                if (atomOneNumber >cellAtomNumber) {
                    insert(cellAtom, d2, dx, dy, dz);
                }  else {
                    if (atomOneNumber <cellAtomNumber) {
                        ((MatrixRow) matrix.search(cellAtomNumber)).insert(atomOne, d2, -dx, -dy, -dz);
                    }
                }
            }
        }

    }

    private void insert(AtomCore cellAtom,
                        double d2, double dx, double dy, double dz) {
        boolean found;
        double dis;
        int    lastEmpty = -1;
        Distance distance;
        AtomCore atom2;
        found = false;
        for (int i = 0; i < size(); i++) {
            distance = get(i);
            if (distance == null) {
                lastEmpty = i;
            } else {
                atom2 = distance.atom2;
                if (cellAtom == atom2) {
                    found = true;
                    break;
                }
            }
        }
        if (!found) {// needs to be inserted
            if (d2 < rMax2) dis = Math.sqrt(d2);
            else dis = Distance.INFINITE_DISTANCE;

            if ((atomOne.status() == AtomStatus.FROZEN) & (cellAtom.status() == AtomStatus.FROZEN))
                distance = new FrozenDistance(atomOne, cellAtom, dx, dy, dz, dis);
            else distance = new Distance(atomOne, cellAtom, dx, dy, dz, dis);
            if (dis < Distance.INFINITE_DISTANCE) {
                distance.mode = DistanceMode.NEW;
                nonBondedRow.add(distance);
                newNonBondedRow.add(distance);
            }
            else    distance.mode = DistanceMode.INFINITE; // That is within the buffer zone.
            if (lastEmpty == -1) add(distance);
            else set(lastEmpty,distance);
        }
    }

    public void test() {
        MolecularSystem atomList = matrix.molecularSystem;
        double x = atomOne.x(),y = atomOne.y() ,z = atomOne.z();
        double dx , dy, dz, d2;


        for (Distance distance :this) {
            if ((distance == null) || (distance.mode() == DistanceMode.INFINITE) || distance.mode() == DistanceMode.BONDED)continue;
            Atom atom2 = distance.atom2();
            dx = x - atom2.x();
            dy = y - atom2.y();
            dz = z - atom2.z();
            d2 = dx*dx + dy * dy + dz*dz;
            if (d2 > rMax2)
                throw new RuntimeException("MatrixRow error #1 "+debugCounter+"\n"+this+"\n"+Math.sqrt(d2)+"   "+distance);
        }
        for (AtomCore atomTwo : atomList) {
            if (atomOne.status().nowhere() || atomTwo.status().nowhere()) continue;
            if (atomTwo.number >= atomOneNumber) continue;
            dx = x - atomTwo.x();
            dy = y - atomTwo.y();
            dz = z - atomTwo.z();
            d2 = dx*dx + dy * dy + dz*dz;
            if (d2 < rMax2) {
                Distance distance = search(atomTwo.number);
                if (distance == null)
                    throw new RuntimeException("MatrixRow error #2 "+debugCounter+"\n"+this+"  "+size()+"\ndistance = "+Math.sqrt(d2)+"\n"+atomOne.atom+"\n"+atomTwo.atom);
            }
        }


    }
}

