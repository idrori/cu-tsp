/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.*;

import java.util.*;

import meshi.util.*;

public class Grid {
    private GridStatus gridStatus;
    public final double edge;
    private final int minSphereRadiusToEdge;
    private final long MAX_GRID_SIZE = 120000;
    private final double SPHERE_COND = 1.3;
    private static int largestGridSize = -999;
    int xSize = 0, ySize = 0, zSize = 0;
    int xSizeNew, ySizeNew, zSizeNew;
    public double minX, minY, minZ, maxX, maxY, maxZ;
    MolecularSystem atoms;
    GridCell[][][] cells;
    GridCell defaultCell;
    double[][] prevCoor;
    ArrayList<AtomCore> movingAtoms;
    int cellCapacity;


    public Grid(double edge, double minSphereRadius) {
        this.edge = edge;
        //setResidue number of the neighboring cells.
        int temp = (int) (minSphereRadius / edge);
        if (temp == (minSphereRadius / edge))
            minSphereRadiusToEdge = temp;
        else
            minSphereRadiusToEdge = temp + 1;
    }

    public Grid(MolecularSystem atoms, double edge, double minSphereRadius) throws UpdateableException {
        if (atoms.size() == 0)
            throw new RuntimeException("No atoms for grid");
        this.atoms = atoms;
        this.edge = edge;
        //setResidue number of the neighboring cells.
        int temp = (int) (minSphereRadius / edge);
        if (temp == (minSphereRadius / edge))
            minSphereRadiusToEdge = temp;
        else
            minSphereRadiusToEdge = temp + 1;
        double edge3 = edge * edge * edge;
        //
        prevCoor = new double[atoms.size()][3];
        for (int i = 0; i < prevCoor.length; i++) {
            prevCoor[i][0] = prevCoor[i][1] = prevCoor[i][2] = -9999.9999;
            movingAtoms = new ArrayList();
            if (edge <= 2) cellCapacity = 2;
            else if (edge > 500) cellCapacity = atoms.size();
            else cellCapacity = round(edge3 / 8);
        }
        gridStatus = build();
        if (gridStatus != GridStatus.OK) throw new RuntimeException(" Failed to build the grid due to " + gridStatus.comment());
        for (int i = 0; i < prevCoor.length; i++) {
            prevCoor[i][0] = prevCoor[i][1] = prevCoor[i][2] = -9999.9999;
            movingAtoms = new ArrayList();
            if (edge <= 2) cellCapacity = 2;
            else if (edge > 500) cellCapacity = atoms.size();
            else cellCapacity = round(edge3 / 8);
        }
    }


    /**
     * Builds the grid. Return false upon failure (say, if the getAtoms are too spread around
     * in a space that is too large to be devided into cells without a memory failure.
     */

    public GridStatus build() throws UpdateableException {
        GridStatus gridStatus;
        //int size = getAtoms.size();

        // First round over all getAtoms.
        // updates min/max values and identify getAtoms that have moved 0.3 * buffer size
        updateMinMaxAndMovingAtoms();
                // end of first round
        xSizeNew = round((maxX - minX) / edge) + 1;
        ySizeNew = round((maxY - minY) / edge) + 1;
        zSizeNew = round((maxZ - minZ) / edge) + 1;

        if ((xSize != xSizeNew) || (ySize != ySizeNew) || (zSize != zSizeNew)) {
            gridStatus = checkNewGridSize();
            if (gridStatus != GridStatus.OK)
                return gridStatus;
            xSize = xSizeNew;
            ySize = ySizeNew;
            zSize = zSizeNew;
            cells = null;
            Runtime.getRuntime().gc();
            cells = new GridCell[xSize][ySize][zSize];
            checkCells();
            for (int X = 0; X < xSize; X++)
                for (int Y = 0; Y < ySize; Y++)
                    for (int Z = 0; Z < zSize; Z++)
                        try {
                            cells[X][Y][Z] = new GridCell(cellCapacity);
                        }
                        catch (Exception ex) {
                            System.out.println(ex);
                            System.out.println("Very weird exception " + X + " " + Y + " " + Z + " " + cells.length +
                                    " " + cells[0].length + " " + cells[0][0].length);
                            throw new UpdateableException();
                        }

        } else { // Grid size did not change.
            for (int X = 0; X < xSize; X++)
                for (int Y = 0; Y < ySize; Y++)
                    for (int Z = 0; Z < zSize; Z++) {
                        try {
                            cells[X][Y][Z].clear();
                        } catch (RuntimeException ex) {
                            System.out.println("cells.length = " + cells.length + " ; xSize = " + xSize + " ; X = " + X);
                            System.out.println("cells[X].length = " + cells[X].length + " ; ySize = " + ySize + " ; Y = " + Y);
                            System.out.println("cells[X][Y].length = " + cells[X][Y].length + " ; zSize = " + zSize + " ; Z = " + Z);
                            for (Object o : atoms)
                                System.out.println(o);
                            throw ex;
                        }
                    }
        }

        for (AtomCore atom : movingAtoms) {
            int X = round((atom.x() - minX) / edge);
            int Y = round((atom.y() - minY) / edge);
            int Z = round((atom.z() - minZ) / edge);

            int xFrom = ((X < minSphereRadiusToEdge) ? 0 : (X - minSphereRadiusToEdge));
            int xTo = ((X >= xSize - minSphereRadiusToEdge) ? xSize : (X + minSphereRadiusToEdge + 1));
            int yFrom = ((Y < minSphereRadiusToEdge) ? 0 : (Y - minSphereRadiusToEdge));
            int yTo = ((Y >= ySize - minSphereRadiusToEdge) ? ySize : (Y + minSphereRadiusToEdge + 1));
            int zFrom = ((Z < minSphereRadiusToEdge) ? 0 : (Z - minSphereRadiusToEdge));
            int zTo = ((Z >= zSize - minSphereRadiusToEdge) ? zSize : (Z + minSphereRadiusToEdge + 1));

            if (minSphereRadiusToEdge <= 1) {
                for (int i = xFrom; i < xTo; i++)
                    for (int j = yFrom; j < yTo; j++)
                        for (int k = zFrom; k < zTo; k++)
                            cells[i][j][k].add(atom);
            } else {
                double cellR = (minSphereRadiusToEdge + 0.5) * (minSphereRadiusToEdge + 0.5);
                for (int i = xFrom; i < xTo; i++)
                    for (int j = yFrom; j < yTo; j++) {
                        if (((i - X) * (i - X) + (j - Y) * (j - Y)) / cellR > SPHERE_COND) continue;
                        for (int k = zFrom; k < zTo; k++) {
                            if (((i - X) * (i - X) + (k - Z) * (k - Z)) / cellR > SPHERE_COND) continue;
                            cells[i][j][k].add(atom);
                        }
                    }
            }
        }
        return GridStatus.OK;
    }

    public GridCell getCell(AtomCore atom) {
        int X = round((atom.x() - minX) / edge);
        int Y = round((atom.y() - minY) / edge);
        int Z = round((atom.z() - minZ) / edge);
        if ((X<0)|(Y<0)|(Z<0))
            throw new RuntimeException("Weird grid cell "+X+" "+Y+" "+Z+"\n"+atom);
        if ((X >= cells.length)||(Y >= cells[0].length) || (Z >= cells[0][0].length))
            throw new RuntimeException("This is weird. \n"+atom.atom + " "+atom.status()+
                                       "\n minX minY minZ = "+minX+" "+minY+" "+minZ+
                                       "\ncells size = "+cells.length+" "+cells[0].length+" "+cells[0][0].length);
        return cells[X][Y][Z];
    }

    public static int round(double d) {
        return (int) (d + 0.5);
    }

    public GridStatus status() {
        return gridStatus;
    }

    private void updateMinMaxAndMovingAtoms(){
        minX = minY = minZ = 10000;
        maxX = maxY = maxZ = -10000;
        double bufferOneThirdSqr = DistanceMatrix.bufferOneThirdSqr;
        movingAtoms.clear();
        for (AtomCore atom : atoms) {
            if (atom.status().activeOrImage()) {
                double x = atom.x();
                double y = atom.y();
                double z = atom.z();
                if (x < minX) minX = x;
                if (y < minY) minY = y;
                if (z < minZ) minZ = z;
                if (x > maxX) maxX = x;
                if (y > maxY) maxY = y;
                if (z > maxZ) maxZ = z;
                int atomNumber = atom.number;
                double dx = x - prevCoor[atomNumber][0];
                double dy = y - prevCoor[atomNumber][1];
                double dz = z - prevCoor[atomNumber][2];
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 > bufferOneThirdSqr) {
                    prevCoor[atomNumber][0] = x;
                    prevCoor[atomNumber][1] = y;
                    prevCoor[atomNumber][2] = z;
                    movingAtoms.add(atom);
                    /// debug ///
//                    if ((atomNumber == 828) || (atomNumber == 16))
//                        System.out.println("Atom " + atom + " has moved more than " + Math.sqrt(bufferOneThirdSqr) + " " + edge);
                    //// debug ////
                }
            }
        }
    }

    private GridStatus checkNewGridSize() {
        if (xSize > 0) { // that is, it is not the first step and very large changes in the grid sizes are very likely to cause problems
            if (xSizeNew > xSize + 10)
                return GridStatus.FAILED_AXIS_CHANGED_TOO_FAST.setComment("A weird event during  grid creation - xSize has changed from " + xSize + " to " + xSizeNew+"\n"+
                        minX+" "+maxX);
            if (ySizeNew > ySize + 10)
                return GridStatus.FAILED_AXIS_CHANGED_TOO_FAST.setComment("A weird event during  grid creation - ySize has changed from " + ySize + " to " + ySizeNew);
            if (zSizeNew > zSize + 10)
                return GridStatus.FAILED_AXIS_CHANGED_TOO_FAST.setComment("A weird event during  grid creation - zSize has changed from " + zSize + " to " + zSizeNew);
        }
        double newVolume = xSizeNew * ySizeNew * zSizeNew;
        // the test below checks whether any of these numbers is NaN
        if (Double.isNaN(newVolume))
            throw new RuntimeException("Weird dimenssions " + xSize + "  " + ySize + "   " + zSize);
        if (newVolume > MAX_GRID_SIZE) {
            for (int i = 0; i < atoms.size(); i++) System.out.println(atoms.get(i));
            return GridStatus.FAILED_GRID_TOO_LARGE.setComment(" An Error wile creating a new grid\n" +
                    "The requested grid size " + xSizeNew + "*" + ySizeNew + "*" + zSizeNew + "=" +
                    (newVolume) +
                    "is larger than MAX_GRID_SIZE=" + MAX_GRID_SIZE + "\n" +
                    minX + " " + minY + " " + minZ + " " + maxX + " " + maxY + " " + maxZ + " ");

        }
        if ((xSizeNew < 0) || (ySizeNew < 0) || (zSizeNew < 0)) {
            System.out.println("************************************************************\n" + " An Error wile creating a new grid");
            for (AtomCore atom : atoms) System.out.println("#### " + atom);
            throw new RuntimeException("Weird grid dimention " + xSize + " " + ySize + " " + zSize + "\n" +
                    minX + " " + minY + " " + minZ + " " + maxX + " " + maxY + " " + maxZ + " ");
        }
        if (newVolume > largestGridSize) {
            largestGridSize = (int) newVolume;
            Utils.println("largest grid size " + largestGridSize);
        }
        return GridStatus.OK;
    }

    private void checkCells()throws UpdateableException{
        if ((cells.length != xSize) |
            (cells[0].length != ySize) |
            (cells[0][0].length != zSize)) {
            System.out.println(" An Error wile creating a new grid" +
                    "The program failed in building the requested "
                    + xSize + "*" + ySize + "*" + zSize + " array.\n" +
                    "probably the array is too large. You may solve this problem by " +
                    "reducing MAX_GRID_SIZE which is currently " + MAX_GRID_SIZE);
            throw new UpdateableException();
        }
    }

}
