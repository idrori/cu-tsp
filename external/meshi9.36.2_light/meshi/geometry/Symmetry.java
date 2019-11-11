/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.AtomList;
import meshi.util.Rms;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 30/05/2010
 * Time: 19:05:50
 * To change this template use File | Settings | File Templates.
 */
public class Symmetry {
    public final double[][] matrix = new double[3][3];
    public final double[] translation = new double[3];

    public Symmetry(double[] translation, double[][] matrix) {
        for (int i = 0; i < 3; i++) {
            this.translation[i] = translation[i];
            for (int j = 0; j < 3; j++)
                this.matrix[i][j] = matrix[i][j];
        }
    }

    public void apply(AtomList atoms) {

    }

    public static Symmetry fromRms(Rms rms) {
        double[] translation = new double[3];
        Coordinates tempCoordinates = rms.getCenterOfMass1();
        translation[0] = tempCoordinates.x();
        translation[1] = tempCoordinates.y();
        translation[2] = tempCoordinates.z();
        return new Symmetry(translation, rms.getMatrix());

    }

    public static void printSymmetry(MeshiWriter writer, Symmetry symmetry, int number) {
        String out;
        for (int i = 1; i <= 3; i++) {
            out = String.format("REMARK 290   SMTRY%d   %d%10.6f%10.6f%10.6f%15.6f", i, number,
                    symmetry.matrix[i - 1][0], symmetry.matrix[i - 1][1], symmetry.matrix[i - 1][2], symmetry.translation[i - 1]);
            writer.println(out);
        }
    }

    public static ArrayList<Symmetry> getSymmetryList(MeshiLineReader reader) {
        ArrayList<Symmetry> out = new ArrayList<Symmetry>();
        StringList symmetryLines = new StringList(reader, SymmetryFilter.filter);
        if (symmetryLines.size() % 3 != 0)
            throw new RuntimeException("Weird number of Symmetry lines " + symmetryLines.size());
        int numberOfSymmetryElements = symmetryLines.size() / 3;
        for (int i = 0; i < numberOfSymmetryElements; i++) {
            out.add(getSymmetry(symmetryLines, i));
        }
        return out;
    }

    public static Symmetry getSymmetry(StringList symmetryLines, int number) {
        String line0 = symmetryLines.get(number * 3);
        String line1 = symmetryLines.get(number * 3 + 1);
        String line2 = symmetryLines.get(number * 3 + 2);

        String[] words0 = line0.split(" ");
        String[] words1 = line1.split(" ");
        String[] words2 = line2.split(" ");

        double[] translation = new double[3];
        translation[0] = Double.valueOf(words0[7]);
        translation[1] = Double.valueOf(words1[7]);
        translation[2] = Double.valueOf(words2[7]);

        double[][] matrix = new double[3][3];
        for (int i = 0; i < 3; i++) {
            matrix[0][i] = Double.valueOf(words0[i + 4]);
            matrix[1][i] = Double.valueOf(words1[i + 4]);
            matrix[2][i] = Double.valueOf(words2[i + 4]);
        }

        return new Symmetry(translation, matrix);
    }

    private static class SymmetryFilter implements Filter {
        public static final SymmetryFilter filter = new SymmetryFilter();

        public boolean accept(Object obj) {
            String line = (String) obj;
            return (line.indexOf("SNTRY") != -1);
        }
    }

}
