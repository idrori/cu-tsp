/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**

 */
public class ProteinInfoListOld extends ArrayList<ProteinInfoOLd> {
    public final String name;

    public ProteinInfoListOld(String name) {
        this.name = name;
    }

    public void print() throws IOException {
        MeshiWriter writer = new MeshiWriter(name + ".pil");
        print(writer);
        writer.close();
    }

    public void filterByOriginalModel() {
        ProteinInfoOLd infoI, infoJ;
        for (int i = 0; i < size(); i++) {
            infoI = get(i);
            if (infoI.on()) {
                for (int j = i + 1; j < size(); j++) {
                    infoJ = get(j);
                    if (infoJ.on()) {
                        if (infoI.originalModel().equals(infoJ.originalModel()))
                            infoJ.turnOff();
                    }
                }
            }
        }
    }

    public void sort() {
        ProteinInfoOLd temp;
        double sumI, sumJ;
        ProteinInfoOLd[] sorted = new ProteinInfoOLd[size()];
        ProteinInfoOLd infoI, infoJ;
        for (int i = 0; i < size(); i++) {
            sorted[i] = get(i);
        }

        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < i; j++) {
                if (sorted[i].sum() > sorted[j].sum()) {
                    temp = sorted[i];
                    sorted[i] = sorted[j];
                    sorted[j] = temp;
                }
            }
        }
        clear();
        for (ProteinInfoOLd info : sorted)
            add(info);
    }

    public void print(MeshiWriter writer) throws IOException {
        if (size() == 0) {
            System.out.println("Empty list " + this);
            return;
        }
        String header = get(0).header();
        writer.println(header);
        for (ProteinInfoOLd proteinInfoOLd : this) {
            if (proteinInfoOLd.on()) {
                if (!proteinInfoOLd.header().equals(header)) {
                    header = proteinInfoOLd.header();
                    writer.println("%" + header);
                }
                writer.println(proteinInfoOLd.values());
                writer.flush();
            }
        }
    }

    public void print(MeshiInfoXMLwriter writer) throws IOException {
        print(writer, 0);
    }

    public void print(MeshiInfoXMLwriter writer, int indentationTabs) throws IOException {
        if (size() == 0) {
            System.out.println("Empty list " + this);
            return;
        }
        writer.println("<ProteinInfoList name=\"" + name + "\">");
        for (ProteinInfoOLd proteinInfoOLd : this) {
            if (proteinInfoOLd.on()) {
                writer.println(proteinInfoOLd.toXML(indentationTabs + 1));
            }
        }
        writer.println("</ProteinInfoList>");
        writer.flush();
    }

    public ProteinInfoListOld toZscore() {
        ProteinInfoOLd sum = new ProteinInfoOLd(this.get(0), "sum", null, "sum of protein info");
        sum.toZero();
        ProteinInfoOLd sum2 = new ProteinInfoOLd(this.get(0), "sum2", null, "sum of protein info squared");
        sum2.toZero();
        for (ProteinInfoOLd proteinInfoOLd : this) {
            sum.addElementValue(proteinInfoOLd);
            sum2.addElementValue2(proteinInfoOLd);
        }

        ProteinInfoOLd avg = new ProteinInfoOLd(sum, "avg", sum.fileName(), "averaged protein info");
        avg.divideBy(size());

        ProteinInfoOLd avg2 = new ProteinInfoOLd(sum2, "avg2", sum2.fileName(), "averaged protein info squared");
        avg2.divideBy(size());

        ProteinInfoOLd avgavg = new ProteinInfoOLd(avg, "avgavg", avg.fileName(), "squared average of protein info");
        avgavg.square();

        ProteinInfoOLd variance = new ProteinInfoOLd(avg2, "variance", avg2.fileName(), "Variance of protein info.");
        variance.subtract(avgavg);

        ProteinInfoOLd std = new ProteinInfoOLd(variance, "std", variance.fileName(), "Standard deviation of protein info.");
        std.sqrt();

        ProteinInfoListOld zScores = new ProteinInfoListOld(name + ".zScores");
        for (ProteinInfoOLd proteinInfoOLd : this) {
            ProteinInfoOLd zScore = new ProteinInfoOLd(proteinInfoOLd, "temp1", proteinInfoOLd.fileName(), "temp2");
            zScore.subtract(avg);
            zScore.divideBy(std);
            zScores.add(zScore);
        }
        return zScores;
    }
}
