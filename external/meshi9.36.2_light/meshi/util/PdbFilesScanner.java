/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.EvaluationException;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;
import meshi.sequences.AlignmentException;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;
import meshi.util.info.ProteinInfoList;
import meshi.util.info.ProteinInfoOLd;
import meshi.util.proteinReaders.ProteinReader;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;

/**

 */
public class PdbFilesScanner {
    private ArrayList<ProteinMetaData>  filesData = new ArrayList();

    public PdbFilesScanner(File directory, FileFilter filter, String target, CommandList commands) {
        if (!directory.exists()) throw new RuntimeException(directory + " does not a exist.");
        if (!directory.isDirectory()) throw new RuntimeException(directory + " is not a directory.");
        File[] files = directory.listFiles(filter);
        System.out.println("Scanning " + files.length + " files from " + directory.getAbsolutePath());
        for (File file : files) {
            if (target.equals("decoys"))
                target = file.getName().substring(0,8);
            ProteinMetaData proteinMetaData =  new ProteinMetaData(directory.getAbsolutePath(), file.getName(), target);
            filesData.add(proteinMetaData);
        }
    }

    public PdbFilesScanner(File directory, FileFilter filter, String target) {
        this(directory, filter, target, null);
    }

    public ProteinInfoList analyze(ProteinAnalyzer[] analyzers, ProteinReader proteinReader) throws UpdateableException,EvaluationException,AlignmentException {
        ProteinInfoList proteinInfoList = new ProteinInfoList();
        ProteinInfo proteinInfo;
        int i = 0;
        for (ProteinMetaData fileData : filesData) {
            System.out.print(".");
            if ((i++) % 50 == 0) System.out.println(" " + i + " ");

            proteinInfo = null;
            for (ProteinAnalyzer analyzer : analyzers) {
                if (proteinInfo == null) {
                    proteinInfo = analyze(fileData, analyzer, proteinReader);
                    if (proteinInfo == null) break;
                } else {
                    ProteinInfo tmp = analyze(fileData, analyzer, proteinReader);
                    if (tmp == null) {
                        proteinInfo = null;
                        break;
                    }
                    proteinInfo.append(tmp);
                }
            }
            if (proteinInfo != null)
                proteinInfoList.add(proteinInfo);
            else {
                System.out.println("Ignoring " + fileData);
            }
        }
        return proteinInfoList;
    }

    public ProteinInfo analyze(ProteinMetaData proteinMetaData, ProteinAnalyzer analyzer, ProteinReader proteinReader) throws UpdateableException, EvaluationException,AlignmentException {
        return analyzer.analyze(proteinReader.readProtein(proteinMetaData));
    }


}
