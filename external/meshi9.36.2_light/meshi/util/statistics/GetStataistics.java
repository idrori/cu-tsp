package meshi.util.statistics;

import meshi.energy.goap.Goap;
import meshi.energy.goap.GoapCreator;
import meshi.energy.goap.GoapInfo;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.ResidueType;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.dssp.DSSP;
import meshi.util.dssp.DsspReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.*;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;
import meshi.util.info.ProteinInfoList;
import meshi.util.proteinReaders.ProteinReader;
import meshi.util.proteinReaders.StandardProteinReader;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 */
public class GetStataistics {

    public static void getStatistics(String[] args, ProteinAnalyzer[] proteinAnalyzers) throws Exception {
        String inputDirectoryName = args[1];
        String outputDirecoryName = args[2];
        String nativDecoysFlag    = args[3];

        DSSPreader    dsspReader       = new DSSPreader(inputDirectoryName,nativDecoysFlag);
        ProteinReader proteinReader    = new StandardProteinReader(dsspReader,args[0]);
        File          rootDirectory    = new File(inputDirectoryName);
        FileFilter    filter           = new PDBfileFilter();
        File[]        dataDirectories  = getDataDirectories(rootDirectory, nativDecoysFlag);
        File          outputDir        = new File(outputDirecoryName);

        if (!outputDir.exists())
            outputDir.mkdir();

        for (File dataDirectory : dataDirectories) {
            PdbFilesScanner scanner  = new PdbFilesScanner(dataDirectory, filter, nativDecoysFlag);
            ProteinInfoList infoList = scanner.analyze(proteinAnalyzers, proteinReader);
            if(nativDecoysFlag.equals("DECOYS")) {
                String outputFileName = outputDirecoryName + "\\" + dataDirectory.getName() + ".xml";
                File testFile = new File(outputFileName);
                if (testFile.exists()) {
                    System.out.println(outputFileName + " exists nothing to do.");
                    continue;
                }
                MeshiWriter writer = new MeshiWriter(outputFileName);

                infoList.writeXml(writer, nativDecoysFlag);
                writer.close();
            }
        }
    }

    private static File[] getDataDirectories(File rootDirectory, String nativDecoysFlag) {
        File[] out;
        if (nativDecoysFlag.equals("NATIVE")) {
            out = new File[1];
            out[0] = new File(rootDirectory+"\\chains");
//            out[0] = new File(rootDirectory+"\\debug");
        }
        else {
            if (nativDecoysFlag.equals("DECOYS")){
                out = rootDirectory.listFiles(new TargetDirectories());
            }
            else throw new RuntimeException("Weird nativDecoysFlag "+nativDecoysFlag);
        }
        return out;
    }

    private static class TargetDirectories implements FileFilter {
        public boolean accept(File file) {
            if (!file.isAbsolute()) return false;
            return (file.getName().startsWith("T0") ||  file.getName().startsWith("R0"));
        }
    }

    public static String getDsspFileName(Protein protein) {
        String dsspFileName;
        String fileName = (String) protein.metaData().get(ProteinMetaData.MetaDataKey.FILE_NAME);
        if (fileName.length() < 20)
            dsspFileName = "dssp/" + protein.name() + ".dssp";
        else dsspFileName = "dssp/" + protein.name().substring(0,8)+"/"+protein.name() + ".dssp";
        return dsspFileName;
    }


    private static class ContactsStatistics implements ProteinAnalyzer {
        enum Types {H_SC_SS(new HydrophobicSideChainSecondaryStructure()), H_SC_COIL(new HydrophobicSideChainCoil()),
                    P_SC_SS(new PolarSideChainSecondaryStructure())      , P_SC_COIL(new PolarSideChainCoil()),
                    BB_SS(new BackboneSecondaryStructure())              , BB_COIL(new BackboneCoil()),
                    OTHER(new OtherFilter());
            Filter filter;
            Types(Filter filter) {
                this.filter = filter;
            }
        }
        double threshold;

        public ContactsStatistics(double threshold) {
            this.threshold = threshold;
        }

        public ProteinInfo analyze(Protein protein) {
            AtomList atoms = protein.atoms();
            Types[] types = new Types [atoms.size()];
            for (int i = 0; i < atoms.size(); i++) types[i] = null;
            int[] typeCounters = new int[Types.values().length];
            int[] typeCounters2 = new int[Types.values().length];

            for (int i = 0; i < Types.values().length; i++) typeCounters[i] = typeCounters2[i] = 0;

            for (int iAtom = 0; iAtom <atoms.size(); iAtom++)
                for (Types type : Types.values())
                    if (type.filter.accept(atoms.get(iAtom))) {
                        if (types[iAtom] == null)
                            types[iAtom] = type;
                        else throw new RuntimeException("Type conflict "+types[iAtom]+"  "+
                                                         type+"\n"+atoms.get(iAtom)+" "+atoms.get(iAtom).type()+" "+atoms.get(iAtom).type().isCarbon());
                        typeCounters[type.ordinal()]++;
                        if (!atoms.get(iAtom).nowhere())
                            typeCounters2[type.ordinal()]++;
                    }

            for (int i = 0; i < Types.values().length; i++)
                    if (typeCounters2[i] < 5) {
                        System.out.println("Protein "+protein.metaData().get(ProteinMetaData.MetaDataKey.NAME)+" "+protein.atoms().size()+" ignored atomType "+Types.values()[i]+" # "+typeCounters[i]);
                        return null;
                    }


            //System.out.println("analysing " + protein.name());

            int[][] contact = getContacts(atoms, types, threshold);

            ArrayList<MeshiInfo> infoList = new ArrayList();
            for (Types type : Types.values()) {
                infoList.add(new MeshiInfo(null, new Integer(typeCounters[type.ordinal()]),type.toString()));
            }
            for (int iType = 0; iType < Types.values().length; iType++)
                for (int jType = 0; jType <= iType; jType++) {
                    int index = iType+((jType+1)*jType)/2;
                    infoList.add(new MeshiInfo(null, new Integer(contact[index][0]),
                                               Types.values()[iType]+"_"+Types.values()[jType]+"_"+threshold));
                    infoList.add(new MeshiInfo(null, new Integer(contact[index][1]),
                            Types.values()[iType]+"_"+Types.values()[jType]+"_"+threshold+"_2"));
                }
            ProteinInfo proteinInfo = new ProteinInfo(protein.metaData(), infoList, protein);
            return proteinInfo;
        }



        public static int[][] getContacts(AtomList atoms, Types[] types, double threshold) {
            double iX, iY, iZ, distance2, dx, dy, dz;
            double threshold2 = threshold*threshold;

            int nContactCombinations = (Types.values().length*(Types.values().length+1))/2;
            int[][] out = new int[nContactCombinations][2];
            int[][] contacts = new int[atoms.size()][Types.values().length];

            for (int iAtom = 0; iAtom <atoms.size(); iAtom++)
                for (int j = 0;  j <Types.values().length; j++)
                    contacts[iAtom][j] = 0;

            for (int iAtom = 1; iAtom < atoms.size(); iAtom++) {
                Atom atomI = atoms.get(iAtom);
                if (atomI.nowhere()) continue;
                Residue residueI = atomI.residue();
                iX = atomI.x();
                iY = atomI.y();
                iZ = atomI.z();
                for (int jAtom = 0; jAtom < iAtom; jAtom++) {
                    Atom atomJ = atoms.get(jAtom);
                    if (atomJ.nowhere()) continue;
                    int iType = types[iAtom].ordinal();
                    Residue residueJ = atomJ.residue();
                    if (residueI == residueJ) continue;
                    dx = iX - atomJ.x();
                    dy = iY - atomJ.y();
                    dz = iZ - atomJ.z();
                    distance2 = dx * dx + dy * dy + dz *dz;
                    if (distance2 > threshold2) continue;
                    int jType = types[jAtom].ordinal();
                    contacts[iAtom][jType]++;
                    contacts[jAtom][iType]++;
                }
            }

            for (int iAtom = 1; iAtom <atoms.size(); iAtom++) {
                int iType = types[iAtom].ordinal();
                int index;
                for (int jType = 0; jType< Types.values().length; jType++) {
                    //Utils.printDebug("XXXXXXXXXX "," "+iAtom+" "+iType);
                    if (iType >= jType)
                         index = jType + ((iType + 1) * iType) / 2;
                    else index = iType + ((jType + 1) * jType) / 2;
                    //Utils.printDebug("xxxxxxxxxx ",""+index+" "+iType+" "+jType);
                    out[index][0] += contacts[iAtom][jType];
                    out[index][1] += contacts[iAtom][jType] * contacts[iAtom][jType];
                }
            }

                return out;
        }
    }

    public static ArrayList<MeshiInfo> getResidueNumbers(Protein protein) {
        InfoType[] types = {InfoType.ALA_NUMBER, InfoType.CYS_NUMBER, InfoType.ASP_NUMBER, InfoType.GLU_NUMBER, InfoType.PHE_NUMBER, InfoType.GLY_NUMBER,
                InfoType.HIS_NUMBER, InfoType.ILE_NUMBER, InfoType.LYS_NUMBER, InfoType.LEU_NUMBER, InfoType.MET_NUMBER, InfoType.ASN_NUMBER,
                InfoType.PRO_NUMBER, InfoType.GLN_NUMBER, InfoType.ARG_NUMBER, InfoType.SER_NUMBER, InfoType.THR_NUMBER, InfoType.VAL_NUMBER,
                InfoType.TRP_NUMBER, InfoType.TYR_NUMBER};

        ArrayList<MeshiInfo> infoList = new ArrayList();
        int[] counts = new int[ResidueType.values().length];
        for (int i = 0; i < counts.length; i++) counts[i] = 0;

        for (Residue residue : protein.residues()) {
            counts[residue.type.ordinal()]++;
        }

        for (int iType = 0; iType < 20; iType++)
            infoList.add(new MeshiInfo(types[iType], new Integer(counts[iType]), "The number of " + ResidueType.values()[iType] + " in protein"));

        return infoList;
    }

    private static class DSSPreader implements DsspReader {
        boolean nativeFlag;
        String  rootDirectoryName;
        public DSSPreader(String rootDirectoryName, String nativeDecoyFlag){
            this.nativeFlag    =  nativeDecoyFlag.equals("NATIVE");
            this.rootDirectoryName = rootDirectoryName;
        }
        public DSSP read(Protein protein) {
            File dsspFile;
            if (nativeFlag) {
                dsspFile = new File(rootDirectoryName+"\\dsspFiles\\"+protein.name()+"." + "dssp");
            }
            else {
                dsspFile = new File(rootDirectoryName+"\\..\\dsspFiles\\"+
                        ((String) protein.metaData().get(ProteinMetaData.MetaDataKey.FILE_NAME)).substring(0,8)+"\\"+
                                    protein.metaData().get(ProteinMetaData.MetaDataKey.FILE_NAME)+".dssp");
            }
            return new DSSP(dsspFile.getAbsolutePath());
        }


    }
}
