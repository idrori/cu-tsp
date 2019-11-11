package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.ExcludedVolCreator;
import meshi.energy.pairwiseNonBondedTerms.ExcludedVol.SwellCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.stretch.Stretch;
import meshi.energy.simpleEnergyTerms.stretch.StretchAllCreator;
import meshi.energy.simpleEnergyTerms.stretch.StretchCreator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 30/04/14
 * Time: 17:27
 * To change this template use File | Settings | File Templates.
 */
public class Unknot extends MeshiProgram{
    private static EnergyCreator[] energyCreators1 = {
            new BondCreator(),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new ExcludedVolCreator()
    };
    private static StretchCreator stretchCreator = new StretchCreator();
    private static EnergyCreator[] energyCreators2 = {
            stretchCreator,
            new StretchAllCreator(),
            new SwellCreator(),
            new BondCreator("simple"),
            new AngleCreator(),
            new PlaneCreator(),
            new OutOfPlaneCreator(),
            new ExcludedVolCreator()
    };

    public static void main(String[] args) throws EvaluationException, UpdateableException, AlignmentException,IOException{
        double meanDistance;
        ArrayList<SummaryElement> summary = new ArrayList();
        initRandom(0);
        CommandList commands = new CommandList(args, new CommandsException("This is weird"));
        String targetName = args[1];
        File   targetDir = new File(targetName);
        if ((!targetDir.exists()) || (!targetDir.isDirectory())) throw new RuntimeException("This is weird");
        MeshiWriter suspects = new MeshiWriter(targetName+"/suspects.txt");
        MeshiWriter logWriter = new MeshiWriter(targetName+"/"+targetName+".findKnots.log");

        File[] files = targetDir.listFiles();
        for (File file : files) {
            if (file.getAbsolutePath().endsWith(".pdb")) {
                System.out.println("Now stretching "+file.getAbsolutePath());
                Protein model = Utils.getProtein(commands, file.getAbsolutePath(), new ResidueExtendedAtomsCreator(), Utils.defaultExceptionHandler);
                model.atoms().sideChains().setNowhere();

                try {
                    Utils.relax(model.atoms(), model, energyCreators1, commands);
                    for (Atom atom : model.atoms()) {
                        if ((atom.name().equals("CB")) ||
                                (atom.name().equals("O")) ||
                                (atom.name().equals("H")))
                            atom.setNowhere();
                     }
                    Utils.relax(model.atoms(), model, energyCreators2, commands);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    suspects.println(file.getAbsolutePath() + " failed to minimize");
                    suspects.flush();
                    continue;
                }
                meanDistance = ((Stretch) stretchCreator.term()).getDistance().distance()/model.chain().numberOfNonDummyResidues();
                if (args.length > 2) {
                    MeshiWriter writer = new MeshiWriter(args[2]);
                    writer.println("End to end distance = " + meanDistance);
                    model.atoms().located().print(writer);
                    writer.close();
                }
                summary.add(new SummaryElement(file.getAbsolutePath(), meanDistance));
                logWriter.println(file.getAbsoluteFile() + " " + meanDistance);
                logWriter.flush();
            }
        }
        logWriter.println("------------------------");
        SummaryElement.addZscore(summary);
        SummaryElement[] sortedSummary = summary.toArray(new SummaryElement[summary.size()]);
        Arrays.sort(sortedSummary);
        for (SummaryElement summaryElement : sortedSummary) {
            logWriter.println(summaryElement);
        }
        if (sortedSummary[0].zScore < -2) {
              for (SummaryElement summaryElement : sortedSummary) {
                  if (summaryElement.zScore < -2)
                      suspects.println(summaryElement);
              }
            suspects.close();
        }

        logWriter.close();
    }

    private static class SummaryElement implements Comparable<SummaryElement> {
        String fileName;
        double endToEnd;
        double zScore;

        public SummaryElement(String fileName, double endToEnd) {
            this.fileName = fileName;
            this.endToEnd = endToEnd;
        }

        public String toString() {
            return fileName+" "+endToEnd+" "+zScore;
        }

        public int compareTo(SummaryElement other) {
            if (zScore > other.zScore) return 1;
            if (zScore == other.zScore) return 0;
            return -1;
        }
        public static void addZscore(ArrayList<SummaryElement> summary) {
            double sum = 0;
            double sum2 = 0;
            for (SummaryElement summaryElement : summary) {
                sum  += summaryElement.endToEnd;
                sum2 += summaryElement.endToEnd * summaryElement.endToEnd;
            }
            double mean  = sum/summary.size();
            double mean2 = sum2/summary.size();
            double std = Math.sqrt(mean2 - mean * mean);
            for (SummaryElement summaryElement : summary) {
                summaryElement.zScore = (summaryElement.endToEnd - mean)/std;
            }
        }
    }
}
