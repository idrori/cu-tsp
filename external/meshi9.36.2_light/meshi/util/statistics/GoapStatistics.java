package meshi.util.statistics;

import meshi.energy.goap.Goap;
import meshi.energy.goap.GoapCreator;
import meshi.energy.goap.GoapInfo;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.ProteinAnalyzer;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;

import java.util.ArrayList;

/**
 * Created by chen on 30/12/2015.
 */
public class GoapStatistics implements ProteinAnalyzer {
    CommandList commands;
    GoapCreator goapCreator;

    public GoapStatistics(CommandList commands) {
        this.commands = commands;
        goapCreator = new GoapCreator();
    }

    public ProteinInfo analyze(Protein protein) {
        ArrayList<MeshiInfo> infoList = new ArrayList();

        Goap goap = (Goap) goapCreator.createEnergyTerm(protein, null, commands);
        GoapInfo tempInfo = (GoapInfo) goap.evaluate();

        ArrayList<MeshiInfo> residueNumbers = GetStataistics.getResidueNumbers(protein);

        infoList.addAll(residueNumbers);
        infoList.add(new MeshiInfo(InfoType.D_FIRE, tempInfo.dFIRE.getValue(), "dFire (Goap) "));
        infoList.add(new MeshiInfo(InfoType.GOAP_AG, tempInfo.goapAG.getValue(), "GOAP_AG (Goap) "));
        infoList.add(new MeshiInfo(InfoType.GOAP, tempInfo.getValue(), "Goap"));
        ProteinInfo proteinInfo = new ProteinInfo(protein.metaData(), infoList, protein);
        return proteinInfo;
    }
}
