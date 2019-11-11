package meshi.scoringFunctions;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

import java.util.ArrayList;

/**
 * Created by chen on 29/04/2016.
 */
public class MetaScore extends GdtScore {
    private final ArrayList<Score> scores;
    public MetaScore(String parametersFileName, String name, ArrayList<Score> scores){
        super(parametersFileName, name);
        this.scores = scores;
    }

    public MeshiInfo score(MeshiInfo energyInfo) {
        ArrayList<MeshiInfo> infos = new ArrayList();
        for (Score score : scores) {
            if (score instanceof MetaScore) continue;
            infos.add(score.score(energyInfo));
        }
        for (MeshiInfo info : infos) {
            energyInfo.getChildren().add(info);
        }
        MeshiInfo out =super.score(energyInfo);
        for (MeshiInfo element : out.flatten())
            if (element.type != InfoType.INTERDECILE) {
                if (((Double) element.getValue()).doubleValue() >= 1) element.setValue(Double.valueOf(1));
                if (((Double) element.getValue()).doubleValue() <= 0.11) element.setValue(Double.valueOf(0.11));
            }
        return out;
    }



}
