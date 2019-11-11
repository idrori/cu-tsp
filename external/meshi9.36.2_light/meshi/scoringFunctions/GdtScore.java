package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.UpdateableException;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 02/02/13
 * Time: 14:47
 * To change this template use File | Settings | File Templates.
 */
public class GdtScore extends CombinedEnergyScore  {
    public GdtScore(String parametersFileName, String name){
        super(parametersFileName, name);
        }
    public MeshiInfo score(MeshiInfo energyInfo) {
        MeshiInfo out =super.score(energyInfo);
        for (MeshiInfo element : out.flatten())
        if (element.type != InfoType.INTERDECILE) {
            if (((Double) element.getValue()) >= 1)    element.setValue(Double.valueOf(1));
            if (((Double) element.getValue()) <= 0.11) element.setValue(Double.valueOf(0.11));
        }
        return out;
    }

    public static ArrayList<Score> getScoreFunctions(CommandList commands) {
        boolean exists = commands.keyExists("selectionScore");
        if (!exists) return null;
        CommandList scoreCommands = commands.firstWordFilter("selectionScore");
        if (scoreCommands.size() == 0) return null;
        ArrayList<Score> out = new ArrayList<Score>(scoreCommands.size());
        for (Command command:scoreCommands) {
            GdtScore gdtScore = new GdtScore(command.secondWord(),command.thirdWord());
            out.add(gdtScore);
        }

        exists = commands.keyExists("metaSelectionScore");
        if (exists) {
            scoreCommands = commands.firstWordFilter("metaSelectionScore");
            if (scoreCommands.size() == 0) return null;
            ArrayList<Score> tmp = new ArrayList<Score>(scoreCommands.size());
            for (Command command : scoreCommands) {
                GdtScore metaScore = new MetaScore(command.secondWord(), command.thirdWord(), out);
                tmp.add(metaScore);
            }

            if (tmp.size()!=0)
                out.addAll(tmp);
        }
        return out;
    }
}
