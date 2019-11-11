package programs;

import meshi.energy.KB2013.KbStatistics;
import meshi.energy.KB2013.Separator;
import meshi.util.UpdateableException;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 11/07/13
 * Time: 17:58
 * To change this template use File | Settings | File Templates.
 */
public class GetKBstatistics {
    public static void main(String[] argv) throws IOException,SAXException,ParserConfigurationException,UpdateableException {
        String fileList = argv[0];
        String outputFile = argv[1];
        String separatorString = argv[2];
        double maxDistance = Double.parseDouble(argv[3]);
        double binSize   = Double.parseDouble(argv[4]);
        Separator separator = Separator.getSeparator(separatorString);
        KbStatistics kbs = new KbStatistics(separator,maxDistance,binSize);
        kbs.learn(fileList);
        kbs.save(argv[1]);
    }
}
