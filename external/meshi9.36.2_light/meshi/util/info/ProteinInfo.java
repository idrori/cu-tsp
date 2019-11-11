package meshi.util.info;

import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;

import java.util.ArrayList;

/**
 * Created by chen on 03/12/2015.
 */
public class ProteinInfo extends MeshiInfoElementList{
    ProteinMetaData proteinMetaData;
    ArrayList<MeshiInfo> infoList;
    Protein protein;

    public ProteinInfo(ProteinMetaData proteinMetaData, ArrayList<MeshiInfo> infoList, Protein protein){
        this.proteinMetaData = proteinMetaData;
        this.infoList        = infoList;
        this.protein         = protein;
    }

    public void append(ProteinInfo other) {
        for (ProteinMetaData.MetaDataKey key : ProteinMetaData.MetaDataKey.values()) {
            if (!proteinMetaData.get(key).equals(other.proteinMetaData.get(key)))
                throw new RuntimeException("Error in append. Incomatible metaData fields.\n" + proteinMetaData + "\n" + other.proteinMetaData);
        }
        infoList.addAll(other.infoList);
    }

    public String toXml(int indentationTabs)  {
        return toXml(indentationTabs, null);
    }

    public String toXml(int indentationTabs, ChainInfo residueInfos)  {
            ProteinInfoHeader header = new ProteinInfoHeader((String) proteinMetaData.get(ProteinMetaData.MetaDataKey.TARGET),
                (String) proteinMetaData.get(ProteinMetaData.MetaDataKey.FILE_NAME)," dummyComment", protein);
        if (residueInfos != null)
            for (ResidueInfo element : residueInfos)
                add(element);
        return super.toXML(indentationTabs, header, "</" + header.type + ">");
    }

    /*{
        String out = indentation+"<ProteinInfo ";
        for (ProteinMetaData.MetaDataKey key : ProteinMetaData.MetaDataKey.values()) {
            Object value = proteinMetaData.get(key);
            out += key+"=\""+value+"\"  ";
        }
        out += ">\n";
        if (infoList != null) {
            for (MeshiInfo info : infoList) {
                out += info.toXml(indentation+"     ")+"\n";
            }
        }
        out += indentation+"</ProteinInfo>\n";
        return out;
    } */
}
