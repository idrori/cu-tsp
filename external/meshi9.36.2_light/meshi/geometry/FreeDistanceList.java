package meshi.geometry;

import meshi.util.Updateable;
import meshi.util.Utils;

import java.util.ArrayList;

/**
 *
 */
public class FreeDistanceList extends ArrayList<FreeDistance> implements Updateable{
    private int numberOfUpdates;

    public FreeDistanceList(){
        numberOfUpdates = 0;
    }
    public void setNumberOfUpdates(int numberOfUpdates){
        this.numberOfUpdates = numberOfUpdates;
    }
    public void update(int numberOfUpdates) {
        Utils.checkUpdate(this.numberOfUpdates,numberOfUpdates,this);
        this.numberOfUpdates++;
        for (FreeDistance d:this)
            d.update();
    }
}
