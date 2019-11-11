package meshi.energy.templateEnergy;
import java.util.*;
import meshi.util.*;

public class TemplateEnergyElementsList extends ArrayList<TemplateEnergyElement> implements Updateable{
    private int numberOfUpdates;
    public TemplateEnergyElementsList() {
	super();
	numberOfUpdates = 0;
    }

    /*------------------------------------------ update --------------------------------------------*/
    /**
     * Updates the distances. 
     **/
    public void update(int numberOfUpdates) throws UpdateableException {
    if (numberOfUpdates == this.numberOfUpdates+1) {
        this.numberOfUpdates++;
        update();
    }
    else if (numberOfUpdates != this.numberOfUpdates)
        throw new UpdateableException();
    }
    public void update() { 
	for (TemplateEnergyElement element:this) element.update();
    }

     public void setNumberOfUpdates(int numberOfUpdates){
        this.numberOfUpdates = numberOfUpdates;
    }
}
