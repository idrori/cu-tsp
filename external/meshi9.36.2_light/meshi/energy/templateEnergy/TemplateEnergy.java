package meshi.energy.templateEnergy;
import meshi.applications.prediction.beautify.BeautifyAttribute;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.util.Updateable;
import meshi.util.UpdateableException;

import java.util.Iterator;

public class TemplateEnergy  extends AbstractEnergy implements Updateable{
    private static TemplateEnergyElementsList elements;
    
    public TemplateEnergy(EnergyInfoElement info) {
	super(toArray(elements = new TemplateEnergyElementsList()), info );
	comment = "TemplateEnergy";
    }
    public void update(){}
    public void update(int i ){
	try {
	    elements.update(i);
	} catch (UpdateableException ex) {throw new RuntimeException("Update problem in TemplateEnergy:\n"+ex);}
    }

    public void ignoreLoopResidues() {
	for (TemplateEnergyElement element:elements) {
	    if (BeautifyAttribute.isLoop(element.residue())) element.off();
	}
    }
    public void ignoreUnhappyResidues(double threshold){
	for (TemplateEnergyElement element:elements) 
	    element.offIfUnhappy(threshold);
    }

    public EnergyInfoElement evaluate() {
	if (!on) return info;
	double e = 0;
	for (TemplateEnergyElement element:elements) 
	    e += element.evaluate();
        info.setValue(e);
	return info;
    }
    
    public void test(TotalEnergy totalEnergy,Atom atom){
        if (! on) System.out.println(""+this +" is off");
         for(TemplateEnergyElement element: elements){
	    if(!element.frozen()){
                element.test(totalEnergy,atom);
	    }
	}
    }
 
    public void evaluateAtoms() throws EvaluationException{
	if (on) {
	    for(TemplateEnergyElement element: elements)
		element.evaluateAtoms();
	}
    }
    
    public void add(ResidueAlignmentColumn column) {
         Residue resBeautify = (Residue)column.cell0().obj;
         Residue resShotgun = (Residue)column.cell1().obj;
         if ((! resShotgun.dummy()) &&
             (! resBeautify.dummy())) {
                Atom caBeautify = resBeautify.ca();
                Atom caShotgun = resShotgun.ca();
                elements.add(new TemplateEnergyElement(caShotgun, caBeautify, weight));
        }
    }

    public void add(ResidueAlignment residueAlignment) {
	    for (Iterator columns = residueAlignment.iterator(); columns.hasNext();) 
		    add((ResidueAlignmentColumn) columns.next());
    }
}

