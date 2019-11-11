package meshi.energy.templateEnergy;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.util.*;
import meshi.sequences.*;
import meshi.util.info.InfoType;

import java.util.*;


public class TemplateEnergyCreator  extends EnergyCreator {
    public final double templateThreshold;
    public  enum Template {FIRST,SECOND};
    public final Template template;
    public final ResidueAlignment residueAlignment;
    public TemplateEnergyCreator(ResidueAlignment residueAlignment, Template template, double templateThreshold) {
	super(InfoType.TEMPLATE_ENERGY);
	this.residueAlignment  = residueAlignment;
	this.template          = template;
	this.templateThreshold = templateThreshold;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
	Residue templateResidue;
	Residue proteinResidue;
	term = new TemplateEnergy(new EnergyInfoElement(infoType,"Template energy"));
	
	for(Iterator columns = residueAlignment.iterator(); columns.hasNext();) {
	    ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
	    if (template == Template.FIRST) {
		templateResidue = (Residue) column.cell(1).obj;
		proteinResidue = (Residue)  column.cell(0).obj;
	    }
	    else {
		templateResidue = (Residue) column.cell(1).obj;
		proteinResidue = (Residue)  column.cell(0).obj;
	    }
	    if (! templateResidue.dummy()) {
		if (templateResidue.ca().distanceFrom(proteinResidue.ca()) < templateThreshold){
		    ((TemplateEnergy)term).add(column);
		}
	    }
	}
    	return term;
    }
}
   
