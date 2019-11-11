package meshi.applications.HHpred;

import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentCell;

public class Qsspr extends AlignmentCell {
	private SecondaryStructure secStruct;

	public Qsspr(SecondaryStructure secStruct, int num) {
		super(secStruct, num);
	}

	@Override
	public boolean gap() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String toString() {
		return "Qsspr [secStruct=" + super.obj.toString() + "]";
	}

}