package meshi.applications.HHpred;

import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentCell;

public class Tsspr extends AlignmentCell {
	private SecondaryStructure secStruct;

	public Tsspr(SecondaryStructure secStruct, int num) {
		super(secStruct, num);
	}

	@Override
	public boolean gap() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String toString() {
		return "Tsspr [secStruct=" + super.obj.toString() + "]";
	}

}
