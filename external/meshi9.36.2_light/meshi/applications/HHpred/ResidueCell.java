package meshi.applications.HHpred;

import meshi.parameters.ResidueType;
import meshi.sequences.AlignmentCell;

public class ResidueCell extends AlignmentCell {

	// private ResidueType restype;
	public ResidueCell(ResidueType restype, int num) {
		super(restype, num);
	}

	@Override
	public boolean gap() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String toString() {
		return "ResidueCell [number=" + number + ", obj=" + obj + ", comment="
				+ comment + ", gap()=" + gap() +
				// ", object()=" + object()
				// + ", number()=" + number() +
				// ", column()=" + column() +
				// ", getClass()=" + getClass() + ", hashCode()=" + hashCode()
				// + ", toString()=" + super.toString() +
				"]";
	}

}
