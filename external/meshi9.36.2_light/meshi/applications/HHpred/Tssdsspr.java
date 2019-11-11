package meshi.applications.HHpred;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentCell;

public class Tssdsspr extends AlignmentCell {

	/**
	 * @param args
	 */

	private SecondaryStructure secStruct;

	public Tssdsspr(SecondaryStructure secStruct, int num) {
		super(secStruct, num);
	}

	@Override
	public boolean gap() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public String toString() {
		return "Tssdsspr [secStruct=" + super.obj.toString() + "]";
	}

}
