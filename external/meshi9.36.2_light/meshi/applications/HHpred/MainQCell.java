package meshi.applications.HHpred;
import meshi.sequences.AlignmentCell;
import meshi.sequences.SequenceAlignmentCell;

public class MainQCell extends SequenceAlignmentCell {

	private AlignmentCell mainQcell;

	public MainQCell(char c, int number) {
		super(c, number);
	}

	public String toString() {
		return "MainQCell [" + super.toString() + "]";
	}

}