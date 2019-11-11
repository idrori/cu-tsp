package meshi.applications.HHpred;

import meshi.sequences.AlignmentCell;
import meshi.sequences.SequenceAlignmentCell;

public class MainTCell extends SequenceAlignmentCell {

	private AlignmentCell mainTcell;

	public MainTCell(char c, int number) {
		super(c, number);
	}

	@Override
	public String toString() {
		return "MainTCell [" + super.toString() + "]";
	}
}
