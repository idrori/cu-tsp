package meshi.applications.HHpred;

import java.util.Arrays;

import meshi.parameters.ResidueType;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentCell;
import meshi.sequences.AlignmentColumn;

public class HhpredAlignmentColumn extends AlignmentColumn {
	public HhpredAlignmentColumn(int numberOfRows) {
		super(numberOfRows);
	}

	public HhpredAlignmentColumn(MainQCell cell0, Q1cell cell1, Q2cell cell2,
			Q3cell cell3, MainTCell cell4, T1cell cell5, T2cell cell6,
			T3cell cell7, Tsspr cell8, Qsspr cell9, Tssdsspr cell10,
			ResidueCell cell11, ResidueCell cell12) {
		super(13);
		cells[0] = cell0;
		cells[1] = cell1;
		cells[2] = cell2;
		cells[3] = cell3;
		cells[4] = cell4;
		cells[5] = cell5;
		cells[6] = cell6;
		cells[7] = cell7;
		cells[8] = cell8;
		cells[9] = cell9;
		cells[10] = cell10;
		cells[11] = cell11;
		cells[12] = cell12;

	}

	public HhpredAlignmentColumn(
			HhpredSequenceAlignmentColumn hhpredSequenceAlignmentColumn,
			String[] seqCh) {
		super(13);
		char tssprch = seqCh[8].charAt(0);
		char qssprch = seqCh[9].charAt(0);
		char tssdsspch = seqCh[10].charAt(0);
		int qssprnum = hhpredSequenceAlignmentColumn.cell(0).number();
		int tssprnum = hhpredSequenceAlignmentColumn.cell(4).number();
		int tssdsspnum = tssprnum;
		SecondaryStructure qss = SecondaryStructure.secondaryStructure(qssprch);
		SecondaryStructure tss = SecondaryStructure.secondaryStructure(tssprch);
		SecondaryStructure tssdss = SecondaryStructure
				.secondaryStructure(tssdsspch);
		Qsspr qssp = new Qsspr(qss, qssprnum);
		Tsspr tssp = new Tsspr(tss, tssprnum);
		Tssdsspr tssdssp = new Tssdsspr(tssdss, tssdsspnum);
		ResidueType qrt = ResidueType.type(hhpredSequenceAlignmentColumn
				.getChar(0));// res type gap?
		ResidueType trt = ResidueType.type(hhpredSequenceAlignmentColumn
				.getChar(4));
		ResidueCell qrtc = new ResidueCell(qrt, hhpredSequenceAlignmentColumn
				.cell(0).number());
		ResidueCell trtc = new ResidueCell(trt, hhpredSequenceAlignmentColumn
				.cell(4).number());
		setColumn(hhpredSequenceAlignmentColumn, tssp, qssp, tssdssp, qrtc,
				trtc);
	}

	private void setColumn(
			HhpredSequenceAlignmentColumn hhpredSequenceAlignmentColumn,
			Tsspr tssp, Qsspr qssp, Tssdsspr tssdssp, ResidueCell qrtc,
			ResidueCell trtc) {

		for (int i = 0; i < hhpredSequenceAlignmentColumn.size(); i++) {
			this.cells[i] = hhpredSequenceAlignmentColumn.cell(i);
		}
		cells[8] = tssp;
		cells[9] = qssp;
		cells[10] = tssdssp;
		cells[11] = qrtc;
		cells[12] = trtc;

	}

	@Override
	public String toString() {
		return "HhpredAlignmentColumn [ cells=" + Arrays.toString(cells)
				+ " ]\n";
	}

	public AlignmentCell getQCell() {
		return cells[11];
	}

	public AlignmentCell getTCell() {
		return cells[12];
	}
}
