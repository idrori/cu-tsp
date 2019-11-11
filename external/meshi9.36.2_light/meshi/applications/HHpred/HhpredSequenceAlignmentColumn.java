package meshi.applications.HHpred;
import java.util.Arrays;

import meshi.sequences.SequenceAlignmentCell;
import meshi.sequences.SequenceAlignmentColumn;

public class HhpredSequenceAlignmentColumn extends SequenceAlignmentColumn {

	public HhpredSequenceAlignmentColumn() {
		super(8);
	}

	public HhpredSequenceAlignmentColumn(int numberOfRows) {
		super(numberOfRows);

	}

	public HhpredSequenceAlignmentColumn(String resNumline, String line,
			int numOfRow, HhpredAlignment hhpredAlignment,
			HhpredSequenceAlignmentColumn lastColumn) {
		super(8);
		String[] seqCh = line.split("	");
		String[] nums = resNumline.split("	");
		int[] resNums = new int[8];
		for (int i = 0; i < 8; i++) {
			resNums[i] = Integer.parseInt(nums[i]) + numOfRow - 1;
		}
		setColumn(seqCh, resNums, lastColumn);
		HhpredAlignmentColumn halc = new HhpredAlignmentColumn(this, seqCh);
		hhpredAlignment.add(halc);
	}

	private void setColumn(String[] seqCh, int[] resNums,
			HhpredSequenceAlignmentColumn lastColumn) {
		SequenceAlignmentCell cell = null;
		final int COLUMNS = 8;
		char resChar = ' ';
		int resNum = 0;
		for (int i = 0; i < COLUMNS; i++) {
			if (i != 0)
				this.add(i - 1, cell);
			cell = null;
			resChar = seqCh[i].charAt(0);
			if (resChar == '.')
				resChar = '-';
			if (isGap(resChar)) {
				if (lastColumn == null)
					resNum = resNums[i] - 1;
				else
					resNum = lastColumn.cell(i).number();
			} else {
				if (lastColumn == null)
					resNum = resNums[i];
				else
					resNum = lastColumn.cell(i).number() + 1;
			}

			switch (i) {
			case 0:
				cell = new MainQCell(resChar, resNum);
				break;
			case 1:
				cell = new Q1cell(resChar, resNum);
				break;
			case 2:
				cell = new Q2cell(resChar, resNum);
				break;
			case 3:
				cell = new Q3cell(resChar, resNum);
				break;
			case 4:
				cell = new MainTCell(resChar, resNum);
				break;
			case 5:
				cell = new T1cell(resChar, resNum);
				break;
			case 6:
				cell = new T2cell(resChar, resNum);
				break;
			case 7:
				cell = new T3cell(resChar, resNum);
				break;
			}
		}
		this.add(7, cell);
		/*
		 * for(int i=0;i<8;i++){ System.out.println(this.cell(i).toString()); }
		 */
	}

	private boolean isGap(char res) {
		boolean ans = false;
		if (res == '-' || res == '.')
			ans = true;
		return ans;
	}

	@Override
	public String toString() {
		return "HhpredSequenceAlignmentColumn [cells=" + Arrays.toString(cells)
				+ "]\n";
	}

}
