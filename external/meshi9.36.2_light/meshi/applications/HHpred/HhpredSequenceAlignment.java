package meshi.applications.HHpred;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import meshi.util.string.StringList;

public class HhpredSequenceAlignment extends
		ArrayList<HhpredSequenceAlignmentColumn> {

	public final StringList comments;

	public HhpredSequenceAlignment() {
		super();
		comments = new StringList();
	}

	public HhpredSequenceAlignment(File file, HhpredAlignment hhpredAlignment)
			throws FileNotFoundException {
		this();
		BufferedReader buffreader;
		Scanner scanner = null;
		buffreader = new BufferedReader(new FileReader(file.getAbsolutePath()));
		linesParsing(scanner, buffreader, hhpredAlignment);
	}

	// *******************************hhpred format
	// example********************************************
	// Model
	// NAME = 2a8e_A
	// qchain qchain1 qchain2 qchain3 tchain tchain1 tchain2 tchain3 t_sspr
	// q_sspr t_ssdssp
	// 15 1 1 1 5 5 9 5
	// Y - - - R G G T C c S
	// L - W F F F F F C c C
	// L - T N T E V T c c C
	// ************************************************************************************************

	private void linesParsing(Scanner scanner, BufferedReader buffreader,
			HhpredAlignment hhpredAlignment) {
		try {
			scanner = new Scanner(buffreader);
			int count = 1;
			final int RESNUMLINE = 4;
			String line = "";
			String resNumline = "";
			while (scanner.hasNextLine()) {
				line = scanner.nextLine();
				if (count == RESNUMLINE - 1) {
					resNumline = scanner.nextLine();
				}
				if (count > RESNUMLINE - 1 && scanner.hasNextLine()) {
					HhpredSequenceAlignmentColumn previousColumn;
					int numOfRow = count - RESNUMLINE + 1;
					if (this.size() != 0)
						previousColumn = this.get(this.size() - 1);
					else
						previousColumn = null;
					HhpredSequenceAlignmentColumn hsac = new HhpredSequenceAlignmentColumn(
							resNumline, line, numOfRow, hhpredAlignment,
							previousColumn);
					add(hsac);
				}

				count++;
			}
		}

		finally {
			if (scanner != null) {
				scanner.close();
			}
		}
	}

}