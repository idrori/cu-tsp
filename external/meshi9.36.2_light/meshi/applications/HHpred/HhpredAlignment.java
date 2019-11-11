package meshi.applications.HHpred;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;

public class HhpredAlignment extends ArrayList<HhpredAlignmentColumn> {

	public HhpredAlignment() {
		super();
		// TODO Auto-generated constructor stub
	}

	public HhpredAlignment(File file) throws FileNotFoundException {
		new HhpredSequenceAlignment(file, this);
	}

	public void addColumn(HhpredAlignmentColumn halc) {
		add(halc);
	}

	@Override
	public String toString() {
		return "HhpredAlignment [\n" + super.toString() + "]";
	}
}