package meshi.applications.HHpred;

import meshi.sequences.AlignmentCell;

public enum Condition {
	Q_RES_AND_TRES_NOT_DUMMY, Q_RES_DUMMY_AND_TRES_NOT_DUMMY, Q_RES_NOT_DUMMY_AND_TRES_DUMMY, Q_RES_AND_TRES_ARE_DUMMY;

	public static Condition getCondition(HhpredAlignmentColumn hhalcol) {
		Condition condition = Q_RES_AND_TRES_NOT_DUMMY;
		AlignmentCell tCell = hhalcol.getTCell();
		AlignmentCell qCell = hhalcol.getQCell();
		String qType = qCell.obj.toString();
		String tType = tCell.obj.toString();
		if (qType.equals("DMY") && !tType.equals("DMY"))
			condition = Condition.Q_RES_DUMMY_AND_TRES_NOT_DUMMY;
		else if (!qType.equals("DMY") && tType.equals("DMY"))
			condition = Condition.Q_RES_NOT_DUMMY_AND_TRES_DUMMY;
		return condition;
	}

}
