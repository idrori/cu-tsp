/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

import meshi.sequences.aligner.IdentityMatrix;
import meshi.sequences.aligner.SubstitutionMatrix;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.util.*;

import java.util.*;

public class ResidueAlignment extends ArrayList<ResidueAlignmentColumn> {
    private ResidueAlignmentColumn lastColumn;
    public final StringList comments;
    private SubstitutionMatrix substitutionMatrix;
    private double score = -1;
    private SequenceAlignment sequenceAlignment;

    

    /**
     * Empty alignment.
     */
    public ResidueAlignment() {
        super();
        lastColumn = null;
        comments = new StringList();
    }



    /**
     * A trivial alignment of two protein objects which are assumed to be two models of the same protein.
     */
    public ResidueAlignment(Chain chain1, String name1,//is this constructor needed?? Tommer 24.9.14
                            Chain chain2, String name2,
                            ResidueAlignmentMethod alignmentMethod) throws AlignmentException{
        this(new ChainList(chain1), name1, new ChainList(chain2), name2, alignmentMethod);

    }

        /**
         * A trivial alignment of two protein objects which are assumed to be two models of the same protein.
         */
    public ResidueAlignment(ChainList chains1, String name1,//is this constructor needed?? Tommer 24.9.14
                            ChainList chains2, String name2,
                            ResidueAlignmentMethod alignmentMethod) throws AlignmentException{
        this();
        if (chains1.size() != chains2.size())
            throw new RuntimeException("Current version can only handle proteins with the same number of chains");
        comments.add(name1);
        comments.add(name2);

        for (int iChain = 0; iChain < chains1.size(); iChain++) {
            Chain chain1 = chains1.get(iChain);
            Chain chain2 = chains2.get(iChain);

            MeshiSequence sequence1 = chain1.sequence();
            MeshiSequence sequence2 = chain2.sequence();

            if (alignmentMethod == ResidueAlignmentMethod.IDENTITY) {
                sequenceAlignment = SequenceAlignment.substitutionAlignment(sequence1, sequence2, new IdentityMatrix(),ResidueAlignmentMethod.IDENTITY);

                score = sequenceAlignment.score();
                for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext(); ) {
                    AlignmentColumn column = (AlignmentColumn) columns.next();
                    AlignmentCell cell1 = column.cell(0);
                    AlignmentCell cell2 = column.cell(1);
                    Residue residue1 = (Residue) cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
                    Residue residue2 = (Residue) cell2.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
                    if ((residue1 != null) & (residue2 != null)) {
                        ResidueAlignmentColumn newColumn = new ResidueAlignmentColumn(residue1, residue2);
                        add(newColumn);
                    }
                }
            } else if (alignmentMethod == ResidueAlignmentMethod.BY_RESIDUE_NUMBER) {
                int low, high, low1 = -9999, low2 = -8888, high1 = -77777, high2 = -66666;
                Residue residue1, residue2;
                for (Residue residue : chain1) {
                    if (!residue.dummy()) {
                        low1 = residue.number();
                        break;
                    }
                }
                high1 = chain1.get(chain1.size() - 1).number();
                for (Residue residue : chain2) {
                    if (!residue.dummy()) {
                        low2 = residue.number();
                        break;
                    }
                }
                high2 = chain2.get(chain2.size() - 1).number();

                if (low1 > low2) low = low1;
                else low = low2;
                if (high1 < high2) high = high1;
                else high = high2;
                for (int i = low; i <= high; i++) {
                    residue1 = chain1.get(i);
                    if (residue1.dummy()) continue;
                    residue2 = chain2.get(i);
                    if (residue2.dummy()) continue;
                    if ((residue1.name.equals(residue2.name)) ||
                            (residue1.name.equals("MET") && residue2.name.equals("MSE")) ||
                            (residue1.name.equals("MSE") && residue2.name.equals("MET"))) {
                        add(new ResidueAlignmentColumn(residue1, residue2));
                    } else throw new AlignmentException("Error while trying to align " + residue1 + "\t" + residue2);
                }
            }
        }
    }
    
    
/*	public ResidueAlignment(Chain chain1, String name1, Chain chain2,
			String name2, ResidueAlignmentMethod alignmentMethod, SubstitutionMatrix substitutionMatrix) //Added by Tommer 1.9.14
			throws AlignmentException {
		this();
		comments.add(name1);
		comments.add(name2);
		MeshiSequence sequence1 = chain1.sequence();
		MeshiSequence sequence2 = chain2.sequence();
		MeshiWriter output;
		this.substitutionMatrix = substitutionMatrix;//Added
       

		if (alignmentMethod == ResidueAlignmentMethod.IDENTITY) {
			SequenceAlignment sa = SequenceAlignment.substitutionAlignment(
					sequence1, sequence2, substitutionMatrix);// Po kavur
																		// haKelev,
																		// Tommer
																		// 31.8.14
			// signature changed Tommer 1.9.14
			//System.out.print(sa); // Added by Tommer 1.9.14
			for (Iterator columns = sa.iterator(); columns.hasNext();) {
				AlignmentColumn column = (AlignmentColumn) columns.next();
				AlignmentCell cell1 = column.cell(0);
				AlignmentCell cell2 = column.cell(1);
				Residue residue1 = (Residue) cell1
						.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
				Residue residue2 = (Residue) cell2
						.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
				if ((residue1 != null) & (residue2 != null)) {
					ResidueAlignmentColumn newColumn = new ResidueAlignmentColumn(
							residue1, residue2);
					add(newColumn);
				}
			}
		} else if (alignmentMethod == ResidueAlignmentMethod.BY_RESIDUE_NUMBER) {
			int low, high, low1 = -9999, low2 = -8888, high1 = -77777, high2 = -66666;
			Residue residue1, residue2;
			for (Residue residue : chain1) {
				if (!residue.dummy()) {
					low1 = residue.number();
					break;
				}
			}
			high1 = chain1.get(chain1.size() - 1).number();
			for (Residue residue : chain2) {
				if (!residue.dummy()) {
					low2 = residue.number();
					break;
				}
			}
			high2 = chain2.get(chain2.size() - 1).number();

			if (low1 > low2)
				low = low1;
			else
				low = low2;
			if (high1 < high2)
				high = high1;
			else
				high = high2;
			for (int i = low; i <= high; i++) {
				residue1 = chain1.get(i);
				if (residue1.dummy())
					continue;
				residue2 = chain2.get(i);
				if (residue2.dummy())
					continue;
				if ((residue1.name.equals(residue2.name))
						|| (residue1.name.equals("MET") && residue2.name
								.equals("MSE"))
						|| (residue1.name.equals("MSE") && residue2.name
								.equals("MET"))) {
					add(new ResidueAlignmentColumn(residue1, residue2));
				} else
					throw new AlignmentException("Error while trying to align "
							+ residue1 + "\t" + residue2);
			}
		}
	}
  */


    public ResidueAlignment(ResidueSequence residueSequence0, ResidueSequence residueSequence1) {
        this(SequenceAlignment.substitutionAlignment(residueSequence0, residueSequence1, new IdentityMatrix(), ResidueAlignmentMethod.IDENTITY));//Tommer 2.9.14
    }

    public ResidueAlignment(ResidueSequence residueSequence0, ResidueSequence residueSequence1, SubstitutionMatrix substitutionMatrix) {
        this(SequenceAlignment.substitutionAlignment(residueSequence0, residueSequence1, substitutionMatrix, ResidueAlignmentMethod.IDENTITY));//Tommer 2.9.14
    }

    public ResidueAlignment(Chain chain0, Chain chain1, SequenceList sequenceList) {
        this(chain0,chain1,sequenceList, new IdentityMatrix());
    }
    public ResidueAlignment(Chain chain0, Chain chain1, SequenceList sequenceList, SubstitutionMatrix substitutionMatrix) {
        this();
        AlignmentColumn column;
        AlignmentCell from;
        AlignmentCell to;
        Residue residue, residue0, residue1;
        AlignmentCell cell0, cell1;

        SequenceList newSequenceList = new SequenceList(sequenceList);
        MeshiSequence sequence0 = chain0.sequence();
        MeshiSequence sequence1 = chain1.sequence();
        MeshiSequence sequence0FromSL = newSequenceList.get(0);
        MeshiSequence sequence1FromSL = newSequenceList.get(1);
        SequenceAlignment sequence0Alignment = SequenceAlignment.substitutionAlignment(sequence0, sequence0FromSL, substitutionMatrix, ResidueAlignmentMethod.IDENTITY);
        SequenceAlignment sequence1Alignment = SequenceAlignment.substitutionAlignment(sequence1, sequence1FromSL, substitutionMatrix, ResidueAlignmentMethod.IDENTITY);

        for (Iterator columns = sequence0Alignment.iterator(); columns.hasNext();) {
            column = (AlignmentColumn) columns.next();
            from = column.cell(0);
            to = column.cell(1);
            residue = (Residue) from.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null) to.addAttribute(residue);
        }
        for (Iterator columns = sequence1Alignment.iterator(); columns.hasNext();) {
            column = (SequenceAlignmentColumn) columns.next();
            from = column.cell(0);
            to = column.cell(1);
            residue = (Residue) from.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null) to.addAttribute(residue);
        }

        SequenceAlignment sequenceAlignment = new SequenceAlignment(newSequenceList);
        for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext();) {
            column = (SequenceAlignmentColumn) columns.next();
            cell0 = column.cell(0);
            cell1 = column.cell(1);
            residue0 = (Residue) cell0.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            residue1 = (Residue) cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if ((residue0 != null) & (residue1 != null) & (!SequenceAlignmentCell.wildcardOrGap((SequenceAlignmentColumn) column))) {
                add(new ResidueAlignmentColumn(residue0, residue1));
            }
        }
    }


    public ResidueAlignment(SequenceAlignment sequenceAlignment) {
        this();
        int key = MeshiAttribute.RESIDUE_ATTRIBUTE;
        for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext();) {
            SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
            SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
            SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
            Residue residue0 = (Residue) cell0.getAttribute(key);
            Residue residue1 = (Residue) cell1.getAttribute(key);
            if ((residue0 != null) & (residue1 != null) & (!SequenceAlignmentCell.wildcardOrGap(column))) {
                ResidueAlignmentCell residueCell0 = new ResidueAlignmentCell(residue0);
                ResidueAlignmentCell residueCell1 = new ResidueAlignmentCell(residue1);
                ResidueAlignmentColumn residueColumn = new ResidueAlignmentColumn(residueCell0, residueCell1);
                add(residueColumn);
            }
        }
    }

    public SequenceAlignment getSequenceAlignment() {
        return sequenceAlignment;
    }

    public boolean add(ResidueAlignmentColumn column) {
        if (lastColumn == null) {
            lastColumn = column;
            return super.add(column);
        }
        if (!column.cell0().gap())
            if (column.residue0().ID().compareTo(lastColumn.residue0().ID()) < 0) return false;
        if (!column.cell1().gap())
            if (column.residue1().ID().compareTo(lastColumn.residue1().ID()) < 0) return false;
        if ((!column.cell0().gap()) &
                (!column.cell1().gap())) lastColumn = column;
        return super.add(column);
    }

    public ResidueAlignment insert(ResidueAlignmentColumn column) {
        ResidueAlignment out = new ResidueAlignment();
        for (String c : comments)
            out.comments.add(c);
        boolean added = false;
        for (Iterator columns = iterator(); columns.hasNext();) {
            ResidueAlignmentColumn current = (ResidueAlignmentColumn) columns.next();
            if (!added) {
                if ((current.residue0().number() < column.residue0().number()) &
                        (current.residue1().number() < column.residue1().number())) {
                    if (!out.add(current)) throw new RuntimeException("weird situation 1");
                } else if ((current.residue0().number() > column.residue0().number()) &
                        (current.residue1().number() > column.residue1().number())) {
                    if (!out.add(column)) throw new RuntimeException("weird situation 2");
                    if (!out.add(current)) throw new RuntimeException("weird situation 3");
                    added = true;
                } else {
                    return (null);
                }
            } else {
                if ((current.residue0().number() <= column.residue0().number()) |
                        (current.residue1().number() <= column.residue1().number())) {
                    throw new RuntimeException("weird situation 4");
                }
                if (!out.add(current)) throw new RuntimeException("weird situation 5");
            }
        }
        if (!added) if (!out.add(column)) throw new RuntimeException("weird situation 6");
        return out;
    }


    public AtomList getCaList(int row) {
        return (new AtomAlignment(toAlignmentArray(this), new CaFilter())).atomList(row);
    }

    public double score() {return score;}
    public String toString() {
        //if (1 == 1) throw new RuntimeException("xxxxxxxxxxxxx");
        return (new SequenceAlignment(this)).toString();
    }


    public ResidueAlignment filter(Filter filter) {
        ResidueAlignment out = new ResidueAlignment();
        for (ResidueAlignmentColumn rc : this)
            if (filter.accept(rc)) out.add(rc);
        return out;
    }

    public void print() {
        for (ResidueAlignmentColumn column : this)
            System.out.println(column);
    }

    public void print(MeshiWriter writer) {
        for (ResidueAlignmentColumn column : this)
            writer.println(column);
    }


    public int countIdentities() {
        int out = 0;
        for (ResidueAlignmentColumn column : this) {
            ResidueAlignmentCell cell0 = (ResidueAlignmentCell) column.cell0();
            ResidueAlignmentCell cell1 = (ResidueAlignmentCell) column.cell1();
            if (cell0.residue().type == cell1.residue().type)
                out++;
        }
        return out;
    }

    public static ResidueAlignment[] toAlignmentArray(ResidueAlignment  residueAlignment) {
        ResidueAlignment[] residueAlignments = new ResidueAlignment[1];
        residueAlignments[0] = residueAlignment;
        return residueAlignments;
    }

}
