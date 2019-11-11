package meshi.applications.prediction.beautify;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.sequences.AlignmentColumn;
import meshi.sequences.ResidueAlignment;
import meshi.util.MeshiAttribute;

import java.util.Iterator;

public class BeautifyAttribute implements MeshiAttribute {
    private boolean loop;
    private boolean problematic;
    private boolean visible;
    private String problematicComment;
    private Residue problematicNeighbor;

    public BeautifyAttribute() {
        loop = problematic = visible = false;
        problematicComment = "";
        problematicNeighbor = null;
    }

    public void setLoop(boolean flag) {
        loop = flag;
        if (visible) visible = flag;
    }

    public boolean isLoop() {
        return loop;
    }

    public static void setLoop(Residue residue, boolean flag) {
        getBeautifyAttribute(residue).setLoop(flag);
    }

    public static boolean isLoop(Residue residue) {
        return getBeautifyAttribute(residue).isLoop();
    }


    public void setProblematic(boolean flag, String comment) {
        problematic = flag;
        problematicComment = comment;
        if (visible) visible = flag;
    }

    public boolean isProblematic() {
        return problematic;
    }

    public String problematicComment() {
        return problematicComment;
    }

    public static void setProblematic(Residue residue, boolean flag, String comment) {
        getBeautifyAttribute(residue).setProblematic(flag, comment);
    }

    public static boolean isProblematic(Residue residue) {
        return getBeautifyAttribute(residue).isProblematic();
    }

    public static String problematicComment(Residue residue) {
        return getBeautifyAttribute(residue).problematicComment();
    }

    public static boolean isVisible(Residue residue) {
        return getBeautifyAttribute(residue).visible;
    }

    public static void setVisible(Residue residue) {
        (getBeautifyAttribute(residue)).visible = true;
    }


    public void setProblematicNeighbor(Residue neighbor) {
        problematicNeighbor = neighbor;
    }

    public boolean hasProblematicNeighbor() {
        return problematicNeighbor != null;
    }

    public Residue problematicNeighbor() {
        return problematicNeighbor;
    }

    public static void setProblematicNeighbor(Residue residue, Residue neighbor) {
        getBeautifyAttribute(residue).setProblematicNeighbor(neighbor);
    }

    public static boolean hasProblematicNeighbor(Residue residue) {
        return getBeautifyAttribute(residue).hasProblematicNeighbor();
    }

    public static Residue problematicNeighbor(Residue residue) {
        return getBeautifyAttribute(residue).problematicNeighbor();
    }

    public static boolean isProblematicOrHasProblematicNeighbor(Residue residue) {
        return (isProblematic(residue) || hasProblematicNeighbor(residue));
    }

    public static boolean isProblematicOrLoop(Residue residue) {
        return (isProblematic(residue) || isLoop(residue));
    }

    public static void copy(Residue source, Residue target) {
        BeautifyAttribute baSource = getBeautifyAttribute(source);
        BeautifyAttribute baTarget = getBeautifyAttribute(target);
        baTarget.setLoop(baSource.isLoop());
        baTarget.setProblematic(baSource.isProblematic(), baSource.problematicComment());
        baTarget.setProblematicNeighbor(baSource.problematicNeighbor);
    }

    public static void copy(ResidueAlignment residueAlignment) {
        for (Iterator columns = residueAlignment.iterator(); columns.hasNext();) {
            AlignmentColumn column = (AlignmentColumn) columns.next();
            Residue modelResidue = (Residue) column.cell0().obj;
            Residue shotgunResidue = (Residue) column.cell1().obj;
            BeautifyAttribute ba = (BeautifyAttribute) shotgunResidue.getAttribute(MeshiAttribute.BEAUTIFY_ATTRIBUTE);
            if (ba != null) {
                modelResidue.addAttribute(ba);
            }
        }
    }

    public static void addTo(Protein protein) {
        for (Iterator residues = protein.residues().iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            residue.addAttribute(new BeautifyAttribute());
        }
    }

    public static BeautifyAttribute getBeautifyAttribute(Residue residue) {
        return (BeautifyAttribute) residue.getAttribute(BEAUTIFY_ATTRIBUTE);
    }

    public String toSrtring() {
        return "BeautifyAttribute " + loop + " " + problematic + " " + problematicComment + " " + problematicNeighbor + " " + visible;
    }

    public int key() {
        return BEAUTIFY_ATTRIBUTE;
    }
}
 
