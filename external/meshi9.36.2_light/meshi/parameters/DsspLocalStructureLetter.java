package meshi.parameters;


/**
 * Created by user on 14/02/2017.
 */
public enum DsspLocalStructureLetter {

    HELIX("H"),
    SHEET("E"),
    BETABRIDGE("B"),
    THREE10HELIX("G"),
    PIHELIX("I"),
    BEND("S"),
    TURN("T"),
    HB_OandN("X"),
    HB_onlyO(">"),
    HB_onlyN("<"),
    BRACKET_3RES("3"),
    BRACKET_4RES("4"),
    BRACKET_5RES("5"),
    CHIRALITY_NEGATIVE("-"),
    CHIRALITY_POSITIVE("+"),
    PARALLEL("a"),
    ANTI_PARALLEL("A"),
    GAP("_"),
    UNK("!"),

    // Reduction letters - additional semantic letters
    //STR2 additional letters
    TWO_SIDED_PARALLEL_BETA("P"),
    TWO_SIDED_ANTI_PARALLEL_BETA("A"),
    TWO_SIDED_MIXED_BETA("M"),
    ONE_SIDED_PARALLEL_BETA("Q"),
    ONE_SIDED_ANTI_PARALLEL_BETA("Z"),

    //DSSP3 additional letters
    COIL("C");

    public String getNameOneLetter() {
        return nameOneLetter;
    }

    private final String nameOneLetter;
    private DsspLocalStructureLetter(String nameOneLetter) {
        this.nameOneLetter = nameOneLetter;
    }



    public static DsspLocalStructureLetter dsspLocalStructureLetter(int loc, char c) {
        char newC = Character.toUpperCase(c);
        DsspLocalStructureLetter dls = null;
        if (loc == 7 || loc == 8) {
            if (Character.isLowerCase(c)) dls = DsspLocalStructureLetter.PARALLEL;
            else if (Character.isUpperCase(c)) dls = DsspLocalStructureLetter.ANTI_PARALLEL;
            else dls = DsspLocalStructureLetter.GAP;
        } else { //if
            for (DsspLocalStructureLetter sec_letter : DsspLocalStructureLetter.values()) {
                if (sec_letter != DsspLocalStructureLetter.PARALLEL &&
                        sec_letter != DsspLocalStructureLetter.ANTI_PARALLEL &&
                        sec_letter.nameOneLetter.charAt(0) == newC)
                    dls = sec_letter;
            }//for
        }//else
        if ( dls!=null && dls.isLegalbyLocation(loc))
            return dls;
        else throw new RuntimeException("Undefined Dssp Local Structure Letter for location " +loc+": " + newC + "\n" + "Please take a look at meshi.parameters.DsspLocalStructureLetter.");
    }

    //retrun true if the letter can be on the first position of the dssp structure definition, else - false
    public boolean isLegalbyLocation(int loc){
        if (    (loc == 1 && isLoc1()) ||
                (loc == 2 && isLoc2_4()) || (loc == 3 && isLoc2_4()) || (loc == 4 && isLoc2_4()) ||
                (loc == 5 && isLoc5()) ||
                (loc ==6 && isLoc6()) ||
                (loc ==7 && isLoc7_8()) ||
                (loc ==8 && isLoc7_8())
                )
            return true;
        else return false;
    }
    public boolean isLoc1(){
        if (this == HELIX ||
                this == SHEET ||
                this == BETABRIDGE ||
                this == THREE10HELIX ||
                this == PIHELIX ||
                this == BEND ||
                this == TURN ||
                this == GAP)
            return true;
        else return false;
    }

    public boolean isLoc2_4(){
        if (this == HB_OandN ||
                this == HB_onlyO ||
                this == HB_onlyN ||
                this == BRACKET_3RES ||
                this == BRACKET_4RES ||
                this == BRACKET_5RES ||
                this == GAP) return true;
        else return false;
    }
    public boolean isLoc5(){
        if (this == BEND || this == GAP)
            return true;
        else return false;
    }
    public boolean isLoc6(){
        if (this == CHIRALITY_NEGATIVE ||
                this == CHIRALITY_POSITIVE ||
                this == GAP) return true;
        else return false;
    }
    public boolean isLoc7_8(){
        if (this == PARALLEL ||
                this == ANTI_PARALLEL ||
                this == GAP) return true;
        else return false;
    }


}
