package meshi.parameters;

/**
 * Created by user on 15/02/2017.
 */
public enum LocalStructureAlphabetType {
    DSSP3("HCE"),
    DSSP7("HCEBGIT"),
    STR2("HCEBGITAPMZQ"),
    DSSP33("Reduction of all LocalStructure of DSSP to 30 top strings with the highmost frequences (70% of data) + HCE"),
    DSSP103("Reduction of all LocalStructure of DSSP to 100 top strings with the highmost frequences (90%) + HCE");

    private final String alphabetDescription;
    private LocalStructureAlphabetType(String alphabetDescription) {this.alphabetDescription = alphabetDescription;}

}
