/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

/**
 * A unique identifier for residues. This class deals
 * with the fact that proteins are not always comprised of a single
 * polypeptide chain and a residue's identity cannot be expressed
 * solely by its number in the polypeptide chain, as it was before.
 * <p/>
 * Field <code>number	</code> are determined
 * at creation and are final; field <code>chain</code> may be changed
 * later.
 *
 * @author Oren
 */
public class ResidueIdentifier implements Comparable {

    public final String chain;
    public final Integer number; // not 'int' so it could use Integer.equals/compareTo
    public final Integer chainNumber;


    public ResidueIdentifier(int number) {
        this(Chain.GENERIC_CHAIN_NAME, number, Chain.GENERIC_CHAIN_NUMBER);
    }

    /**
     * XmlPatchLibrary constructor of this class.
     */
    public ResidueIdentifier(String chain, int number, int chainNumber) {
        this.chain = chain;
        this.number = number;
        this.chainNumber = chainNumber;
    }


    /**
     * Compares <code>this</code> ResidueIdentifier to another object.
     * The order is lexicographic: first the chains are compared and
     * if they are equal the numbers are compared.
     *
     * @throws ClassCastException if the given object is not a
     *                            ResidueIdentifier.
     */
    public int compareTo(Object obj) {
        if (obj == null) throw new NullPointerException();

        ResidueIdentifier other = (ResidueIdentifier) obj;

        if (!chain.equals(other.chain))
            return chain.compareTo(other.chain);
        else
            return number.compareTo(other.number);
    }

    /**
     * Returns <code>true</code> if and only if <code>this</code>
     * ResidueIdentifier is equal to some other object.
     */
    public boolean equals(Object other) {
        if (!(other instanceof ResidueIdentifier))
            return false;
        else
            return (compareTo(other) == 0);
    }

    public final int number() {
        return number;
    }

    public final String chain() {return chain;}
}
