/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.geometry.Coordinates;
import meshi.geometry.putH.NotEnoughBoundAtomsException;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.parameters.BBatom;
import meshi.parameters.ResidueMode;
import meshi.parameters.ResidueType;

public class ResidueExtendedAtoms extends Residue {
    public static final ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();

    public final Atom H, N, CA, CB, C, O, OXT;


    protected ResidueExtendedAtoms(String name) {
        super(name);
        H = N = CA = CB = C = O = OXT = null;
    }

    /**
     * Do not use this constructor to instantiate a creator object.
     */
    public ResidueExtendedAtoms(ResidueType type,
                                AtomList atomList,
                                ResidueIdentifier id,
                                ResidueMode mode, MolecularSystem molecularSystem) {
        super(type.nameThreeLetters(), type, id, new AtomList(molecularSystem), mode);
        Atom old;
        double tf;


        CA = getAtom("CA", type.caType(), atomList,this,molecularSystem);

        old = find(atomList, BBatom.N);
        if (old == null) tf = new Double(0);
        else tf = old.temperatureFactor();
        N = new Atom("N", this, type.nType(), new Coordinates(old), 1,tf,molecularSystem);
        if (mode == ResidueMode.NTER) {
            N.setType(AtomType.TRN);
            H = null;
        } else {
            if (type == ResidueType.PRO) H = null;
            else {
                old = find(atomList, BBatom.H);
                if (old == null) tf = new Double(0);
                else tf = old.temperatureFactor();
                H = new Atom("H", this, type.hType(), new Coordinates(old), 1,tf,molecularSystem);
            }
        }

        if (type != ResidueType.GLY) {
            old = find(atomList, BBatom.CB);
            if (old == null) tf = new Double(0);
            else tf = old.temperatureFactor();
            CB = new Atom("CB", this, type.cbType(), new Coordinates(old), 1,tf,molecularSystem);
        } else CB = null;

        old = find(atomList, BBatom.C);
        if (old == null) tf = new Double(0);
        else tf = old.temperatureFactor();
        C = new Atom("C", this, type.cType(), new Coordinates(old),1, tf,molecularSystem);

        old = find(atomList, BBatom.O);
        if (old == null) tf = new Double(0);
        else tf = old.temperatureFactor();
        O = new Atom("O", this, type.oType(), new Coordinates(old), 1,tf,molecularSystem);

        if (mode != ResidueMode.CTER) {
            OXT = null;
        } else {
            C.setType(AtomType.TRC);
            O.setType(AtomType.TRO);
            old = find(atomList, "OXT");
            if (old == null) tf = new Double(0);
            else tf = old.temperatureFactor();
            OXT = new Atom("OXT", this, AtomType.TRO, new Coordinates(old),1,  tf,molecularSystem);
        }
        addAtomsAndBonds();
    }


    private void addAtomsAndBonds() {
        if (H != null) getAtoms().add(H);
        getAtoms().add(N);
        getAtoms().add(CA);
        if (CB != null) getAtoms().add(CB);
        getAtoms().add(C);
        getAtoms().add(O);
        if (OXT != null) getAtoms().add(OXT);
        if (H != null) bonds.add(H.bond(N));
        bonds.add(N.bond(CA));
        bonds.add(CA.bond(C));
        bonds.add(C.bond(O));
        if (CB != null) bonds.add(CA.bond(CB));
        if (OXT != null) bonds.add(C.bond(OXT));
    }


    public static Atom getAtom(String name, AtomType type, AtomList atomList, Residue residue,
                               MolecularSystem molecularSystem) {
        Atom old = find(atomList, name);
        double tf;
        if (old == null) tf = new Double(0);
        else tf = old.temperatureFactor();
        return new Atom(name, residue, type, new Coordinates(old), 1, tf,molecularSystem);
    }

    public Atom head() {
        return C;
    }

    public Atom tail() {
        return N;
    }


    private static String nameConverter(String name) {
        if (name.equals("1HD2")) return "HD21";
        if (name.equals("2HD2")) return "HD22";
        if (name.equals("1HE2")) return "HE21";
        if (name.equals("2HE2")) return "HE22";
        return name;
    }

    public final Atom ca() {
        return CA;
    }

    public final Atom cb() {
        return CB;
    }

    public final Atom c() {
        return C;
    }

    public final Atom o() {
        return O;
    }

    public final Atom n() {
        return N;
    }

    public final Atom h() {
        return H;
    }

    public final Atom amideN() {
        return N;
    }

    // ***this is Guy's code***

    //a)
    /*      AA replacing methods        */

    /*public Ala toAla() {
        Arg out;
        ResidueMode mode = getMode();
        return new Ala(getgetAtoms(), ID(), mode);
    }

    public Cys toCys() {
        ResidueMode mode = getMode();
        return new Cys(getgetAtoms(), ID(), mode);
    }

    public Asp toAsp() {
        ResidueMode mode = getMode();
        return new Asp(getgetAtoms(), ID(), mode);
    }

    public Glu toGlu() {
        ResidueMode mode = getMode();
        return new Glu(getgetAtoms(), ID(), mode);
    }

    public His toHis() {
        ResidueMode mode = getMode();
        return new His(getgetAtoms(), ID(), mode);
    }

    public Ile toIle() {
        ResidueMode mode = getMode();
        return new Ile(getgetAtoms(), ID(), mode);
    }

    public Lys toLys() {
        ResidueMode mode = getMode();
        return new Lys(getgetAtoms(), ID(), mode);
    }

    public Leu toLeu() {
        ResidueMode mode = getMode();
        return new Leu(getgetAtoms(), ID(), mode);
    }

    public Met toMet() {
        ResidueMode mode = getMode();
        return new Met(getgetAtoms(), ID(), mode);
    }

    public Asn toAsn() {
        ResidueMode mode = getMode();
        return new Asn(getgetAtoms(), ID(), mode);
    }

    public Pro toPro() {
        ResidueMode mode = getMode();
        return new Pro(getgetAtoms(), ID(), mode);
    }

    public Gln toGln() {
        ResidueMode mode = getMode();
        return new Gln(getgetAtoms(), ID(), mode);
    }

    public Arg toArg() {
        ResidueMode mode = getMode();
        return new Arg(getgetAtoms(), ID(), mode);
    }

    public Ser toSer() {
        ResidueMode mode = getMode();
        return new Ser(getgetAtoms(), ID(), mode);
    }

    public Thr toThr() {
        ResidueMode mode = getMode();
        return new Thr(getgetAtoms(), ID(), mode);
    }

    public Val toVal() {
        ResidueMode mode = getMode();
        return new Val(getgetAtoms(), ID(), mode);
    }

    public Trp toTrp() {
        ResidueMode mode = getMode();
        return new Trp(getgetAtoms(), ID(), mode);
    }

    public Tyr toTyr() {
        ResidueMode mode = getMode();
        return new Tyr(getgetAtoms(), ID(), mode);
    }
      */
    public void addHydrogens(BondParametersList bondParametersList,
                             AngleParametersList angleParametersList) {
        if ((H != null) &&
                H.nowhere() &&
                (!N.nowhere()) &&
                (!CA.nowhere())) {
                PutHposLog log = PutHpos.pos(H, bondParametersList, angleParametersList);
                if ((log != null) && (log.numberOfLegandsWithCoordinatesOfHeavy == 1)) {/* The orientation of the hydrogen is random (as the
									coordinates of the carbonyl atom of the previous residue
									are not setResidue.*/
                    H.resetCoordinates();
                }
        }
    }

    public static void addHydrogens(Protein protein, BondParametersList bondParametersList,
                                    AngleParametersList angleParametersList) {
        for (Residue residue : protein.residues()) {
            ((ResidueExtendedAtoms) residue).addHydrogens(bondParametersList, angleParametersList);
        }
    }
}

