/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation.hydrogenBonds;

import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.util.Utils;
import meshi.util.mathTools.Sigma;

/**
 * A specific hydrogen bond description following: Dahiyat et al. Protein Sci. 1996 May;5(5):895-903.
 * <p/>
 * The hydrogen bond value has the following dependence:
 * hbVal = sigmoidDis(Rd-a)*sigmoidAng(angle1)*sigmoidAng(angle2)
 * <p/>
 * Where Rd-a is the distance between the donnor and acceptor, and sigmoidDis(Rd-a) is the sidmoid function
 * described in meshi.util.mathTools.Sigma.java  . The parameters for this sigmoid are given in the
 * constructor (p1,p2,end,valAtp1,valAtp2).
 * <p/>
 * The two angles in the HB are made of 4 getAtoms are defined as follows:
 * <p/>
 * a1~~~~~~a2
 * .
 * .
 * a3~~~~~~a4
 * <p/>
 * ...  -  the HB
 * ~~~  -  covelant bonding
 * Angle 1 - between getAtoms 1,2,3.
 * Angle 2 - between getAtoms 2,3,4.
 * There could be two cases:
 * 1) The hydrogen atom in the HB is not known explicitly (NoH). In this case a2 and a3 are the polar getAtoms,
 * and a1 and a4 are the heavy getAtoms to which they are attached (base getAtoms).
 * 2) The hydrogen atom in the HB is known explicitly (WithH). Either a2 or a3 must be a hydrogen. The other atom
 * COVELANTLY bonded to the hydrogen must be the polar atom.
 * <p/>
 * The angular score of each angle is a sigmoidAng that equals 0.0 bellow a certain threshold (sigmoidBegins) and 1.0
 * above a certain threshold (sigmoidEnds). Between the thresholds it raises smoothly. The HB angular score is the
 * product of the two sigmoidAns's of angle 1 and 2.
 * The sigmoid thresholds are obtained as constants from the interface "DahiyatImplementationConstants":
 * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when the hydrogen
 * in the HB is defined in MESHI.
 * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB sigmoids where the hydrogen in the HB is present.
 * <p/>
 * <p/>
 * Note:
 * 1) The constructor assume that the atom list it get as a parameter has the donor and acceptor in the first
 * two places (see documentation of mother class). The third place in the list is the atom that is covalently bonded
 * to the first atom. The forth place in the list is the atom that is covalently bonded to the second atom.
 * 2) The sidmoids of angle1 and 2 has the same form.
 */
public class HydrogenBondDahiyat extends AbstractHydrogenBond implements DahiyatImplementationConstants {

    // Parameters
    private double p1;
    private double p2;
    private double end;
    private double valAtp1;
    private double valAtp2;
    private double angSigmoidBegins; // The cosine of the begining
    private double angSigmoidEnds; // The cosine of the end

    // Variables used to calculate the angle sigmoids product
    // The zeroDerivative field is 'true' if the angles are in their smooth regimn, and this happens most of
    // the time. Some calculations can be skiped by using this field.
    private final Sigma sigma = new Sigma();
    private double hbAngScore;
    private boolean zeroDerivative;
    private double DhbAngScoreDx1;
    private double DhbAngScoreDy1;
    private double DhbAngScoreDz1;
    private double DhbAngScoreDx2;
    private double DhbAngScoreDy2;
    private double DhbAngScoreDz2;
    private double DhbAngScoreDx3;
    private double DhbAngScoreDy3;
    private double DhbAngScoreDz3;
    private double DhbAngScoreDx4;
    private double DhbAngScoreDy4;
    private double DhbAngScoreDz4;
    private double sigmCosAng1;
    private boolean derivativeAng1Zero;
    private double DsigmCosAng1Dx1;
    private double DsigmCosAng1Dy1;
    private double DsigmCosAng1Dz1;
    private double DsigmCosAng1Dx2;
    private double DsigmCosAng1Dy2;
    private double DsigmCosAng1Dz2;
    private double DsigmCosAng1Dx3;
    private double DsigmCosAng1Dy3;
    private double DsigmCosAng1Dz3;
    private double sigmCosAng2;
    private boolean derivativeAng2Zero;
    private double DsigmCosAng2Dx2;
    private double DsigmCosAng2Dy2;
    private double DsigmCosAng2Dz2;
    private double DsigmCosAng2Dx3;
    private double DsigmCosAng2Dy3;
    private double DsigmCosAng2Dz3;
    private double DsigmCosAng2Dx4;
    private double DsigmCosAng2Dy4;
    private double DsigmCosAng2Dz4;

    // Variables used by the distance sigmoid
    private double hbDisScore;
    private double DhbDisScoreDx1;
    private double DhbDisScoreDy1;
    private double DhbDisScoreDz1;
    private double DhbDisScoreDx2;
    private double DhbDisScoreDy2;
    private double DhbDisScoreDz2;
    private double DhbDisScoreDx3;
    private double DhbDisScoreDy3;
    private double DhbDisScoreDz3;
    private double DhbDisScoreDx4;
    private double DhbDisScoreDy4;
    private double DhbDisScoreDz4;

    // These distances are used in the calculations of the angles and distance sigmoids.
    // The numbering of the getAtoms (1 to 4) is as describe in the top of the class.
    Distance disAng1atoms12;
    Distance disAng1atoms32;
    Distance disAng2atoms23;
    Distance disAng2atoms43;
    Distance disDonorAcceptor;
    Distance hbDistance;


    // The 4 getAtoms in atom list a referenced by these fields, according to the numbering of the
    // getAtoms (1 to 4) that is described at the top of the class.
    private Atom a1;
    private Atom a2;
    private Atom a3;
    private Atom a4;

    // These fields tell how to convert from the 1..4 atom numbering to the initial atom list numbering
    private int a1toAtomList;
    private int a2toAtomList;
    private int a3toAtomList;
    private int a4toAtomList;


    public HydrogenBondDahiyat() {
        throw new RuntimeException("\nERROR: without parameters the hydrogen bonds cannot be formed.\n");
    }

    public HydrogenBondDahiyat(DistanceMatrix dm, AtomList atomList,
                               double p1, double p2, double end, double valAtp1, double valAtp2) {
        super(dm, atomList);
        this.p1 = p1;
        this.p2 = p2;
        this.end = end;
        this.valAtp1 = valAtp1;
        this.valAtp2 = valAtp2;

        // Finding the conversion between the atomList order and the 1..4 notation of getAtoms at the top.
        atomListTo1_4numbering();

        if (a2.type().isHydrogen() || a3.type().isHydrogen()) {// There is explicit hydrogen here
            angSigmoidBegins = Math.cos(sigmoidEndsWithH * Math.PI / 180.0); // note the switch end-begin. cos is a monotonicly DECREASING function
            angSigmoidEnds = Math.cos(sigmoidBeginsWithH * Math.PI / 180.0);
        } else {  // No explicit hydrogen
            angSigmoidBegins = Math.cos(sigmoidEndsNoH * Math.PI / 180.0); // note the switch end-begin. cos is a monotonicly DECREASING function
            angSigmoidEnds = Math.cos(sigmoidBeginsNoH * Math.PI / 180.0);
        }

        disAng1atoms12 = null;
        disAng1atoms32 = null;
        disAng2atoms23 = null;
        disAng2atoms43 = null;
        disDonorAcceptor = null;
        updateHBvalueAndDerivatives();
    }

    /**
     * This method activates all the code for updating the Dahiyat H-bond.
     */
    public void updateHBvalueAndDerivatives() {
        // Assigning distance variables
        Distance temp;
        if ((disAng1atoms12 == null) || disAng1atoms12.dead()) {
            disAng1atoms12 = dm.distance(a1, a2);
            if (disAng1atoms12 == null) disAng1atoms12 = new DistanceMirror(dm.distance(a2, a1));
        }

        if ((disAng1atoms32 == null) || disAng1atoms32.dead()) {
            disAng1atoms32 = dm.distance(a3, a2);
            if (disAng1atoms32 == null) {
                temp = dm.distance(a2, a3);
                if (temp == null) return;
                disAng1atoms32 = new DistanceMirror(temp);
            }
        }



        if ((disAng2atoms23 == null) || disAng2atoms23.dead()) {
            disAng2atoms23 = dm.distance(a2, a3);
            if (disAng2atoms23 == null) disAng2atoms23 = new DistanceMirror(dm.distance(a3, a2));
        }

//        if (disAng2atoms23.distance() < 2) {
//            Atom atom1 = disAng2atoms23.atom1();
//            Atom atom2 = disAng2atoms23.atom2();
////            if ((atom1.type() == AtomType.KH)  || (atom2.type() == AtomType.KH))
////                   meshi.util.Utils.printDebug(this,"????? "+disAng2atoms23);
//        }


        if ((disAng2atoms43 == null) || disAng2atoms43.dead()) {
            disAng2atoms43 = dm.distance(a4, a3);
            if (disAng2atoms43 == null) disAng2atoms43 = new DistanceMirror(dm.distance(a3, a4));
        }

        if ((disAng1atoms32.mode() != DistanceMode.MIRROR))
            hbDistance = disAng1atoms32;
        else hbDistance = disAng2atoms23;
        if (hbDistance.mode() == DistanceMode.INFINITE) {
            returnZero();
            return;
        }
        if ((disDonorAcceptor == null) || disDonorAcceptor.dead() ) {
            disDonorAcceptor = dm.distance(atomList.atomAt(0), atomList.atomAt(1));
            if (disDonorAcceptor == null) {
                temp = dm.distance(atomList.atomAt(0), atomList.atomAt(1));
                if (temp == null || temp.mode() == DistanceMode.INFINITE) {
                    returnZero();
                    return;
                }
                try {
                    disDonorAcceptor = new DistanceMirror(temp);
                }
                catch (RuntimeException ex) {
                    throw new RuntimeException("Cannot create a DistanceMirror between:\n" +
                            atomList.atomAt(0) + " and \n" + atomList.atomAt(1) + "\n" +
                            "their real distance is: " + (new FreeDistance(atomList.atomAt(0), atomList.atomAt(1))).distance() + "\n" +
                            "and that coased " + ex);
                }
            }
        }
        updateAngProduct1234();
        updateDis();
        hbVal = hbDisScore * hbAngScore;
        derivatives[a1toAtomList][0] = hbDisScore * DhbAngScoreDx1 + hbAngScore * DhbDisScoreDx1;
        derivatives[a1toAtomList][1] = hbDisScore * DhbAngScoreDy1 + hbAngScore * DhbDisScoreDy1;
        derivatives[a1toAtomList][2] = hbDisScore * DhbAngScoreDz1 + hbAngScore * DhbDisScoreDz1;

        derivatives[a2toAtomList][0] = hbDisScore * DhbAngScoreDx2 + hbAngScore * DhbDisScoreDx2;
        derivatives[a2toAtomList][1] = hbDisScore * DhbAngScoreDy2 + hbAngScore * DhbDisScoreDy2;
        derivatives[a2toAtomList][2] = hbDisScore * DhbAngScoreDz2 + hbAngScore * DhbDisScoreDz2;

        derivatives[a3toAtomList][0] = hbDisScore * DhbAngScoreDx3 + hbAngScore * DhbDisScoreDx3;
        derivatives[a3toAtomList][1] = hbDisScore * DhbAngScoreDy3 + hbAngScore * DhbDisScoreDy3;
        derivatives[a3toAtomList][2] = hbDisScore * DhbAngScoreDz3 + hbAngScore * DhbDisScoreDz3;

        derivatives[a4toAtomList][0] = hbDisScore * DhbAngScoreDx4 + hbAngScore * DhbDisScoreDx4;
        derivatives[a4toAtomList][1] = hbDisScore * DhbAngScoreDy4 + hbAngScore * DhbDisScoreDy4;
        derivatives[a4toAtomList][2] = hbDisScore * DhbAngScoreDz4 + hbAngScore * DhbDisScoreDz4;
    }

    public void returnZero() {
            hbVal = 0;
            derivatives[a1toAtomList][0] = 0;
            derivatives[a1toAtomList][1] = 0;
            derivatives[a1toAtomList][2] = 0;
            derivatives[a2toAtomList][0] = 0;
            derivatives[a2toAtomList][1] = 0;
            derivatives[a2toAtomList][2] = 0;
            derivatives[a3toAtomList][0] = 0;
            derivatives[a3toAtomList][1] = 0;
            derivatives[a3toAtomList][2] = 0;
            derivatives[a4toAtomList][0] = 0;
            derivatives[a4toAtomList][1] = 0;
            derivatives[a4toAtomList][2] = 0;
            return;
    }

    public void debug() {
        Utils.printDebug(this,"disDonorAcceptor = "+disDonorAcceptor);
        Utils.printDebug(this,"disDonorAcceptor.distance() = "+disDonorAcceptor.distance());
        Utils.printDebug(this," end, p1, p2, valAtp1, valAtp2 = "+end+" "+p1+" "+ p2+" "+ valAtp1+" "+ valAtp2);
        sigma.sigma(disDonorAcceptor.distance(), end, p1, p2, valAtp1, valAtp2);
        Utils.printDebug(this,"sigma.s() = "+sigma.s());
        Utils.printDebug(this, "sigmCosAng1 = "+sigmCosAng1+" sigmCosAng2 = "+sigmCosAng2);
        double d12x = disAng1atoms12.dDistanceDx();
        double d12y = disAng1atoms12.dDistanceDy();
        double d12z = disAng1atoms12.dDistanceDz();
        double d32x = disAng1atoms32.dDistanceDx();
        double d32y = disAng1atoms32.dDistanceDy();
        double d32z = disAng1atoms32.dDistanceDz();
        Utils.printDebug(this,"distances12&32 "+disAng1atoms12+" ; "+disAng1atoms32);
        Utils.printDebug(this,"Atoms "+disAng1atoms12.atom1()+" "+disAng1atoms12.atom2());
        Utils.printDebug(this,"ddddd  "+d12x+" "+d12y+" "+d12z+" "+d32x+" "+d32y+" "+d32z);
        double  cosAng = d12x*d32x+d12y*d32y+d12z*d32z;
                Utils.printDebug(this, "\ncosAng = "+ cosAng+" "+Math.acos(cosAng) +"\n");

    }
    /**
     * This method update the distance sigmoid and all the relevent derivatives.
     * The numbering of the getAtoms (1 to 4) is as describe in the top of the class.
     */
    private final void updateDis() {
        sigma.sigma(disDonorAcceptor.distance(), end, p1, p2, valAtp1, valAtp2);
        hbDisScore = sigma.s();
        if (!a2.type().isHydrogen() && !a3.type().isHydrogen()) {  // case 1:   base---O...O---base    or    base---O...N---base
            DhbDisScoreDx1 = 0;
            DhbDisScoreDy1 = 0;
            DhbDisScoreDz1 = 0;
            DhbDisScoreDx2 = sigma.s_tag() * disDonorAcceptor.dDistanceDx();
            DhbDisScoreDy2 = sigma.s_tag() * disDonorAcceptor.dDistanceDy();
            DhbDisScoreDz2 = sigma.s_tag() * disDonorAcceptor.dDistanceDz();
            DhbDisScoreDx3 = -DhbDisScoreDx2;
            DhbDisScoreDy3 = -DhbDisScoreDy2;
            DhbDisScoreDz3 = -DhbDisScoreDz2;
            DhbDisScoreDx4 = 0;
            DhbDisScoreDy4 = 0;
            DhbDisScoreDz4 = 0;
        } else if (!a2.type().isHydrogen() && a3.type().isHydrogen()) { // case 2:   base---O...H---N
            DhbDisScoreDx1 = 0;
            DhbDisScoreDy1 = 0;
            DhbDisScoreDz1 = 0;
            DhbDisScoreDx2 = sigma.s_tag() * disDonorAcceptor.dDistanceDx();
            DhbDisScoreDy2 = sigma.s_tag() * disDonorAcceptor.dDistanceDy();
            DhbDisScoreDz2 = sigma.s_tag() * disDonorAcceptor.dDistanceDz();
            DhbDisScoreDx3 = 0;
            DhbDisScoreDy3 = 0;
            DhbDisScoreDz3 = 0;
            DhbDisScoreDx4 = -DhbDisScoreDx2;
            DhbDisScoreDy4 = -DhbDisScoreDy2;
            DhbDisScoreDz4 = -DhbDisScoreDz2;
        } else if (a2.type().isHydrogen() && !a3.type().isHydrogen()) { // case 3:   N---H...O---base
            DhbDisScoreDx1 = sigma.s_tag() * disDonorAcceptor.dDistanceDx();
            DhbDisScoreDy1 = sigma.s_tag() * disDonorAcceptor.dDistanceDy();
            DhbDisScoreDz1 = sigma.s_tag() * disDonorAcceptor.dDistanceDz();
            DhbDisScoreDx2 = 0;
            DhbDisScoreDy2 = 0;
            DhbDisScoreDz2 = 0;
            DhbDisScoreDx3 = -DhbDisScoreDx1;
            DhbDisScoreDy3 = -DhbDisScoreDy1;
            DhbDisScoreDz3 = -DhbDisScoreDz1;
            DhbDisScoreDx4 = 0;
            DhbDisScoreDy4 = 0;
            DhbDisScoreDz4 = 0;
        } else if (a2.type().isHydrogen() && a3.type().isHydrogen()) { // case 4:   N---H...H---N
            DhbDisScoreDx1 = sigma.s_tag() * disDonorAcceptor.dDistanceDx();
            DhbDisScoreDy1 = sigma.s_tag() * disDonorAcceptor.dDistanceDy();
            DhbDisScoreDz1 = sigma.s_tag() * disDonorAcceptor.dDistanceDz();
            DhbDisScoreDx2 = 0;
            DhbDisScoreDy2 = 0;
            DhbDisScoreDz2 = 0;
            DhbDisScoreDx3 = 0;
            DhbDisScoreDy3 = 0;
            DhbDisScoreDz3 = 0;
            DhbDisScoreDx4 = -DhbDisScoreDx1;
            DhbDisScoreDy4 = -DhbDisScoreDy1;
            DhbDisScoreDz4 = -DhbDisScoreDz1;
        } else
            throw new RuntimeException("If this line is reached, something is definitely wrong");
    }


    /**
     * This method update the product: sigmoidAng(angle1)*sigmoidAng(angle2)
     * and all the relevent derivatives. The numbering of the getAtoms (1 to 4) is as describe in the top
     * of the class.
     */
    private final void updateAngProduct1234() {

        updateAng123();
        updateAng234();

        hbAngScore = sigmCosAng1 * sigmCosAng2;
        if (derivativeAng1Zero && derivativeAng2Zero) {
            zeroDerivative = true;
            DhbAngScoreDx1 = DhbAngScoreDy1 = DhbAngScoreDz1 =
                    DhbAngScoreDx2 = DhbAngScoreDy2 = DhbAngScoreDz2 =
                            DhbAngScoreDx3 = DhbAngScoreDy3 = DhbAngScoreDz3 =
                                    DhbAngScoreDx4 = DhbAngScoreDy4 = DhbAngScoreDz4 = 0.0;
        } else {
            zeroDerivative = false;
            DhbAngScoreDx1 = DsigmCosAng1Dx1 * sigmCosAng2;
            DhbAngScoreDy1 = DsigmCosAng1Dy1 * sigmCosAng2;
            DhbAngScoreDz1 = DsigmCosAng1Dz1 * sigmCosAng2;

            DhbAngScoreDx2 = DsigmCosAng1Dx2 * sigmCosAng2 + DsigmCosAng2Dx2 * sigmCosAng1;
            DhbAngScoreDy2 = DsigmCosAng1Dy2 * sigmCosAng2 + DsigmCosAng2Dy2 * sigmCosAng1;
            DhbAngScoreDz2 = DsigmCosAng1Dz2 * sigmCosAng2 + DsigmCosAng2Dz2 * sigmCosAng1;
            DhbAngScoreDx3 = DsigmCosAng1Dx3 * sigmCosAng2 + DsigmCosAng2Dx3 * sigmCosAng1;
            DhbAngScoreDy3 = DsigmCosAng1Dy3 * sigmCosAng2 + DsigmCosAng2Dy3 * sigmCosAng1;

            DhbAngScoreDz3 = DsigmCosAng1Dz3 * sigmCosAng2 + DsigmCosAng2Dz3 * sigmCosAng1;
            DhbAngScoreDx4 = DsigmCosAng2Dx4 * sigmCosAng1;
            DhbAngScoreDy4 = DsigmCosAng2Dy4 * sigmCosAng1;
            DhbAngScoreDz4 = DsigmCosAng2Dz4 * sigmCosAng1;
        }
    }

    // Auxilary method to 'updateAngProduct1234'

    private final void updateAng123() {
        double cosAng = disAng1atoms12.dDistanceDx() * disAng1atoms32.dDistanceDx() +
                disAng1atoms12.dDistanceDy() * disAng1atoms32.dDistanceDy() +
                disAng1atoms12.dDistanceDz() * disAng1atoms32.dDistanceDz();

        sigma.sigma(1.0 + cosAng, 2, 1.0 + angSigmoidBegins, 1.0 + angSigmoidEnds, 1.0, 0.0); //Chen - take a look
        sigmCosAng1 = sigma.s();
        if (sigma.s_tag() != 0.0) {
            derivativeAng1Zero = false;
            DsigmCosAng1Dx1 = sigma.s_tag() * disAng1atoms12.invDistance() * (disAng1atoms32.dDistanceDx() - cosAng * disAng1atoms12.dDistanceDx());
            DsigmCosAng1Dy1 = sigma.s_tag() * disAng1atoms12.invDistance() * (disAng1atoms32.dDistanceDy() - cosAng * disAng1atoms12.dDistanceDy());
            DsigmCosAng1Dz1 = sigma.s_tag() * disAng1atoms12.invDistance() * (disAng1atoms32.dDistanceDz() - cosAng * disAng1atoms12.dDistanceDz());
            DsigmCosAng1Dx3 = sigma.s_tag() * disAng1atoms32.invDistance() * (disAng1atoms12.dDistanceDx() - cosAng * disAng1atoms32.dDistanceDx());
            DsigmCosAng1Dy3 = sigma.s_tag() * disAng1atoms32.invDistance() * (disAng1atoms12.dDistanceDy() - cosAng * disAng1atoms32.dDistanceDy());
            DsigmCosAng1Dz3 = sigma.s_tag() * disAng1atoms32.invDistance() * (disAng1atoms12.dDistanceDz() - cosAng * disAng1atoms32.dDistanceDz());
            DsigmCosAng1Dx2 = -DsigmCosAng1Dx1 - DsigmCosAng1Dx3;
            DsigmCosAng1Dy2 = -DsigmCosAng1Dy1 - DsigmCosAng1Dy3;
            DsigmCosAng1Dz2 = -DsigmCosAng1Dz1 - DsigmCosAng1Dz3;
        } else {
            derivativeAng1Zero = true;
            DsigmCosAng1Dx1 = DsigmCosAng1Dy1 = DsigmCosAng1Dz1 = DsigmCosAng1Dx3 = DsigmCosAng1Dy3 =
                    DsigmCosAng1Dz3 = DsigmCosAng1Dx2 = DsigmCosAng1Dy2 = DsigmCosAng1Dz2 = 0.0;
        }
    }

    // Auxilary method to 'updateAngProduct1234'

    private final void updateAng234() {

        double cosAng = disAng2atoms23.dDistanceDx() * disAng2atoms43.dDistanceDx() +
                        disAng2atoms23.dDistanceDy() * disAng2atoms43.dDistanceDy() +
                        disAng2atoms23.dDistanceDz() * disAng2atoms43.dDistanceDz();

        if (cosAng < -1) {
            throw new RuntimeException("cosAng < -1  "+cosAng+"\nthis is " + this + "\n"+
                    disAng2atoms23+ "  "+ disAng2atoms23.mode()+"  "+(disAng2atoms23==hbDistance) +" "+hbDistance.mode()+"\n"+disAng2atoms23.atom1()+"\n"+disAng2atoms23.atom2()+"\n"+disAng2atoms43+"\n"+disDonorAcceptor+"\n"+disDonorAcceptor.atom1()+"\n"+disDonorAcceptor.atom2());
        }
        sigma.sigma(1.0 + cosAng, 2, 1.0 + angSigmoidBegins, 1.0 + angSigmoidEnds, 1.0, 0.0);
        sigmCosAng2 = sigma.s();
        if (sigma.s_tag() != 0.0) {
            derivativeAng2Zero = false;
            DsigmCosAng2Dx2 = sigma.s_tag() * disAng2atoms23.invDistance() * (disAng2atoms43.dDistanceDx() - cosAng * disAng2atoms23.dDistanceDx());
            DsigmCosAng2Dy2 = sigma.s_tag() * disAng2atoms23.invDistance() * (disAng2atoms43.dDistanceDy() - cosAng * disAng2atoms23.dDistanceDy());
            DsigmCosAng2Dz2 = sigma.s_tag() * disAng2atoms23.invDistance() * (disAng2atoms43.dDistanceDz() - cosAng * disAng2atoms23.dDistanceDz());
            DsigmCosAng2Dx4 = sigma.s_tag() * disAng2atoms43.invDistance() * (disAng2atoms23.dDistanceDx() - cosAng * disAng2atoms43.dDistanceDx());
            DsigmCosAng2Dy4 = sigma.s_tag() * disAng2atoms43.invDistance() * (disAng2atoms23.dDistanceDy() - cosAng * disAng2atoms43.dDistanceDy());
            DsigmCosAng2Dz4 = sigma.s_tag() * disAng2atoms43.invDistance() * (disAng2atoms23.dDistanceDz() - cosAng * disAng2atoms43.dDistanceDz());
            DsigmCosAng2Dx3 = -DsigmCosAng2Dx2 - DsigmCosAng2Dx4;
            DsigmCosAng2Dy3 = -DsigmCosAng2Dy2 - DsigmCosAng2Dy4;
            DsigmCosAng2Dz3 = -DsigmCosAng2Dz2 - DsigmCosAng2Dz4;
        } else {
            derivativeAng2Zero = true;
            DsigmCosAng2Dx4 = DsigmCosAng2Dy4 = DsigmCosAng2Dz4 = DsigmCosAng2Dx3 = DsigmCosAng2Dy3 =
                    DsigmCosAng2Dz3 = DsigmCosAng2Dx2 = DsigmCosAng2Dy2 = DsigmCosAng2Dz2 = 0.0;
        }
    }

    private final void atomListTo1_4numbering() {
        if (!atomList.atomAt(2).type().isHydrogen() && !atomList.atomAt(3).type().isHydrogen()) {  // case 1:   base---O...O---base    or    base---O...N---base
            a1 = atomList.atomAt(2);
            a2 = atomList.atomAt(0);
            a3 = atomList.atomAt(1);
            a4 = atomList.atomAt(3);
            a1toAtomList = 2;
            a2toAtomList = 0;
            a3toAtomList = 1;
            a4toAtomList = 3;
        } else if (!atomList.atomAt(2).type().isHydrogen() && atomList.atomAt(3).type().isHydrogen()) { // case 2:   base---O...H---N
            a1 = atomList.atomAt(2);
            a2 = atomList.atomAt(0);
            a3 = atomList.atomAt(3);
            a4 = atomList.atomAt(1);
            a1toAtomList = 2;
            a2toAtomList = 0;
            a3toAtomList = 3;
            a4toAtomList = 1;
        } else if (atomList.atomAt(2).type().isHydrogen() && !atomList.atomAt(3).type().isHydrogen()) { // case 3:   N---H...O---base
            a1 = atomList.atomAt(0);
            a2 = atomList.atomAt(2);
            a3 = atomList.atomAt(1);
            a4 = atomList.atomAt(3);
            a1toAtomList = 0;
            a2toAtomList = 2;
            a3toAtomList = 1;
            a4toAtomList = 3;
        } else if (atomList.atomAt(2).type().isHydrogen() && atomList.atomAt(3).type().isHydrogen()) { // case 4:   N---H...H---N
            a1 = atomList.atomAt(0);
            a2 = atomList.atomAt(2);
            a3 = atomList.atomAt(3);
            a4 = atomList.atomAt(1);
            a1toAtomList = 0;
            a2toAtomList = 2;
            a3toAtomList = 3;
            a4toAtomList = 1;
        } else
            throw new RuntimeException("If this line is reached, something is definitely wrong");
    }

    public boolean zeroDerivative() {
        return zeroDerivative; }
}


  
