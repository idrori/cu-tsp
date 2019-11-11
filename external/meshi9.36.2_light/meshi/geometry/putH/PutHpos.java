/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry.putH;

import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;

import java.util.*;

/**
 * This procedure is was copied and manually translated to java, from Ron Elber's puth program available in the MOIL package.
 * The original fortarn (both code and comments) is kept as comments.
 */
/*subroutine pos(success,ih)
  c
  c Placing hydrogens
  c We made use of known coordinates and the data stored in the connectivity block
  c Only bond lengths and bond angles are used (i.e. NO torsions are employed)
  c In case of missing data, the direction are picked at random
  c (e.g. water molecule)
  c No information isOn possible hydrogen bonds candidate is taken into
  c account in this admittedly oversimplifed placement
  c
  c 
                                   include 'COMMON/LENGTH.BLOCK'
                   include 'COMMON/CONNECT.BLOCK'
                   include 'COMMON/UNITS.BLOCK'
                   include 'COMMON/COORD.BLOCK'
                   include 'COMMON/LINE.BLOCK'
                   include 'COMMON/DEBUG.BLOCK'*/

public class PutHpos {
    private static Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

    public static PutHposLog pos(Atom hydrogenAtom,
                                 BondParametersList bondParametersList,
                                 AngleParametersList angleParametersList)  {

        /* integer ih
                     logical success
                     character*3 name
                     integer namel
                     integer level
                     integer i,j,nbondj,iseed
                     integer ibj(maxptbn)
                     real ran2
                     double precision reqhj
                     double precision cos1,cos2,cos3
                     double precision sin1,sin2,sin3
                     double precision gamma,delta
                     double precision e1(3),e2(3),e3(3),e4(3)
                     double precision pi
                     data iseed/-112233/
                     save iseed
                     double pi = 4.d0*Math.atan(1.d0);
                     debug  = .false.
                     name   = 'pos'
                     namel  = 3
                     success = .false. */
        /* c
         c Find out a bonded particle to the neig particle (a neighbor to the pt)
         c It is assumed that hydrogen has only one bond (i.e. no Boron)
         c */
        /*do 1 i=1,nb
                       if (ib1(i).eq.ih) then
                       if (coor(1,ib2(i)).gt.9998.d0) then
                       level = 1
                       call alert(name,namel,'H bonded to undefined crd',25,level)
                       return
                       end if
                       j     = ib2(i)
                       reqhj = req(i)
                       go to 2
                       else if (ib2(i).eq.ih) then
                       if (coor(1,ib1(i)).gt.9998.d0) then
                       level = 1
                       call alert(name,namel,'H bonded to undefined crd',25,level)
                       return
                       end if
                       j     = ib1(i)
                       reqhj = req(i)
                       go to 2
                       end if
                       1	continue */
        Atom iHydrogen = hydrogenAtom;
        if (!iHydrogen.type().isHydrogen())
            throw new RuntimeException("Weird argument to pos " + iHydrogen + " " + iHydrogen.type());
        if (iHydrogen.bonded().size() != 1)
            throw new RuntimeException("Weird number of getAtoms bonded to " + iHydrogen + " " + iHydrogen.bonded().size());
        Atom jNeighbor = iHydrogen.bonded().get(0);
        if (!iHydrogen.nowhere())
            throw new RuntimeException("Hydrogen atom already has coordinates.\n" + iHydrogen + "\n" + iHydrogen.core);
        if (jNeighbor.type().isHydrogen())
            throw new RuntimeException("Weird hydrogen " + jNeighbor + " bound to " + iHydrogen);
        if (jNeighbor.nowhere())
            return null;
        BondParameters bondParameters = (BondParameters) bondParametersList.getParameters(new BondParameters(jNeighbor.type(),
                iHydrogen.type()));
        if (bondParameters == null) throw new RuntimeException("Cannot find parameters for the bond of " +
                iHydrogen + " and " + jNeighbor);
        double targetDistance = bondParameters.target;
        /*level = 1
                        call alert(name,namel,'H bonded to nothing !?',22,level)
                        return
                        2	continue */
        /* c
         c Find out all particles bonded to j (excluding the hydrogen)
         c count only bonds to a particle with defined coordinates
         c */
        /*nbondj = 0
                      do 3 i=1,nb
                      if ((ib1(i).eq.j .and. ib2(i).ne.ih) .or.
                      1		  (ib2(i).eq.j .and. ib1(i).ne.ih)) then
                      if (ib1(i).eq.j) then
                      if (coor(1,ib2(i)).gt.9998.d0) go to 3
                      nbondj = nbondj + 1
                      ibj(nbondj)   = ib2(i)
                      else
                      if (coor(1,ib1(i)).gt.9998.d0) go to 3
                      nbondj = nbondj + 1
                      ibj(nbondj)   = ib1(i)
                      end if
                      end if
                      3	continue*/
        /* c
         c Now check how many particles are bonded to the j-th particle
         c */

        /*if (nbondj.eq.0) then
                        c generating position for the hydrogen connected to an isolated center
                        c use correct distance but random direction
                        c
                        e1(1) = ran2(iseed)
                        e1(2) = ran2(iseed)
                        e1(3) = ran2(iseed)
                        cos1  = e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3)
                        cos1  = 1.d0/dsqrt(cos1)
                        e1(1) = e1(1)*cos1
                        e1(2) = e1(2)*cos1
                        e1(3) = e1(3)*cos1
                        coor(1,ih) = coor(1,j) + reqhj*e1(1)
                        coor(2,ih) = coor(2,j) + reqhj*e1(2)
                        coor(3,ih) = coor(3,j) + reqhj*e1(3)
                      */
        if (jNeighbor.bonded().size() == 2) {
//            throw new NotEnoughBoundAtomsException("An error while trying to assign coordinates to " +
//                    iHydrogen + "\n" +
//                    "The neighbor heavy atom " + jNeighbor + "\n" +
//                    "has only " + jNeighbor.bonded().size() +
//                    " boundAtoms\n" +
//                    "Hydrogen neighbor with less than 3 bonds is not supported yet\n" +
//                    "This needs to be changed in order to support all-atom-model");
            if (!jNeighbor.nowhere()) {
                iHydrogen.setXYZ(jNeighbor.x()+MeshiProgram.randomNumberGenerator().nextDouble(),
                        jNeighbor.y()+MeshiProgram.randomNumberGenerator().nextDouble(),
                        jNeighbor.z()+MeshiProgram.randomNumberGenerator().nextDouble());
            }

            return null;
        }
        if ((jNeighbor.bonded().size() >= 4) || (jNeighbor.bonded().size() <= 1))
            throw new RuntimeException("An error while trying to assign coordinates to " +
                    iHydrogen + "\n" +
                    "The neighbor heavy atom " + jNeighbor + "\n" +
                    "has " + jNeighbor.bonded().size() +
                    " boundAtoms\n" +
                    "Hydrogen neighbor with more than 3 bonds is not supported yet\n" +
                    "This needs to be changed in order to support all-atom-model");
        if (jNeighbor.bonded().size() == 3) {
            int nNowhere = 0;
            for (Iterator bondedIter = jNeighbor.bonded().iterator(); bondedIter.hasNext();) {
                Atom atom = (Atom) bondedIter.next();
                if (atom.nowhere()) nNowhere++;
            }
            if (nNowhere == 0) throw new RuntimeException("Non of the neighbors for the heavy neighbor of \n" +
                    "the hydrogen has coordinates.");
            if (nNowhere == 3) throw new RuntimeException("This is really weird 4  "+hydrogenAtom); //Cannot assign coordinates
            if (nNowhere == 2) { // The hydrogen position is not uniquely defined.
                int counter = 0;
                while (true) { /* It is expected that we will get out of this loop by the return statement
				  most probably during the first pass. Otherwise, an exception will be thrown when the 
				  loop counter reaches 100*/
                    if (counter >= 100) throw new RuntimeException("failed to generate non-lineary-dependent-vectors");
                    counter++;

                    Atom kAtom = null;
                    for (Iterator bondedIter = jNeighbor.bonded().iterator(); bondedIter.hasNext();) {
                        Atom atom = (Atom) bondedIter.next();
                        if (!atom.nowhere()) kAtom = atom;
                    }
                    if (kAtom == null) throw new RuntimeException("weird 8");
                    AngleParameters key = new AngleParameters(iHydrogen.type(),
                            jNeighbor.type(),
                            kAtom.type());
                    AngleParameters angleParameters = (AngleParameters) angleParametersList.getParameters(key);
                    if (angleParameters == null) throw new RuntimeException("Cannot find andle parameters for " +
                            iHydrogen + " , " + jNeighbor + " and " +
                            kAtom);
                    double targetAngle = angleParameters.target;

                    if (Math.abs(targetAngle - Math.PI) < 0.000001) throw new RuntimeException("Flat angle");
                    /*else if (nbondj.eq.1) then
                         c search for the angle correspoding to ih-neigh-ibj(1)
                         do 4 i=1,nangl
                         if (iangl2(i).eq.j) then
                         if ((ih.eq.iangl1(i) .and. ibj(1).eq.iangl3(i)) .or.
                         1	      (ih.eq.iangl3(i) .and. ibj(1).eq.iangl1(i))) then
                         c get the cosine of the angle between the OTHER bond to j and the bond between
                         c j and the hydrogen
                         if (dabs(angleq(i)-pi).lt.1.d-6
                         1	      .or. dabs(angleq(i)).lt.1.d-6) then
                         level = 1
                         call alert(name,namel,'LINEAR ANGLE',12,level)
                         end if*/
                    double cos1 = Math.cos(targetAngle);
                    /*cos1    = dabs(cos(angleq(i)))
                          else
                          goto 4
                          end if
                        */
                    double sin1 = Math.sqrt(1 - cos1 * cos1);
                    /*sin1    = dsqrt(1.d0-cos1*cos1) */
                    //c generate a unit vector along the OTHER bond
                    double[] e2 = new double[3];
                    try {
                        e2[0] = jNeighbor.x() - kAtom.x();
                        e2[1] = jNeighbor.y() - kAtom.y();
                        e2[2] = jNeighbor.z() - kAtom.z();
                    } catch (RuntimeException ex) {
                        System.out.println("A problem in e2 vector creation\n" +
                                "jNeighbor = " + jNeighbor + " nowhere = " + jNeighbor.nowhere() + "\n" +
                                "kAtom     = " + kAtom + " nowhere = " + kAtom.nowhere() + "\n" +
                                "iHydrogen = " + iHydrogen + " nowhere = " + iHydrogen.nowhere() + "\n");
                        ex.printStackTrace();
                        throw ex;
                    }
                    double cos2 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
                    cos2 = 1 / Math.sqrt(cos2);
                    e2[0] = e2[0] * cos2;
                    e2[1] = e2[1] * cos2;
                    e2[2] = e2[2] * cos2;
                    /* e2(1)   = coor(1,j) - coor(1,ibj(1))
                           e2(2)   = coor(2,j) - coor(2,ibj(1))
                           e2(3)   = coor(3,j) - coor(3,ibj(1))
                           cos2    = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
                           cos2    = 1.d0/dsqrt(cos2)
                           e2(1)   = e2(1)*cos2
                           e2(2)   = e2(2)*cos2
                           e2(3)   = e2(3)*cos2*/
                    /* c generate a vector in a random direction, to initiate sampling in the
                c allowed circle */
                    double[] e3 = new double[3];
                    e3[0] = randomNumberGenerator.nextDouble();
                    e3[1] = randomNumberGenerator.nextDouble();
                    e3[2] = randomNumberGenerator.nextDouble();
                    double cos3 = e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2];
                    cos3 = 1 / Math.sqrt(cos3);
                    e3[0] = e3[0] * cos3;
                    e3[1] = e3[1] * cos3;
                    e3[2] = e3[2] * cos3;
                    cos2 = e2[0] * e3[0] + e2[1] * e3[1] + e2[2] * e3[2];
                    if (Math.abs(cos2 - 1) < 0.000001) {
                        System.out.println("PutHpos.pos warrning: Vectors linearly dependent " + counter);
                        /*e3(1)   = ran2(iseed)
                           e3(2)   = ran2(iseed)
                           e3(3)   = ran2(iseed)
                           cos3    = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
                           cos3    = 1.d0/dsqrt(cos3)
                           e3(1)   = e3(1)*cos3
                           e3(2)   = e3(2)*cos3
                           e3(3)   = e3(3)*cos3
                           cos2    = e2(1)*e3(1) +  e2(2)*e3(2) + e2(3)*e3(3)
                           if ((dabs(cos2-1.d0).lt.1.d-6) .or.
                           1	     (dabs(cos2+1.d0).lt.1.d-6)) then
                           level = 0
                           call alert(name,namel,'Vectors linearly dependent',27,level)
                           end if */
                    } else {
                        /* c orthonormalize e3 with respect to the existing bond direction */
                        e3[0] = e3[0] - cos2 * e2[0];
                        e3[1] = e3[1] - cos2 * e2[1];
                        e3[2] = e3[2] - cos2 * e2[2];
                        cos3 = e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2];
                        cos3 = 1 / Math.sqrt(cos3);
                        e3[0] = e3[0] * cos3;
                        e3[1] = e3[1] * cos3;
                        e3[2] = e3[2] * cos3;
                        /* e3(1)   = e3(1) - cos2*e2(1)
                          e3(2)   = e3(2) - cos2*e2(2)
                          e3(3)   = e3(3) - cos2*e2(3)
                          cos3    = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
                          cos3    = 1.d0/dsqrt(cos3)
                          e3(1)   = e3(1)*cos3
                          e3(2)   = e3(2)*cos3
                          e3(3)   = e3(3)*cos3*/
                        /*c generate the hydrogen position. the displacement vector
                c as compared to the j atom coordinateis a linear combination
                c of the two orthonormalized vectors e2 and e3*/
                        /* coor(1,ih)   = coor(1,j) + reqhj*(sin1*e3(1)+cos1*e2(1))
                          coor(2,ih)   = coor(2,j) + reqhj*(sin1*e3(2)+cos1*e2(2))
                          coor(3,ih)   = coor(3,j) + reqhj*(sin1*e3(3)+cos1*e2(3))*/
                        iHydrogen.setXYZ(jNeighbor.x() + targetDistance * (sin1 * e3[0] + cos1 * e2[0]),
                                jNeighbor.y() + targetDistance * (sin1 * e3[1] + cos1 * e2[1]),
                                jNeighbor.z() + targetDistance * (sin1 * e3[2] + cos1 * e2[2]));
                        iHydrogen.setStatus(AtomStatus.NORMAL);
//			iHydrogen.setTemperatureFactor(jNeighbor.temperatureFactor());
                        if (Math.abs(targetDistance - iHydrogen.distanceFrom(jNeighbor)) > 0.05)
                            throw new RuntimeException("Failed to assign coordinates to " + iHydrogen + " its distance from " +
                                    jNeighbor + " is " + iHydrogen.distanceFrom(jNeighbor) + "\n" +
                                    "Should have been " + targetDistance);
                        return new PutHposLog(iHydrogen, jNeighbor, 3, 2);
                        /* go to 5
                          end if
                          4	 continue
                          level = 1
                          call alert(name,namel,'ANGLE NOT FOUND',15,level)
                          5	 continue*/
                    }
                }
            } else { // the position of the hydrogen is uniquely (up to mirror symmetry) defined.
                /*else if (nbondj.eq.2) then */
                /*c only two bonds are used in the present verison. A danger that does not exist
              c in the present force field (extended atom model) but does exist
              c for all atom case is
              c the possibility of inverted improper torsion.*/
                /*cos1 = 999.d0
                         cos2 = 999.d0
                         cos3 = 999.d0

                                         do 6 i=1,nangl
                          if (iangl2(i).eq.j) then
                          if (dabs(angleq(i)-pi).lt.1.d-6 .or.
                          1	    dabs(angleq(i)).lt.1.d-6) then
                          level = 1
                          call alert(name,namel,'LINEAR ANGLE',12,level)
                          end if
                          if ((iangl1(i).eq.ih .and. iangl3(i).eq.ibj(1)) .or.
                          1		(iangl1(i).eq.ibj(1) .and. iangl3(i).eq.ih)) then
                          cos1 = cos(angleq(i))
                          else if ((iangl1(i).eq.ih .and. iangl3(i).eq.ibj(2)) .or.
                          1		(iangl1(i).eq.ibj(2) .and. iangl3(i).eq.ih)) then
                          cos2 = cos(angleq(i))
                          else if ((iangl1(i).eq.ibj(1) .and. iangl3(i).eq.ibj(2)) .or.
                          1		(iangl1(i).eq.ibj(2) .and. iangl3(i).eq.ibj(1))) then
                          cos3 = cos(angleq(i))
                          end if
                          end if
                          6	 continue*/
                Atom kAtom, lAtom;
                if (jNeighbor.bonded().get(0) == iHydrogen) {
                    kAtom = jNeighbor.bonded().get(1);
                    lAtom = jNeighbor.bonded().get(2);
                } else if (jNeighbor.bonded().get(1) == iHydrogen) {
                    kAtom = jNeighbor.bonded().get(0);
                    lAtom = jNeighbor.bonded().get(2);
                } else if (jNeighbor.bonded().get(2) == iHydrogen) {
                    kAtom = jNeighbor.bonded().get(0);
                    lAtom = jNeighbor.bonded().get(1);
                } else throw new RuntimeException("This is really weird 3");
                AngleParameters parametersIJK = (AngleParameters) angleParametersList.getParameters(new AngleParameters(iHydrogen.type(),
                        jNeighbor.type(),
                        kAtom.type()));
                if (parametersIJK == null) throw new RuntimeException("Cannot find angle parameters for " +
                        iHydrogen + " , " + jNeighbor + " and " +
                        kAtom);
                AngleParameters parametersIJL = (AngleParameters) angleParametersList.getParameters(new AngleParameters(iHydrogen.type(),
                        jNeighbor.type(),
                        lAtom.type()));
                if (parametersIJL == null) throw new RuntimeException("Cannot find angle parameters for " +
                        iHydrogen + " , " + jNeighbor + " and " +
                        lAtom);
                AngleParameters parametersKJL = (AngleParameters) angleParametersList.getParameters(new AngleParameters(kAtom.type(),
                        jNeighbor.type(),
                        lAtom.type()));
                if (parametersKJL == null) throw new RuntimeException("Cannot find angle parameters for " +
                        kAtom + " , " + jNeighbor + " and " +
                        lAtom);
                double targetIJK = parametersIJK.target;
                double targetIJL = parametersIJL.target;
                double targetKJL = parametersKJL.target;
                double cos1ijk = Math.cos(targetIJK);
                double cos2ijl = Math.cos(targetIJL);
                double cos3kjl = Math.cos(targetKJL);
                /*if (cos1.gt.1.d0 .or. cos2.gt.1.d0 .or. cos3.gt.1d0) then
                                  level = 1
                                  call alert(name,namel,'Undefined cosine values',23,level)
                                  end if*/
                double[] e1jk = new double[3], e2jl = new double[3], e3ji = new double[3];
                e1jk[0] = kAtom.x() - jNeighbor.x();
                e1jk[1] = kAtom.y() - jNeighbor.y();
                e1jk[2] = kAtom.z() - jNeighbor.z();
                /* e1(1) = coor(1,ibj(1)) - coor(1,j)
                                  e1(2) = coor(2,ibj(1)) - coor(2,j)
                                  e1(3) = coor(3,ibj(1)) - coor(3,j)*/
                /* c sin1 is used here as a buffer for normalization */
                double sin1 = e1jk[0] * e1jk[0] + e1jk[1] * e1jk[1] + e1jk[2] * e1jk[2];
                sin1 = 1 / Math.sqrt(sin1);
                e1jk[0] *= sin1;
                e1jk[1] *= sin1;
                e1jk[2] *= sin1;
                /* sin1  = e1(1)*e1(1) + e1(2)*e1(2) + e1(3)*e1(3)
                                sin1  = 1.d0/dsqrt(sin1)
                                e1(1) = e1(1)*sin1
                                e1(2) = e1(2)*sin1
                                e1(3) = e1(3)*sin1 */
                e2jl[0] = lAtom.x() - jNeighbor.x();
                e2jl[1] = lAtom.y() - jNeighbor.y();
                e2jl[2] = lAtom.z() - jNeighbor.z();
                /* e2(1) = coor(1,ibj(2)) - coor(1,j)
                                e2(2) = coor(2,ibj(2)) - coor(2,j)
                                e2(3) = coor(3,ibj(2)) - coor(3,j)*/
                /* c sin2 is used here as a buffer for normalization */
                double sin2 = e2jl[0] * e2jl[0] + e2jl[1] * e2jl[1] + e2jl[2] * e2jl[2];
                sin2 = 1 / Math.sqrt(sin2);
                e2jl[0] *= sin2;
                e2jl[1] *= sin2;
                e2jl[2] *= sin2;
                /* sin2  = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
                                   sin2  = 1.d0/dsqrt(sin2)
                                   e2(1) = e2(1)*sin2
                                   e2(2) = e2(2)*sin2
                                   e2(3) = e2(3)*sin2*/
                /* c orthonormalize e2 with respect to e1 */
                sin1 = e1jk[0] * e2jl[0] + e1jk[1] * e2jl[1] + e1jk[2] * e2jl[2];
                e2jl[0] -= sin1 * e1jk[0];
                e2jl[1] -= sin1 * e1jk[1];
                e2jl[2] -= sin1 * e1jk[2];
                /* sin1  = e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3)
                                    e2(1) = e2(1) - sin1*e1(1)
                                    e2(2) = e2(2) - sin1*e1(2)
                                    e2(3) = e2(3) - sin1*e1(3)*/
                //c normalize the new e2
                sin2 = e2jl[0] * e2jl[0] + e2jl[1] * e2jl[1] + e2jl[2] * e2jl[2];
                sin2 = 1 / Math.sqrt(sin2);
                e2jl[0] *= sin2;
                e2jl[1] *= sin2;
                e2jl[2] *= sin2;
                /* sin2  = e2(1)*e2(1) + e2(2)*e2(2) + e2(3)*e2(3)
                                      sin2  = 1.d0/dsqrt(sin2)
                                      e2(1) = e2(1)*sin2
                                      e2(2) = e2(2)*sin2
                                      e2(3) = e2(3)*sin2*/
                //c generate the third unit vector by a vector products of e1 and e2
                //c
                e3ji[0] = e1jk[1] * e2jl[2] - e1jk[2] * e2jl[1];
                e3ji[1] = e1jk[2] * e2jl[0] - e1jk[0] * e2jl[2];
                e3ji[2] = e1jk[0] * e2jl[1] - e1jk[1] * e2jl[0];
                /* e3(1) = e1(2)*e2(3) - e1(3)*e2(2)
                                       e3(2) = e1(3)*e2(1) - e1(1)*e2(3)
                                       e3(3) = e1(1)*e2(2) - e1(2)*e2(1)*/
                // c check that e3 is normalized
                double sin3 = e3ji[0] * e3ji[0] + e3ji[1] * e3ji[1] + e3ji[2] * e3ji[2];
                if (Math.abs(sin3) < 0.000001) throw new RuntimeException("Fishy vector product");
                /* sin3 = e3(1)*e3(1) + e3(2)*e3(2) + e3(3)*e3(3)
                                       if (dabs(sin3).lt.1.d-6) then
                                       level = 1
                                       call alert(name,namel,'Fishy vector product',20,level)
                                       end if */
                /* c A unit vector along the bond between j and ih is expanded in terms
               c of e1,e2 and e3 as follows - R(ih)-R(j) = cos1 e1 + gamma e2 + delta e3
               c where gamma=(cos2 -cos1*cos3)/sqrt(1-cos3^2) and delta=sqrt(1-gamma^2-cos1^2) */
                double gamma = (cos2ijl - cos1ijk * cos3kjl) / Math.sqrt(1 - cos3kjl * cos3kjl);
                double delta = gamma * gamma + cos1ijk * cos1ijk;
                if (1 - delta < 0.0000000001) delta = 0;
                else delta = Math.sqrt(1 - delta);
                /* gamma = (cos2-cos3*cos1)/dsqrt(1.d0 - cos3*cos3)
                                       delta = gamma*gamma+cos1*cos1
                                       if (1.d0-delta.lt.1.d-10) then
                                       delta = 0.d0
                                       else
                                       delta = dsqrt(1.d0-delta)
                                       end if */
                iHydrogen.setXYZ(jNeighbor.x() + targetDistance * (cos1ijk * e1jk[0] + gamma * e2jl[0] + delta * e3ji[0]),
                        jNeighbor.y() + targetDistance * (cos1ijk * e1jk[1] + gamma * e2jl[1] + delta * e3ji[1]),
                        jNeighbor.z() + targetDistance * (cos1ijk * e1jk[2] + gamma * e2jl[2] + delta * e3ji[2]));
                iHydrogen.setStatus(AtomStatus.NORMAL);
//		iHydrogen.setTemperatureFactor(jNeighbor.temperatureFactor());
                /* coor(1,ih)=coor(1,j)+
                                       reqhj*(cos1*e1(1)+gamma*e2(1)+delta*e3(1))
                                       coor(2,ih)=coor(2,j)+
                                       reqhj*(cos1*e1(2)+gamma*e2(2)+delta*e3(2))
                                       coor(3,ih)=coor(3,j)+
                                       reqhj*(cos1*e1(3)+gamma*e2(3)+delta*e3(3))*/
            }
        }
        if (Math.abs(targetDistance - iHydrogen.distanceFrom(jNeighbor)) > 0.05)
            throw new RuntimeException("Failed to assign coordinates to " + iHydrogen + " its distance from " +
                    jNeighbor + " is " + iHydrogen.distanceFrom(jNeighbor) + "\n" +
                    "Should have been " + targetDistance);
        return new PutHposLog(iHydrogen, jNeighbor, 3, 3);
    }
}
/*
if (jNeighbor.bonded().size == 4) {
	                                                         // else if (nbondj.eq.3) then
	    //c
	    //c calculate eigenvectors along the exisiting bonds pointing to the central
	    //c bond.
	    //c
	    e1
	    e1(1) = coor(1,j) - coor(1,ibj(1))
	 e1(2) = coor(2,j) - coor(2,ibj(1))
	 e1(3) = coor(3,j) - coor(3,ibj(1))

	 e2(1) = coor(1,j) - coor(1,ibj(2))
	 e2(2) = coor(2,j) - coor(2,ibj(2))
	 e2(3) = coor(3,j) - coor(3,ibj(2))

	 e3(1) = coor(1,j) - coor(1,ibj(3))
	 e3(2) = coor(2,j) - coor(2,ibj(3))
	 e3(3) = coor(3,j) - coor(3,ibj(3))

c calculate a unit vector in the direction of the missing bond. It is
c assumed that it is in the direction of the average of the e1 to e3
c eigenvectors which (hopefully) is not too bad
c
	 e4(1) = e1(1) + e2(1) + e3(1)
	 e4(2) = e1(2) + e2(2) + e3(2)
	 e4(3) = e1(3) + e2(3) + e3(3)
	 delta = 1.d0/dsqrt(e4(1)*e4(1)+e4(2)*e4(2)+e4(3)*e4(3))
	 e4(1) = e4(1)*delta
	 e4(2) = e4(2)*delta
	 e4(3) = e4(3)*delta
	 coor(1,ih) = coor(1,j) + reqhj*e4(1)
	 coor(2,ih) = coor(2,j) + reqhj*e4(2)
	 coor(3,ih) = coor(3,j) + reqhj*e4(3)

	else
	 write(stdo,*)' Nbondj = ',nbondj
	 call alert('pos',3,' Impossible # of bonds to H',26,1)
	end if
	success = .true.
	return
	end
*/

