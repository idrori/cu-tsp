package meshi.molecularElements.loops;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Random;

import meshi.energy.EnergyCreator;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.angle.AngleParameters;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParameters;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.plane.PlaneParameters;
import meshi.energy.simpleEnergyTerms.plane.PlaneParameters.Isomer;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.Coordinates;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import meshi.util.KeyWords;
import meshi.parameters.AtomType;
import meshi.parameters.MeshiPotential;

public class AtomFinding {
	public final static double PI = Math.PI;
    private static boolean debug = false;

	/*
	 * public enum AtomFindingRoutine{//Tommer 22.10.14 XO_LIKE, // this is XH ,
	 * YOH, NOD, QOE, DO. Tommer 30.11.14 //FD_LIKE, //this is FD, YD. Tommer
	 * 30.11.14 SECOND_FD_LIKE, FIRST_FD_LIKE //this is a prototype for adding
	 * based on two getAtoms and an angle }
	 */
	

	/*private enum ORDER {
		A, B, G, D, E, Z, H; // should this be private? Tommer 29.11.14
		public static ORDER stringToEnum(String s) {
			switch (s.charAt(0)) {
			case 'A':
				return A;
			case 'B':
				return B;
			case 'G':
				return G;
			case 'D':
				return D;
			case 'E':
				return E;
			case 'Z':
				return Z;
			case 'H':
				return H;
			default:
				throw new RuntimeException(
						"error finding lexicographical distance from backbone");

			}
		}
	}

	private class AtomLexicographicComparator implements Comparator<Atom> {
		public int compare(Atom atom1, Atom atom2) {
			if (!((atom1.type().toString().endsWith("CA")
					|| atom1.type().toString().endsWith("CB")
					|| atom1.type().toString().endsWith("CG")
					|| atom1.type().toString().endsWith("CD")
					|| atom1.type().toString().endsWith("CE")
					|| atom1.type().toString().endsWith("CZ")
					|| atom1.type().toString().endsWith("CH"))
					&&( atom2.type().toString().endsWith("CA")
					|| atom2.type().toString().endsWith("CB")
					|| atom2.type().toString().endsWith("CG")
					|| atom2.type().toString().endsWith("CD")
					|| atom2.type().toString().endsWith("CE")
					|| atom2.type().toString().endsWith("CZ") 
					|| atom2.type().toString().endsWith("CH")))) {
				
				if (atom1.residueName().equals("PRO")){
				//	System.out.println("Proline temporarily treated as start of loop\n"+
						// "atom1.type(): "+atom1.type()+"res: "+atom1.residueName());
					return -1;
				}else{
					if (atom2.residueName().equals("PRO")){
						//System.out.println("Proline temporarily treated as start of loop\n"+
							//	"atom2.type(): "+atom2.type()+"res: "+atom2.residueName());
						return 1;
					}else{
						throw new RuntimeException(
								"input atom is not recognized as an ordered carbon"+
										"\natom1: "+atom1+atom1.type()+"\natom2: "+ atom2+atom2.type());
					}
				}		
			}
			
			ORDER letFirstAtom, letSecondAtom;
			
			//System.out.println("atomType: " + atom1.type().toString()
				//	+ "\n substring: " + atom1.type().toString().substring(1));
			// Check. Tommer 3.1.15
		//	System.out.println("atomType2: " + atom2.type().toString()
			//		+ "\n substring: " + atom2.type().toString().substring(1));
					
			letFirstAtom = ORDER.stringToEnum(atom1.type().toString()
					.substring(2));
			letSecondAtom = ORDER.stringToEnum(atom2.type().toString()
					.substring(2));

			if (letFirstAtom.ordinal() > letSecondAtom.ordinal()) {
				return 1;
			} else {
				if (letFirstAtom.ordinal() < letSecondAtom.ordinal()) {
					return -1;
				}
				return 0;
			}

		}
	}*/

	public static double getBondDistance(Atom atom1, Atom atom2, BondParametersList parametersList) {
		AtomPair atomPair = new AtomPair(atom1, atom2);
		BondParameters pairPar = (BondParameters) parametersList
				.parameters(atomPair);
		/*
		 * System.out.println("null pairPar? "+(pairPar==null)+
		 * ", null atomPair "+(atomPair==null));//check Tommer 27.10.14
		 */// Not Needed! 25.12.14
		return pairPar.target;

	}

	private static double getAngle(Atom atom1, Atom atom2, Atom atom3,
                                   AngleParametersList angleParametersList) {
		AngleParameters key = new AngleParameters(atom1.type(), atom2.type(),
				atom3.type());// what about the order? Tommer 29.11.14
		AngleParameters angleParameters = (AngleParameters) angleParametersList
				.getParameters(key);
		if (angleParameters == null)
			throw new RuntimeException("Cannot find angle parameters for "
					+ atom1 + " , " + atom2 + " and " + atom3);
		double targetAngle = angleParameters.target;
		// throw new RuntimeException("angle is "+targetAngle); Tommer check
		// 30.11.14
		return targetAngle;
	}
	
	
	private static Coordinates atomLocation2Bonded(double alpha,
			double gamma, Coordinates farCoord, Coordinates nearCoord,
			Coordinates secondNeighborCoord, double a, double b, double topAngle){
		double c = nearCoord.distanceFrom(secondNeighborCoord);
		double c0 = thirdTriangleEdge(a, b, topAngle);
		double beta = angleInTriangle(a, b, c0);
		double sinTheta= Coordinates.sineOfOut_of_PlaneAngle(PI-alpha, PI-beta, gamma);
		double gamma2= Math.acos(Coordinates.cosineOfgamma2(sinTheta, PI - alpha));
		double cosTheta= Math.sqrt(1-sinTheta*sinTheta);
		Coordinates xNorm = farCoord.differenceVector(nearCoord).normalize();
		Coordinates zNorm = xNorm.vectorMult(nearCoord.differenceVector(secondNeighborCoord)).normalize();
		Coordinates yNorm = zNorm.vectorMult(xNorm);
		Coordinates result = nearCoord.vectorAdd(zNorm.scalarMult(b*sinTheta));
		result = result.vectorAdd(xNorm.scalarMult(-b*cosTheta*Math.cos(gamma2)));
		result = result.vectorAdd(yNorm.scalarMult(b*Math.sin(gamma2)));
		if(result.x()!=result.x())
			throw new RuntimeException("alpha: "+alpha+"\nbeta: "+beta+"\ngamma: "+gamma+
					"cosineOfgamma2(sinTheta, alpha): "+Coordinates.cosineOfgamma2(sinTheta, alpha)+
					"\ngamma2: "+gamma2+"\nsinTheta: "+sinTheta+"\ncosTheta: "+cosTheta+
					"\na: "+a+" b: "+b+" c, c0: "+ c+"  "+c0+"\nxNorm: "+xNorm+
					"\nyNorm: "+yNorm+"\nresult: "+result);
		
		/*if (1==1)
			throw new RuntimeException("c: "+c+" beta: "+beta+" sinTheta: "+sinTheta+" gamma2: "+gamma2+
					" cosTheta: "+cosTheta+ "\nresult: "+result);*/
		return result;
		
	}
	
	private static double angleInTriangle(double a, double b, double c){
		double cosAngle= (b*b+c*c-a*a)/(2*b*c);
		
		/*if(!((a<4.9357)&&(a>4.9355)))
			throw new RuntimeException("a: "+a+" b: "+b+" c: "+c+" cosAngle: "+
					cosAngle+"\nacos(cos): "+Math.acos(cosAngle));*/
		return Math.acos(cosAngle);
		//returns the angle facing a. Tommer 24.1.15
	}
	
	private static double thirdTriangleEdge(double a, double b, double angle){
		return Math.sqrt(a*a+b*b-2*a*b*Math.cos(angle));
	}
	
	private static PlaneParameters.Isomer planeAngleType(AtomType type1, AtomType type2,
			AtomType type3, AtomType type4, PlaneParametersList parametersList){
		//this function needs to be given the getAtoms in the right order. Tommer 24.1.15

		PlaneParameters parameters=(PlaneParameters) parametersList.parameters(type1, type2, type3, type4);
		if(parameters==null) return null;//does this do the right thing?? Tommer 9.1.15
		return parameters.isomer;
	}
	
	
	private static PlaneParameters.Isomer planeAngleTypeBackAndForth(AtomType type1, AtomType type2,
			AtomType type3, AtomType type4, PlaneParametersList parametersList){
		PlaneParameters.Isomer ans;
		if ((ans=planeAngleType(type1, type2, type3, type4, parametersList))!=null){
			return ans;
		}
		else{
			ans=planeAngleType(type4, type3, type2, type1, parametersList);
			return ans;
		}
	}
	
	/*private static PlaneParameters.Isomer planeAngleType24(AtomType type1, AtomType type2,
			AtomType type3, AtomType type4, CommandList commands, ArrayList<AtomType> orderedTypes){
		ArrayList<AtomType> typeList=new ArrayList<AtomType>();
		typeList.add(type1);
		typeList.add(type2);
		typeList.add(type3);
		typeList.add(type4);
		
		PlaneParameters.Isomer ans;
		for (int ind1=0; ind1<4; ind1++){
			for (int ind2=0;ind2<4;ind2++){
				if (ind2==ind1) continue;
				for(int ind3=0;ind3<4;ind3++){
					if((ind3==ind1)||(ind3==ind2)) continue;
					for (int ind4=0;ind4<4;ind4++){
						if((ind4==ind1)||(ind4==ind2)||ind4==ind3) continue;
						if ((ans=planeAngleType(typeList.get(ind1), typeList.get(ind2),
								typeList.get(ind3), typeList.get(ind4), commands))!=null){
							orderedTypes.add(typeList.get(ind1));
							orderedTypes.add(typeList.get(ind2));
							orderedTypes.add(typeList.get(ind3));
							orderedTypes.add(typeList.get(ind4));
							return ans;
							
						}
					}
				}
			}
		}
		return null;
	}*/
	

	
	private static AtomList locatedSublistWithoutGivenAtom(AtomList list, Atom queryAtom) {
		// Tommer 21.12.14
		AtomList locatedBondedAtoms = new AtomList(list.molecularSystem);
		for (Atom atom : list) {
			if (atom.equals(queryAtom)||!atom.residueName().equals(queryAtom.residueName())) 
				//the returned list does not include queryAtom.
				//we only want getAtoms from the same residue. Tommer 4.2.15
				continue;
			if (!atom.nowhere()) {
				locatedBondedAtoms.add(atom);
			}
		}
		return locatedBondedAtoms;
	}
	

	/*
	 * private static AtomList findNeededAtomsForXO_LIKE(Protein model,Atom
	 * atom) { Atom bondedAtom; if (atom.bonded().size()!=1){ throw new
	 * RuntimeException("unloacted atom has an unexpected number of" +
	 * "bonded getAtoms (not 1)"); } bondedAtom=atom.bonded().get(0); if
	 * (bondedAtom.bonded().size()!=3){
	 * System.out.println(bondedAtom.bonded().size()+" bonded getAtoms: \n"+
	 * bondedAtom.bonded().get(0).name+", "+bondedAtom.bonded().get(1).name+
	 * ", atom:"
	 * +atom.name+", atomNumber:"+atom.number()+", numOfAtoms:"+Atom.numberOfAtoms
	 * ()+", protein: "+ model.name()+", ProteinSize"+
	 * model.getAtoms().size()+", residueNumber:"+atom.residueNumber()); throw new
	 * RuntimeException("the atom bonded to the unlocated atom has an" +
	 * " unexpected number of bonded getAtoms (not 3)"); }
	 * 
	 * 
	 * AtomList atomsInPlane=new AtomList(bondedAtom); for(int
	 * atomInd=0;atomInd<3;atomInd++){
	 * if(bondedAtom.bonded().get(atomInd).equals(atom)) continue;
	 * atomsInPlane.add(bondedAtom.bonded().get(atomInd)); } return
	 * atomsInPlane;//three getAtoms that share a plane with given atom are
	 * returned. }
	 * 
	 * private static AtomList findNeededAtomsForFD_LIKE(Protein model,Atom
	 * atom) { Atom bondedAtom, bondedAtom1, bondedAtom2;
	 * 
	 * if (atom.bonded().size()!=2){ throw new
	 * RuntimeException("atom: "+atom+", type: "+atom.type()+
	 * "unloacted atom has an unexpected number of" +
	 * "bonded getAtoms (not 2). \n bonded getAtoms: "+atom.bonded().size()
	 * +" getAtoms, "+atom.bonded()); } bondedAtom1=atom.bonded().atomAt(0);
	 * bondedAtom2=atom.bonded().atomAt(1);
	 * 
	 * AtomLexicographicComparator lexComparator = new AtomFinding().new
	 * AtomLexicographicComparator(); //the above line is WEIRD but i can't seem
	 * to find a way around it... TOMMER 29.11.14
	 * if(lexComparator.compare(bondedAtom1, bondedAtom2)==-1){
	 * bondedAtom=bondedAtom1; } else{ if(lexComparator.compare(bondedAtom1,
	 * bondedAtom2)==1){ bondedAtom=bondedAtom2; } else{ throw new
	 * RuntimeException("error comparing getAtoms by lexicographical order"); } }
	 * 
	 * AtomList atomsInPlane=new AtomList(bondedAtom); for(int
	 * atomInd=0;atomInd<3;atomInd++){
	 * if(bondedAtom.bonded().get(atomInd).equals(atom)) continue;
	 * atomsInPlane.add(bondedAtom.bonded().get(atomInd)); } return
	 * atomsInPlane;//three getAtoms that share a plane with given atom are
	 * returned. }
	 * 
	 * private static AtomList findNeededAtoms(Protein model, Atom atom,
	 * AtomFinding.AtomFindingRoutine routine) { switch (routine){ case XO_LIKE:
	 * return findNeededAtomsForXO_LIKE(model, atom); case SECOND_FD_LIKE:
	 * return findNeededAtomsForFD_LIKE(model, atom); case FIRST_FD_LIKE: return
	 * findNeededAtomsForFD_LIKE(model, atom); default: throw new
	 * RuntimeException("the function 'findNeededAtoms' cannot yet handle" +
	 * "AtomFindingRoutines which are not XO_LIKE, FD_LIKE or YD_LIKE"); } }
	 */

	/*
	 * private static boolean isOfType(Atom atom, AtomFinding.AtomFindingRoutine
	 * routine, Protein model){ ArrayList<Atom> neighbors; switch (routine){
	 * case XO_LIKE: return
	 * atom.backboneH()||atom.backboneO()||(atom.type()==AtomType.NOD)||
	 * (atom.type()==AtomType.QOE)||(atom.type()==AtomType.DO);
	 * 
	 * 
	 * 
	 * case SECOND_FD_LIKE:
	 * if(!((atom.type()==AtomType.FCD)||(atom.type()==AtomType.YCD))) return
	 * false; neighbors =findNeededAtomsForFD_LIKE(model, atom); return
	 * (!neighbors.get(1).nowhere())&&(!neighbors.get(2).nowhere()); //throw new
	 * RuntimeException("this routine ("+routine+") is too specific.\n"+ //
	 * "this method only checks for atomTypes and therefore cannot determine the answer."
	 * );
	 * //isOfType=((atom.type()==AtomType.FCD))||((atom.type()==AtomType.YCD));
	 * //isOfType=(isOfType&&((!neighbors.get(1).nowhere())&&(!neighbors.get(2).
	 * nowhere())));
	 * 
	 * case FIRST_FD_LIKE:
	 * if(!((atom.type()==AtomType.FCD)||(atom.type()==AtomType.YCD))) return
	 * false; neighbors =findNeededAtomsForFD_LIKE(model, atom); return
	 * neighbors.get(1).nowhere()^neighbors.get(2).nowhere(); //throw new
	 * RuntimeException("this routine ("+routine+") is too specific.\n"+ //
	 * "this method only checks for atomTypes and therefore cannot determine the answer."
	 * );
	 * //isOfType=((atom.type()==AtomType.FCD))||((atom.type()==AtomType.YCD));
	 * /
	 * /isOfType=(isOfType&&(neighbors.get(1).nowhere()^neighbors.get(2).nowhere
	 * ())); //this is still problematic. what if the wrong neighbor is
	 * missing?!
	 * 
	 * default: throw new RuntimeException("unhandleable routine "+routine); } }
	 */

	
	
	private static boolean giveCoordinatesFor2BondedAtoms(Atom atom, Atom bonded1,
			Atom bonded2 ,BondParametersList bondParametersList, AngleParametersList angleParametersList, MeshiWriter outWriter){
		AtomList neighbors = locatedSublistWithoutGivenAtom(bonded1.bonded(), atom);
		Coordinates result;
		Atom nearAtom, farAtom, secondNeighborAtom;
		if (neighbors.size()>0){
			nearAtom = bonded1; farAtom=neighbors.get(0); secondNeighborAtom = bonded2;
		}
		else{
			neighbors = locatedSublistWithoutGivenAtom(bonded2.bonded(), atom);
			if (neighbors.size()>0){
				nearAtom = bonded2; farAtom=neighbors.get(0); secondNeighborAtom = bonded1;
			}
			else{
				return false;
			}
		}
		double alpha = getAngle(atom, nearAtom, farAtom, angleParametersList);
		Coordinates secondNeighborCoord = new Coordinates(secondNeighborAtom);
		Coordinates nearCoord = new Coordinates(nearAtom);
		Coordinates farCoord = new Coordinates(farAtom);
		double a = getBondDistance(secondNeighborAtom, atom, bondParametersList);
		double b = getBondDistance(nearAtom, atom, bondParametersList);
		double c = nearCoord.distanceFrom(secondNeighborCoord);
		double topAngle = getAngle(secondNeighborAtom, atom, nearAtom, angleParametersList);
		double c0 = thirdTriangleEdge(a, b, topAngle);
		double d = nearCoord.distanceFrom(farCoord);
		double e = farCoord.distanceFrom(secondNeighborCoord);
		double gamma = angleInTriangle(e, d, c);
		result = atomLocation2Bonded(alpha,
				 gamma, farCoord, nearCoord, secondNeighborCoord, a, b, topAngle);

		if(((c>c0*1.2)||(c<c0/1.2))&&(!atom.residueName().equals("PRO"))) return false;
        if (((c>c0*1.1)||(c<c0*0.9))&&(!atom.residueName().equals("PRO"))) return false;
        double distNear = new Coordinates(result).distanceFrom(new Coordinates(nearAtom));//for check Tommer 31.1.15
        double distNear0 = getBondDistance(atom, nearAtom, bondParametersList);//for check Tommer 31.1.15
        if (((distNear>distNear0*1.1)||(distNear<distNear0*0.9))&&(!atom.residueName().equals("PRO")))
            return false;
        double distSec = new Coordinates(atom).distanceFrom(new Coordinates(secondNeighborAtom));//for check Tommer 31.1.15
        double distSec0 = getBondDistance(atom, nearAtom, bondParametersList);//for check Tommer 31.1.15
        if (((distSec>distSec0*1.1)||(distSec<distSec0*0.9))&&(!atom.residueName().equals("PRO")))
            return false;
        atom.setXYZ(result);
		return true;
	}
	
	private static boolean giveCoordinatesFor1BondedAtom(Atom atom,
			                                             Atom bondedAtom, BondParametersList bondParametersList, AngleParametersList angleParametersList, PlaneParametersList planeParametersList, MeshiWriter outWriter)
			throws IOException {
		// this function tries to give coordinates to an atom with 1 located
		// bonded atom. upon success it
		// returns "true" and upon failure it returns "false". Tommer 12.12.14

		AtomList neighborsOfBonded = bondedAtom.bonded();
        AtomList locatedNeighborsOfBonded = null;
		locatedNeighborsOfBonded = locatedSublistWithoutGivenAtom(neighborsOfBonded, atom);
		int numLocatedNeighborsOfBonded = locatedNeighborsOfBonded.size();

		switch (numLocatedNeighborsOfBonded) {
		case 0:
			if (debug) {
				outWriter.println("the atom "+ atom+
						"has no located NeighborsOfBonded.\n");
			}

			return false;// the atom cannot be located
		case 1://here should come a function that checks CIS/TRANS. Tommer 15.1.15
			if (checkForCISorTRANS(atom, bondedAtom, locatedNeighborsOfBonded.get(0),bondParametersList, angleParametersList, planeParametersList)){
//				outWriter.println("\n\n INNER CASE 1 Atom, CISorTRANS: " + atom +"\n"+atom.core.status()
//						+ "\n bondedAtom: " + bondedAtom
//						+ "\n neighbor of bonded: "
//						+ locatedNeighborsOfBonded.get(0) + "\n");
				double dist = new Coordinates(atom).distanceFrom(new Coordinates(bondedAtom));//for check Tommer 31.1.15
				double dist0 = getBondDistance(atom, bondedAtom, bondParametersList);//for check Tommer 31.1.15
				if ((dist>dist0*1.1)||(dist<dist0*0.9)){
					throw new RuntimeException("weird distance between getAtoms:\n "+
				"atom: "+atom+ "\nbondedAtom: "+bondedAtom+"\ncase: 1 INNER CASE: 1"+
				"\ndistance: "+dist+"  distance0: "+dist0);
				}//check Tommer 31.1.15
				return true;
				//if (!(roundForDEBUG.endsWith("9")||roundForDEBUG.endsWith("10")))throw new RuntimeException("is this reached? 22.1.15 round: "+ roundForDEBUG);
			}
			
			giveCoordinatesFor1NeighborOfBonded(atom, bondedAtom,
					locatedNeighborsOfBonded.get(0), bondParametersList, angleParametersList);

			if (debug) outWriter.println("\n\n INNER CASE 1 Atom: " + atom +"\n"+atom.core.status()
					+ "\n bondedAtom: " + bondedAtom
					+ "\n neighbor of bonded: "
					+ locatedNeighborsOfBonded.get(0) + "\n");

			return true;
		case 2:
			giveCoordinatesFor2NeighborsOfBondedAtom(atom, bondedAtom,
					locatedNeighborsOfBonded.get(0),
					locatedNeighborsOfBonded.get(1), bondParametersList, angleParametersList);
			// if (atom.number()==90){
			if(debug) outWriter.println("\n\n INNER CASE 2 Atom: " + atom
					+"\n"+atom.core.status()
					+ "\n bondedAtom: " + bondedAtom
					+ "\n 1st neighbor of bonded: "
					+ locatedNeighborsOfBonded.get(0)
					+ "\n 2nd neighbor of bonded: "
					+ locatedNeighborsOfBonded.get(1) + "\n");

			double dist = new Coordinates(atom).distanceFrom(new Coordinates(bondedAtom));//for check Tommer 31.1.15
			double dist0 = getBondDistance(atom, bondedAtom, bondParametersList);//for check Tommer 31.1.15
			if ((dist>dist0*1.1)||(dist<dist0*0.9)){
				throw new RuntimeException("weird distance between getAtoms:\n "+
			"atom: "+atom+ "\nbondedAtom: "+bondedAtom+"\ncase: 1 INNER CASE: 2"+
			"\ndistance: "+dist+"  distance0: "+dist0);
			}//check Tommer 31.1.15
			
			return true;
		default:

			throw new RuntimeException(
					"there are too many located neighbors for the bonded atom (more than 2)");
		}

	}

	private static boolean checkForCISorTRANS(Atom atom1, Atom atom2,
			Atom atom3, BondParametersList bondParametersList, AngleParametersList angleParametersList, PlaneParametersList planeParametersList){
		AtomList lastAtoms=locatedSublistWithoutGivenAtom(atom3.bonded(),atom2);
	
		switch(lastAtoms.size()){ 
		case 0:
			return false;
		case 1: 
			PlaneParameters.Isomer angleType1= planeAngleTypeBackAndForth(atom1.type(), atom2.type()
					, atom3.type(),lastAtoms.get(0).type(),planeParametersList);
			if(angleType1==null) return false; //do we need a message here?? Tommer 17.1.15
				
			giveCoordinatesForCISorTRANS(lastAtoms.get(0), atom3, atom2, atom1,bondParametersList, angleParametersList , angleType1);
			return true;
		case 2:
			PlaneParameters.Isomer angleType2= planeAngleTypeBackAndForth(atom1.type(), atom2.type()
					, atom3.type(),lastAtoms.get(0).type(),planeParametersList);
			int index;
			if(angleType2!=null) index=0;
			else{
				angleType2= planeAngleTypeBackAndForth(atom1.type(), atom2.type()
						, atom3.type(),lastAtoms.get(1).type(),planeParametersList);
				if(angleType2!=null) index=1;
				else return false;//do we need a message here? Tommer 17.1.15
	  		}
			angleType2= planeAngleTypeBackAndForth(atom1.type(), atom2.type()
					, atom3.type(),lastAtoms.get(index).type(),planeParametersList);
			
			giveCoordinatesForCISorTRANS(lastAtoms.get(index), atom3, atom2, atom1,  bondParametersList, angleParametersList, angleType2);
			
			return true;

		default:

			throw new RuntimeException(
					"there are too many located last getAtoms (more than 2)");
		}
	}
	
	/*private static boolean checkForImproperTRANS(Atom atom1, Atom atom2,
			Atom middleAtom, Atom atom4ToBePlaced, CommandList commands){
		//getAtoms 1,2 should be the getAtoms not in the middle and with coordinates, regardless of order.
		//atom 3 should be the middle atom. atom 4 should be the atom to be placed.
		if((planeAngleType(atom1.type(), atom2.type(), middleAtom.type(), atom4ToBePlaced.type(), commands)!=null)||
			(planeAngleType(atom2.type(), atom1.type(), middleAtom.type(), atom4ToBePlaced.type(), commands)!=null)){
			return true;
		}
		else{
			return false;
		}
	}*/

	
	

	private static boolean giveCoordinatesForCISorTRANS(Atom atom1, Atom atom2,
			Atom atom3, Atom atom4toBeFound, BondParametersList bondParametersList, AngleParametersList angleParametersList, Isomer CISorTRANS){

		double thirdToFourthDistance= getBondDistance(atom3, atom4toBeFound, bondParametersList);
		Coordinates thirdAtom;
		Coordinates XdirNorm = new Coordinates(atom2).differenceVector
				(thirdAtom= new Coordinates(atom3)).normalize();

		Coordinates firstToSecond= new Coordinates(atom1).differenceVector(new Coordinates(atom2));

		Coordinates YdirNorm= (XdirNorm).vectorMult(firstToSecond.vectorMult(XdirNorm)).normalize(); 

		double angle = getAngle(atom2, atom3, atom4toBeFound, angleParametersList);

		Coordinates thirdToFourthAtom;

		if(CISorTRANS.equals(Isomer.CIS)){
			thirdToFourthAtom = XdirNorm.scalarMult(thirdToFourthDistance*
					Math.cos(PI-angle)).vectorAdd(YdirNorm.scalarMult(-thirdToFourthDistance*
							Math.sin(angle)));
		}
		else{
			if(CISorTRANS.equals(Isomer.TRANS)||CISorTRANS.equals(Isomer.CIS_TRANS)){
				thirdToFourthAtom = XdirNorm.scalarMult(thirdToFourthDistance*
						Math.cos(PI-angle)).vectorAdd(YdirNorm.scalarMult
								(thirdToFourthDistance*Math.sin(angle)));
			}
			else
				throw new RuntimeException("neither CIS nor TRANS nor CIS_TRANS");
		}


		Coordinates fourthCoord = thirdAtom.vectorAdd(thirdToFourthAtom);

		atom4toBeFound.setXYZ(fourthCoord);	
		double dist = new Coordinates(atom3).distanceFrom(new Coordinates(atom4toBeFound));//for check Tommer 31.1.15
		double dist0 = getBondDistance(atom3, atom4toBeFound, bondParametersList);//for check Tommer 31.1.15
		if ((dist>dist0*1.1)||(dist<dist0*0.9)){
			throw new RuntimeException("weird distance between getAtoms:\n "+
		"atom: "+atom4toBeFound+ "\nbondedAtom: "+atom3+"\ncase: 1 INNER CASE: CISorTRANS"+
					"\ndistance: "+dist+"  distnance0: "+dist0+"\nangleType: "+CISorTRANS+
					"\nthirdToFourthAtom norm: "+thirdToFourthAtom.norm()+
					"\nthirdToFourthDistance: "+thirdToFourthDistance);
		}//check Tommer 31.1.15
		return true;//this is somewhat redundant... Tommer 17.1.15
	}
	
	/*private static boolean giveCoordinatesForImproperTRANS(Atom atom1, Atom atom2, Atom middleAtom, Atom atom4ToBePlaced,
			CommandList commands){
		Coordinates middle = new Coordinates (middleAtom);
		Coordinates xNorm = new Coordinates(atom1).differenceVector(middle).normalize();
		Coordinates zNorm = xNorm.vectorMult(new Coordinates(atom2)).differenceVector(middle).normalize();
		Coordinates yNorm = zNorm.vectorMult(xNorm);
		double distToFourth = middleAtom.distanceFrom(atom4ToBePlaced);
		double angle134 = getAngle(atom1,middleAtom, atom4ToBePlaced, commands);
		Coordinates middleToFourth = xNorm.scalarMult(-distToFourth*Math.cos(angle134)).
				vectorAdd(yNorm.scalarMult(distToFourth*Math.sin(angle134)));
		Coordinates fourthCoord = middle.vectorAdd(middleToFourth);
		
		atom4ToBePlaced.setXYZ(fourthCoord);	
		return true;
	}*/
	
	public static void finalFindAtoms(Protein model, BondParametersList bondParametersList,
                                      AngleParametersList angleParametersList, PlaneParametersList planeParametersList) throws IOException {
        MeshiWriter debugWriter = null;
		if (debug)
           debugWriter = new MeshiWriter("debug\\LogFinalFindAtoms.txt");

        AtomList exceptionalProlineAtoms = new AtomList(model.molecularSystem);
		//the exceptional Proline getAtoms will be added in the end
	
		for (Atom atom : model.atoms()) {
			if (!atom.nowhere()) continue;

			if ((atom.type().equals(AtomType.PCG)) &&
                    (atom.bonded().size()!= atom.bonded().located().size())) {
				    //this takes care of Proline ring. Tommer 31.1.15
				    exceptionalProlineAtoms.add(atom);
				    continue;
			}
				
			int numLocatedBonded = 0;
			AtomList bondedAtoms, locatedBondedAtoms = null;


			bondedAtoms = atom.bonded();
			locatedBondedAtoms = locatedSublistWithoutGivenAtom(bondedAtoms, atom);
			numLocatedBonded = locatedBondedAtoms.size();
			
			
			switch (numLocatedBonded) {
			case 0:

                if (debug)
                    debugWriter.println("this atom has no located neigbors: " + atom);

				break;
			case 1:// this is the main case i have so far considered

				giveCoordinatesFor1BondedAtom(atom, locatedBondedAtoms.get(0),
						bondParametersList, angleParametersList, planeParametersList, debugWriter);

                if (debug)
                    debugWriter.println("found case 1 Atom: " + atom
                        + "\n" + atom.core.status()
                        + "\n bondedAtom: " + locatedBondedAtoms.get(0));
				break;
			case 3:
			case 2:
				Atom bondedAtom;
				bondedAtom=locatedBondedAtoms.get(0);
			    if(!giveCoordinatesFor2BondedAtoms(atom, locatedBondedAtoms.get(0),
							locatedBondedAtoms.get(1), bondParametersList, angleParametersList, debugWriter))
						if(!giveCoordinatesFor1BondedAtom(atom, bondedAtom, bondParametersList, angleParametersList,planeParametersList,
								debugWriter)){
							bondedAtom=locatedBondedAtoms.get(1);
							if(!giveCoordinatesFor1BondedAtom(atom, bondedAtom, bondParametersList, angleParametersList, planeParametersList,
									debugWriter)){
								debugWriter.println("found case 2 or 3 Atom, Coordinates NOT GIVEN: " + atom
                                        + "\n bondedAtom: " + bondedAtom
                                        + "\n [this is a temporary message structure]");
								break;
							}
						}

				if (debug) {
                    if (numLocatedBonded == 2)
                        debugWriter.println("\nfound case 2 Atom: " + atom
                                + "\n bondedAtom: " + bondedAtom
                                + "\n [this is a temporary message structure]");
                    else
                        debugWriter.println("\nfound case 3 Atom: " + atom
                                + "\n bondedAtom: " + bondedAtom
                                + "\n [this is a temporary message structure]");
                }
				break;
		

			default:
				throw new RuntimeException(
						"there are too many located bonded getAtoms (more than 3)");
			}

		}
		for (Atom atom : exceptionalProlineAtoms) {
			AtomList locatedBonded = atom.bonded().located();
			if (locatedBonded.size()!=2) 
		        continue;
			if(!giveCoordinatesFor2BondedAtoms(atom, locatedBonded.get(0),
							locatedBonded.get(1), bondParametersList, angleParametersList, debugWriter))
				Utils.println("unlocatable Proline CG (PCG) - location failed"
						+"atom: "+atom);
			if (debug)
                debugWriter.println("\n\nProline PCG added: " + atom
                    + "\n bondedAtom1: " + locatedBonded.get(0) +
                    "\n bondedAtom2: " + locatedBonded.get(1));

		}
		
		if (debug)
            debugWriter.close();

	}
	
	public static void finalFindAtomsIterative(Protein model,
                                               BondParametersList bondParametesList,
                                               AngleParametersList angleParametersList,
                                               PlaneParametersList planeParametersList, int round) throws IOException {
		
		MeshiWriter outWriter = new MeshiWriter("..\\MeshiTatiana\\tommerOut\\"
				+ "getAtoms for routine finalFindAtomsIterative round "+round+".pdb");
		int residueNumber=-1;
	
		for (Atom atom : model.atoms()) {
			if (!atom.nowhere()||atom.residueNumber()==residueNumber) continue;
			
			if ((atom.type().equals(AtomType.PCG))
					&&(atom.bonded().size()!=
					atom.bonded().located().size())) {
				//this takes care of Proline ring. Tommer 31.1.15
				continue;
			}
				
			int numLocatedBonded = 0;
			AtomList bondedAtoms, locatedBondedAtoms = null;


			bondedAtoms = atom.bonded();
			locatedBondedAtoms = locatedSublistWithoutGivenAtom(bondedAtoms, atom);
			numLocatedBonded = locatedBondedAtoms.size();
			
			
			switch (numLocatedBonded) {
			case 0:
		
				outWriter.println("this atom has no located neigbors: " + atom);
				
				break;
			case 1:

				giveCoordinatesFor1BondedAtom(atom, locatedBondedAtoms.get(0),
						bondParametesList, angleParametersList, planeParametersList, outWriter);
				
				outWriter.println("found case 1 Atom: " + atom
						+"\n"+atom.core.status()
						+"\n bondedAtom: " + locatedBondedAtoms.get(0));
				break;
			case 3:
			case 2:
				Atom bondedAtom;
				bondedAtom=locatedBondedAtoms.get(0);
				try{
					if(!giveCoordinatesFor2BondedAtoms(atom, locatedBondedAtoms.get(0),
							locatedBondedAtoms.get(1), bondParametesList, angleParametersList, outWriter))
						if(!giveCoordinatesFor1BondedAtom(atom, bondedAtom, bondParametesList, angleParametersList, planeParametersList,
								outWriter)){
							bondedAtom=locatedBondedAtoms.get(1);
							if(!giveCoordinatesFor1BondedAtom(atom, bondedAtom, bondParametesList, angleParametersList, planeParametersList,
									outWriter)){
								outWriter.println("found case 2 or 3 Atom, Coordinates NOT GIVEN: " + atom
										+ "\n bondedAtom: " +bondedAtom
										+ "\n [this is a temporary message structure]");
								break;
							}
						}

				}catch(Exception ex){
					try{
						model.printAtomsToFile("..\\MeshiTatiana\\tommerOut\\proteinAfterException.pdb");
					}catch(Exception IOex){
						IOex.printStackTrace();
					}
					ex.printStackTrace();
					throw new RuntimeException("program thrown. this is usually due to a "
							+ "geometric problem.");
				}

				if(numLocatedBonded==2)
					outWriter.println("\nfound case 2 Atom: " + atom
							+ "\n bondedAtom: " +bondedAtom
							+ "\n [this is a temporary message structure]");
				else
					outWriter.println("\nfound case 3 Atom: " + atom
							+ "\n bondedAtom: " +bondedAtom
							+ "\n [this is a temporary message structure]");

				break;
		

			default:
				throw new RuntimeException(
						"there are too many located bonded getAtoms (more than 3)");
			}
			residueNumber=atom.residueNumber();

		}
	
		outWriter.close();

	}
	
	

	/*
	 * private static ArrayList<Atom> find1stFDlikeAtoms(Protein model) throws
	 * IOException{ ArrayList<Atom> targetAtoms = new ArrayList<Atom>();
	 * MeshiWriter outWriter = new
	 * MeshiWriter("..\\MeshiTatiana\\tommerOut\\"+"getAtoms for routine YD-LIKE .pdb"
	 * ); boolean residueVisited=false;//this is to make sure we only take up to
	 * 1 atom per residue, Tommer 12.12.14 int residueNumber=0; for(Atom atom :
	 * model.getAtoms()){ if(residueVisited){
	 * if(atom.residueNumber()!=residueNumber){ residueVisited=false; } else{
	 * continue; } } if(isOfType(atom,
	 * AtomFinding.AtomFindingRoutine.FIRST_FD_LIKE, model) && atom.nowhere()){
	 * ArrayList<Atom> neighbors = findNeededAtoms(model, atom,
	 * AtomFinding.AtomFindingRoutine.FIRST_FD_LIKE); try{ boolean
	 * locatedNeighbors=true;
	 * 
	 * boolean firstAtom=true; if (neighbors.get(0).nowhere()){
	 * outWriter.print(atom
	 * .toString()+"	closest used neighbor: "+neighbors.get(0)
	 * +"the coordinates for this neighbor are unknown.\n");
	 * locatedNeighbors=false; } AtomLexicographicComparator lexComparator = new
	 * AtomFinding().new AtomLexicographicComparator(); //the above line is
	 * WEIRD but i can't seem to find a way around it... TOMMER 29.11.14
	 * 
	 * for (int ind=1;ind<3;ind++){//THIS NEEDS MODIFICATION! TOMMER 22.11.14 if
	 * (neighbors.get(ind).nowhere() &&
	 * lexComparator.compare(neighbors.get(ind), atom)<0){//there is a problem
	 * outWriter.print(atom.toString()+"	neighbor: "+neighbors.get(ind)
	 * +"the coordinates for this neighbor are unknown. Routine YD-LIKE.\n");
	 * locatedNeighbors=false; } if (!neighbors.get(ind).nowhere() &&
	 * lexComparator.compare(neighbors.get(ind), atom)==0){ firstAtom=false; }
	 * 
	 * } if(locatedNeighbors && firstAtom){ //this is the main case
	 * outWriter.print
	 * ("YD_LIKE: "+atom.toString()+"\n neighbor: "+neighbors.get(0)+
	 * "\n neighbor: "+neighbors.get(1)+"\n neighbor: "+neighbors.get(2)+"\n");
	 * targetAtoms.add(atom); //here the atom is added to the designated list
	 * residueNumber=atom.residueNumber(); residueVisited=true;
	 * System.out.println("found YD_LIKE Atom: "+ atom); } if(locatedNeighbors
	 * && !firstAtom){ //this case is like the previous cases. we do not need to
	 * take it to the HD List. // } firstAtom=true; locatedNeighbors=true;
	 * 
	 * 
	 * } catch(RuntimeException ex){ if(ex.getMessage()==null){
	 * outWriter.close(); throw ex; } if(ex.getMessage().contains("(not 3)")){
	 * outWriter.print(atom.toString()+"this  atom does not have 3 neighbors. "+
	 * "this usually indicates that it is in the first residue. residue number: "
	 * + atom.residueNumber()+"\n"); continue; } else { outWriter.close(); throw
	 * ex; } } }
	 * 
	 * } outWriter.close(); //throw new
	 * RuntimeException("this is the end of findYDlikeAtoms"); return
	 * targetAtoms; }
	 * 
	 * private static ArrayList<Atom> findXO_LIKEor2ndFDlikeAtoms (Protein
	 * model, AtomFinding.AtomFindingRoutine routine) throws IOException{
	 * ArrayList<Atom> targetAtoms = new ArrayList<Atom>(); MeshiWriter
	 * outWriter = new MeshiWriter("..\\MeshiTatiana\\tommerOut\\"+
	 * "getAtoms for routine basicRoutine.pdb"); for(Atom atom : model.getAtoms()){
	 * if(isOfType(atom, routine, model) && atom.nowhere()){//it looks like the
	 * "nowhere" check is aptly placed. Tommer 29.11.14 ArrayList<Atom>
	 * neighbors = findNeededAtoms(model, atom, routine); try{ boolean
	 * locatedNeighbors=true;
	 * 
	 * for (int ind=0;ind<3;ind++){//THIS IS SEMI-REPLICATED in another
	 * function. TOMMER 22.11.14 if (neighbors.get(ind).nowhere()){
	 * outWriter.print(atom.toString()+"	neighbor: "+neighbors.get(ind)
	 * +"the coordinates for this neighbor are unknown.\n");
	 * locatedNeighbors=false; } } if(locatedNeighbors){//this is the desired
	 * case
	 * outWriter.print("XOorFDlike: "+atom.toString()+"\n neighbor: "+neighbors
	 * .get(0)+
	 * "\n neighbor: "+neighbors.get(1)+"\n neighbor: "+neighbors.get(2)+"\n");
	 * targetAtoms.add(atom); //here the atom is added to the designated list }
	 * locatedNeighbors=true;//reset for next iteration
	 * 
	 * } catch(RuntimeException ex){ if(ex.getMessage()==null){
	 * outWriter.close(); throw ex; } if(ex.getMessage().contains("(not 3)")){
	 * outWriter.print(atom.toString()+"this  atom does not have 3 neighbors. "+
	 * "this usually indicates that it is in the first residue. residue number: "
	 * + atom.residueNumber()+"\n"); continue; } else { outWriter.close(); throw
	 * ex; } }
	 * 
	 * } } outWriter.close(); return targetAtoms; }//this method only finds
	 * getAtoms that have 3 neighbors with known coordinates and can therefore be
	 * placed.
	 */

	private static Coordinates getCoordinates(Atom atom) {// Tommer 23.10.14
		return new Coordinates(atom.x(), atom.y(), atom.z());
	}

	private static void giveCoordinatesFor1NeighborOfBonded(Atom atom,
			Atom bondedAtom, Atom farAtom, BondParametersList bondParametersList, AngleParametersList angleParametersList) {
		// this is the case where we only have two located getAtoms to use

		if (bondedAtom.nowhere() || farAtom.nowhere()) {
			System.out.println("atom: " + atom + "\n bondedAtom: " + bondedAtom
					+ "\n farAtom: " + farAtom);
			throw new RuntimeException("(at least) one of the getAtoms is not located");
		}
		
		double nearAtomToDestinationDistance = AtomFinding.getBondDistance(atom, bondedAtom, bondParametersList);
		double angle = getAngle(atom, bondedAtom, farAtom, angleParametersList);
		// this order is fine. Tommer 3.1.15
		
		/*
		 * next the distance should be decomposed into parallel and
		 * perpendicular components. the parallel component is easily added and
		 * the perpendicular should be added at random to one of the possible
		 * directions.
		 */
		
		Coordinates nearCoords = getCoordinates(bondedAtom);
		Coordinates farCoords = getCoordinates(farAtom);
		Coordinates axisDirNorm = farCoords.differenceVector(nearCoords).normalize();//zNorm
		// right direction. Tommer 3.1.15

		double parallelDist = nearAtomToDestinationDistance
				* (-Math.cos(angle));// Distance to the middle of the circle
		double radius = nearAtomToDestinationDistance * (Math.sin(angle));
		//circle radius

		Coordinates parallelVector = axisDirNorm.scalarMult(parallelDist);

		Coordinates middleOfCircle = nearCoords.vectorAdd(parallelVector);
		
		Coordinates xNorm = new Coordinates(
				axisDirNorm.y(), -axisDirNorm.x(), 0).normalize();
		// the perpendicular plane is depicted with two axes: x, y.
		Coordinates yNorm = axisDirNorm.vectorMult(xNorm);



		/*if(isFirstTry){//This does not work... Tommer 3.1.15
			Coordinates smallX_AxisVec=xNorm.scalarMult(radius/100);
			atom.setXYZ(middleOfCircle.x()+smallX_AxisVec.x(),
					middleOfCircle.y()+smallX_AxisVec.y(),
					middleOfCircle.z()+smallX_AxisVec.z());
			return;
		}*/
		//MeshiProgram.initRandom(0);//Necessary only when "working SOLO"! 		
		//Tommer 28.12.14
		try{//This is fishy programming... Tommer 28.12
			Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();
	
			double phi = randomNumberGenerator.nextDouble() * 2 * PI;
			
			Coordinates perpendicularVector = xNorm.scalarMult(Math
					.cos(phi) * radius);//this is the x projection
			perpendicularVector = perpendicularVector.vectorAdd(yNorm
					.scalarMult(Math.sin(phi) * radius)); //this adds the y projection

			atom.setXYZ(middleOfCircle.x() + perpendicularVector.x(),
					middleOfCircle.y() + perpendicularVector.y(),
					middleOfCircle.z() + perpendicularVector.z());
			
		}catch(RuntimeException e){
				MeshiProgram.initRandom(0);
				Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

				double phi = randomNumberGenerator.nextDouble() * 2 * PI;
				Coordinates perpendicularVector = xNorm.scalarMult(Math
						.cos(phi) * radius);//this is the x projection
				perpendicularVector = perpendicularVector.vectorAdd(yNorm
						.scalarMult(Math.sin(phi) * radius)); //this adds the y projection

				atom.setXYZ(middleOfCircle.x() + perpendicularVector.x(),
						middleOfCircle.y() + perpendicularVector.y(),
						middleOfCircle.z() + perpendicularVector.z());
				
			}
	}

	private static void giveCoordinatesFor2NeighborsOfBondedAtom(Atom atom,
			                                                     Atom triangleTipAtom,
                                                                 Atom triangleBaseAtom1,
			                                                    Atom triangleBaseAtom2,
                                                                BondParametersList bondParametersList,
                                                                AngleParametersList angleParametersList) {
		if (triangleTipAtom.nowhere() || triangleBaseAtom1.nowhere()
				|| triangleBaseAtom2.nowhere())
			throw new RuntimeException("one (or more) of the getAtoms is not located");

		Coordinates triangleTipCoord = getCoordinates(triangleTipAtom);
		Coordinates triangleBaseCoord1 = getCoordinates(triangleBaseAtom1);
		Coordinates triangleBaseCoord2 = getCoordinates(triangleBaseAtom2);
		double alpha = getAngle(atom, triangleTipAtom, triangleBaseAtom1,angleParametersList);
		double beta = getAngle(atom, triangleTipAtom, triangleBaseAtom2,angleParametersList);
		double gamma = getAngle(triangleBaseAtom1, triangleTipAtom, triangleBaseAtom2, angleParametersList);
		double sinTheta = Coordinates.sineOfOut_of_PlaneAngle(alpha, beta, gamma);
		double gamma2 = Math.acos(Coordinates.cosineOfgamma2(sinTheta, alpha));
		double distance = getBondDistance(atom, triangleTipAtom, bondParametersList);

		Coordinates vec1 = triangleBaseCoord1
				.differenceVector(triangleTipCoord);
		Coordinates vec2 = triangleBaseCoord2
				.differenceVector(triangleTipCoord);
		Coordinates normalX = vec1.normalize(), normalZ = normalX.vectorMult(vec2).normalize();
		Coordinates normalY= normalZ.vectorMult(normalX);//changed 31.1.15, there were geometric errors.
		// this direction seems to yield better results...
		// Still Indeterminate! Tommer 28.12.14 nearly no effect on sample...

		//if(atom.residueNumber()==65)// CHECK Tommer 28.12.14
			//throw new RuntimeException("residue 65 by giveCoordinatesFor2NeighborsOfBondedAtom");
		double dX = -distance * Math.cos(alpha);
		double dY = dX*Math.tan(gamma2);
		double dZ = distance*sinTheta;
		
		Coordinates designatedAtomCoords = triangleTipCoord.vectorAdd(normalX
				.scalarMult(dX).vectorAdd(normalY.scalarMult(dY)
						.vectorAdd(normalZ.scalarMult(dZ))));

		if (Double.isNaN(designatedAtomCoords.x())) {
			System.out.println("normal 1: " + normalX + "\n normalY: "
					+ normalY + "\n normalZ: " + normalZ
					+ "\n alpha: " + alpha + ", beta: " + beta + ", gamma: "
					+ gamma + ", distance: " + distance + "\n sinTheta: "
					+ sinTheta);
			throw new RuntimeException(
					"we have a NaN value. this is not supposed to happen.");
		}
		atom.setXYZ(designatedAtomCoords.x(), designatedAtomCoords.y(),
				designatedAtomCoords.z());
	}
	
	public static Protein copyOfProtein(Protein templateProtein){
		AtomList tempList= new AtomList(templateProtein.atoms().get(0));
		for(int index=1;index<templateProtein.atoms().size();index++){
			tempList.add(templateProtein.atoms().get(index));
		}
		Protein resProtein= new Protein(tempList, new ResidueExtendedAtomsCreator());
		return resProtein;
	}
	
}



/*
 * private static void giveCoordinatesTo1stFD_LIKE(Atom atom, Atom nearAtom,
 * Atom farAtom, Atom farAtom2, CommandList commands){//this is the case where
 * we only have two located neighbors AtomLexicographicComparator lexComparator
 * = new AtomFinding().new AtomLexicographicComparator(); //the above line is
 * WEIRD but i can't seem to find a way around it... TOMMER 29.11.14
 * if(lexComparator.compare(farAtom, farAtom2)>0){ farAtom=farAtom2; } else{
 * if(lexComparator.compare(farAtom, farAtom2)==0){ throw new
 * RuntimeException("the two farAtoms are equal by letter. "+
 * "this is not supposed to happen."); } }
 * 
 * double nearAtomToDestinationDistance=AtomFinding.getBondDistance(atom,
 * nearAtom, commands); double angle=getAngle(atom, nearAtom, farAtom,
 * commands);//Descending order? Tommer 30.11.14
 * 
 * //* next the distance should be decomposed into parallel and perpendicular
 * components. // * the parallel component is easily added and the perpendicular
 * should be added // * at random to one of the possible directions. // *
 * Coordinates nearCoords= getCoordinates(nearAtom); Coordinates farCoords=
 * getCoordinates(farAtom); Coordinates axisDirVec =
 * farCoords.differenceVector(nearCoords);//SHOULD BE right direction, Tommer
 * 30.11.14
 * 
 * 
 * double parallelDist=nearAtomToDestinationDistance*(-Math.cos(angle)); double
 * perpendicularDist=nearAtomToDestinationDistance*(Math.sin(angle));
 * 
 * Coordinates parallelVector = new
 * Coordinates(axisDirVec.x()*parallelDist/axisDirVec.norm(),
 * axisDirVec.y()*parallelDist/axisDirVec.norm(),
 * axisDirVec.z()*parallelDist/axisDirVec.norm());
 * 
 * 
 * Coordinates middleOfCircle=nearCoords.vectorAdd(parallelVector); Coordinates
 * normalZ_AxisVec=parallelVector.normalize(); Coordinates normalX_AxisVec= new
 * Coordinates( normalZ_AxisVec.y()/(Math.pow(normalZ_AxisVec.x(),
 * 2)+Math.pow(normalZ_AxisVec.y(), 2)),
 * -normalZ_AxisVec.x()/(Math.pow(normalZ_AxisVec.x(),
 * 2)+Math.pow(normalZ_AxisVec.y(), 2)),0); Coordinates normalY_AxisVec =
 * normalZ_AxisVec.vectorMult(normalX_AxisVec); Random randomNumberGenerator =
 * MeshiProgram.randomNumberGenerator(); double phi =
 * randomNumberGenerator.nextDouble()*2*PI; Coordinates perpendicularVector =
 * normalX_AxisVec.scalarMult(Math.cos(phi)*perpendicularDist);
 * perpendicularVector
 * =perpendicularVector.vectorAdd(normalY_AxisVec).scalarMult(
 * Math.sin(phi)*perpendicularDist);
 * 
 * atom.setXYZ(middleOfCircle.x()+perpendicularVector.x(),middleOfCircle.y()+
 * perpendicularVector.y(), middleOfCircle.z()+perpendicularVector.z());
 * 
 * }
 * 
 * 
 * private static void giveCoordinatesToXO_LIKEor2ndFD_LIKE(Atom atom, Atom
 * triangleTipAtom, Atom triangleBaseAtom1, Atom triangleBaseAtom2, CommandList
 * commands){ Coordinates triangleTipCoord= getCoordinates(triangleTipAtom);
 * Coordinates triangleBaseCoord1= getCoordinates(triangleBaseAtom1);
 * Coordinates triangleBaseCoord2= getCoordinates(triangleBaseAtom2);
 * 
 * Coordinates Mcoordinates =
 * triangleBaseCoord1.MiddleCoordinates(triangleBaseCoord2); Coordinates
 * dirVecCoordinates = Mcoordinates.differenceVector(triangleTipCoord);
 * 
 * double
 * middleToSecondAtomDistance=Mcoordinates.distanceFrom(triangleTipCoord);
 * double secondAtomToDestinationDistance=AtomFinding.getBondDistance(atom,
 * triangleTipAtom, commands); double
 * multParam=middleToSecondAtomDistance/secondAtomToDestinationDistance
 * ;//parameter to multiply vector
 * 
 * atom.setXYZ(triangleTipCoord.x()+dirVecCoordinates.x()*multParam,
 * triangleTipCoord.y()+dirVecCoordinates.y()*multParam,
 * triangleTipCoord.z()+dirVecCoordinates.z()*multParam); }
 * 
 * public static void giveCoordinates(Protein model, ArrayList<Atom>
 * atomsToBePlaced,AtomFindingRoutine routine, CommandList commands) throws
 * IOException { MeshiWriter outWriter = new
 * MeshiWriter("..\\MeshiTatiana\\tommerOut\\"+routine.toString()+".pdb");
 * for(Atom atom: atomsToBePlaced){//Tommer changed variable name 22.10.14 try{
 * AtomList neededAtoms = AtomFinding.findNeededAtoms( model, atom, routine);
 * ///
 * /////////////////////////////////////////////////////////////////////////////
 * //// H 1. Find M - middle of C-Ca // | 2. Find the HM line equation // | 3.
 * Find H coordinates outside CNCa triangle // N |NH|=0.92 // /\ // / \ // / \
 * // C CA // // // Oh // | Tyrosine // | 1. Find M - middle of CE1-CE2 // CZ 2.
 * Find the OM line equation // /\ 3. Find H coordinates outside CE1CzCE2
 * triangle // / \ |CzO|=1.37 // CE1 CE2 // | | // | | // \ / // \ / // | // R
 * //
 * ////////////////////////////////////////////////////////////////////////////
 * ///
 * 
 * Atom triangleTipAtom=neededAtoms.get(0);//this is the atom at the tip of the
 * triangle Atom triangleBaseAtom1 = neededAtoms.get(1) ,
 * triangleBaseAtom2=neededAtoms.get(2);//distinction between the remaining
 * getAtoms is irrelevant
 * 
 * 
 * if
 * ((triangleBaseAtom1==null)||(triangleTipAtom==null)||(triangleTipAtom==null
 * )){ throw new
 * RuntimeException("error recognizing getAtoms: uninitialized atom found (null)");
 * } if ((triangleBaseAtom1.equals(triangleTipAtom))||(triangleTipAtom.equals(
 * triangleBaseAtom2))|| (triangleTipAtom.equals(triangleBaseAtom2))){ throw new
 * RuntimeException("error recognizing getAtoms: identical getAtoms assigned " +
 * "more than once"); }
 * 
 * switch (routine){ case XO_LIKE: giveCoordinatesToXO_LIKEor2ndFD_LIKE(atom,
 * triangleTipAtom, triangleBaseAtom1, triangleBaseAtom2, commands); break; case
 * SECOND_FD_LIKE: giveCoordinatesToXO_LIKEor2ndFD_LIKE(atom, triangleTipAtom,
 * triangleBaseAtom1, triangleBaseAtom2, commands); break; case
 * FIRST_FD_LIKE://this is the case where we only have two located neighbors
 * giveCoordinatesTo1stFD_LIKE(atom, triangleTipAtom, triangleBaseAtom1,
 * triangleBaseAtom2, commands); break; default: throw new
 * RuntimeException("can't give coordinates to getAtoms with this AtomFindingRoutine: "
 * + routine);
 * 
 * } } catch(NullPointerException ex){ ex.printStackTrace(); outWriter.close();
 * throw new RuntimeException("giveCoordinates has failed");
 * 
 * } }
 * 
 * 
 * for(Atom atom : atomsToBePlaced){ outWriter.print(atom+"\n"); }
 * outWriter.close(); }
 * 
 * public static ArrayList<Atom> findAtomsOfType (Protein model,
 * AtomFinding.AtomFindingRoutine routine) throws IOException{ if
 * (routine.equals
 * (AtomFinding.AtomFindingRoutine.XO_LIKE)||routine.equals(AtomFinding
 * .AtomFindingRoutine.SECOND_FD_LIKE)) return
 * findXO_LIKEor2ndFDlikeAtoms(model, routine); if
 * (routine.equals(AtomFinding.AtomFindingRoutine.FIRST_FD_LIKE)) return
 * find1stFDlikeAtoms(model); throw new RuntimeException(
 * "this function can only take care of the BASIC or YD-LIKE AtomFindingRoutine types."
 * ); }
 */

