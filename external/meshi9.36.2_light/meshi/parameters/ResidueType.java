package meshi.parameters;
public enum ResidueType {
    ALA("ALA","A",AtomType.AH,AtomType.AN,AtomType.ACA,AtomType.AC,AtomType.AO,AtomType.ACB,true,null,true,115),     //0
	CYS("CYS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,135), //1
	ASP("ASP","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true,150), //2
	GLU("GLU","E",AtomType.EH,AtomType.EN,AtomType.ECA,AtomType.EC,AtomType.EO,AtomType.ECB,true,190.0), //3
	PHE("PHE","F",AtomType.FH,AtomType.FN,AtomType.FCA,AtomType.FC,AtomType.FO,AtomType.FCB,true,null,true,210), //4
	GLY("GLY","G",AtomType.GH,AtomType.GN,AtomType.GCA,AtomType.GC,AtomType.GO,null        ,true,75), //5
	HIS("HIS","H",AtomType.HH,AtomType.HN,AtomType.HCA,AtomType.HC,AtomType.HO,AtomType.HCB,true,null,true,195), //6
	ILE("ILE","I",AtomType.IH,AtomType.IN,AtomType.ICA,AtomType.IC,AtomType.IO,AtomType.ICB,true,null,true,175), //7
	LYS("LYS","K",AtomType.KH,AtomType.KN,AtomType.KCA,AtomType.KC,AtomType.KO,AtomType.KCB,true,200), //8
	LEU("LEU","L",AtomType.LH,AtomType.LN,AtomType.LCA,AtomType.LC,AtomType.LO,AtomType.LCB,true,null,true,170), //9
	MET("MET","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true,"MSE",true,185), //10
	ASN("ASN","N",AtomType.NH,AtomType.NN,AtomType.NCA,AtomType.NC,AtomType.NO,AtomType.NCB,true,160), //11
	PRO("PRO","P",null,       AtomType.PN,AtomType.PCA,AtomType.PC,AtomType.PO,AtomType.PCB,true,145), //12
	GLN("GLN","Q",AtomType.QH,AtomType.QN,AtomType.QCA,AtomType.QC,AtomType.QO,AtomType.QCB,true,180),
	ARG("ARG","R",AtomType.RH,AtomType.RN,AtomType.RCA,AtomType.RC,AtomType.RO,AtomType.RCB,true,225),
	SER("SER","S",AtomType.SH,AtomType.SN,AtomType.SCA,AtomType.SC,AtomType.SO,AtomType.SCB,true,115),
	THR("THR","T",AtomType.TH,AtomType.TN,AtomType.TCA,AtomType.TC,AtomType.TO,AtomType.TCB,true,140),
	VAL("VAL","V",AtomType.VH,AtomType.VN,AtomType.VCA,AtomType.VC,AtomType.VO,AtomType.VCB,true,null,true,155),
	TRP("TRP","W",AtomType.WH,AtomType.WN,AtomType.WCA,AtomType.WC,AtomType.WO,AtomType.WCB,true,null,true,255),
	TYR("TYR","Y",AtomType.YH,AtomType.YN,AtomType.YCA,AtomType.YC,AtomType.YO,AtomType.YCB,true,null,true,230),
	    // Modified residues often found in PDB files
    // chen 25/3/12 this residue type is handled as an alternative name to MET
	//MSE("MSE","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true), //selenium methionine
	//
	YMC("YMC","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true,-1), //??
	SLN("SLN","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true,-1), //selenium methionine
 	CSE("CSE","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //selenium cysteine
 	CSW("CSW","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //CYSTEINE-S-DIOXIDE
 	CSZ("CSZ","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //selenium cysteine
 	OCS("OCS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //CYSTEINESULFONIC ACID
 	SMC("SMC","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //CYSTEINE RESIDUE METHYLATED IN S-POSITION
 	SYG("SYG","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //GLUTAMYL-S-CYSTEINE
 	CME("CME","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //S,S-(2-HYDROXYETHYL)THIOCYSTEINE
 	CCS("CCS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //????
 	CSX("CSX","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true,-1), //????
        SEP("SEP","S",AtomType.SH,AtomType.SN,AtomType.SCA,AtomType.SC,AtomType.SO,AtomType.SCB,true,-1), //PHOSPHOSERINE
	ASQ("ASQ","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true,-1), //PHOSPHOASPARTATE
	PHD("PHD","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true,-1), //PHOSPHORYLATION
	PCA("PCA","E",AtomType.EH,AtomType.EN,AtomType.ECA,AtomType.EC,AtomType.EO,AtomType.ECB,true,-1), //PYROGLUTAMIC ACID
	ARO("ARO","R",AtomType.RH,AtomType.RN,AtomType.RCA,AtomType.RC,AtomType.RO,AtomType.RCB,true,-1), //???
	HHS("HHS","H",AtomType.HH,AtomType.HN,AtomType.HCA,AtomType.HC,AtomType.HO,AtomType.HCB,true,-1), //???
	DMY(" - ","-",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false,-1), /* a dummy residue used
													   as a place-holder 
													   or as a search key.*/
    CREATOR("XXX","X",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false,-1),
    // for helium cluster minimization
    HEL("HEL","H",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false,-1);

    private final String nameThreeLetters;
    private final String nameOneLetter;
    private final AtomType caType,hType,nType,cType,oType,cbType;
    private final String[] alternativeNames;
    public final boolean isHydrophobic;
    /*
       in       square Angsroms
       From     http://prowl.rockefeller.edu/aainfo/volume.htm
       Based on C. Chotia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14.
     */
    public final double maxExposedSurfaceArea;

    private final boolean standard;
    private ResidueType(String nameThreeLetters, String nameOneLetter,
			AtomType hType, AtomType nType, AtomType caType,
			AtomType cType, AtomType oType, AtomType cbType, boolean standard, double maxExposedSurfaceArea) {
        this(nameThreeLetters, nameOneLetter,
			hType, nType, caType,
			cType, oType, cbType, standard,null,false, maxExposedSurfaceArea);
    }
    private ResidueType(String nameThreeLetters, String nameOneLetter,
			AtomType hType, AtomType nType, AtomType caType, 
			AtomType cType, AtomType oType, AtomType cbType, boolean standard,
            String alternativeName,boolean isHydrophobic, double maxExposedSurfaceArea) {
	    this.nameThreeLetters = nameThreeLetters;
	    this.nameOneLetter = nameOneLetter;
	    this.hType  = hType;
	    this.nType  = nType;
	    this.caType = caType;
	    this.cType  = cType;
	    this.oType  = oType;
	    this.cbType = cbType;
	    this.standard = standard;
        this.isHydrophobic = isHydrophobic;
        if (alternativeName == null)this.alternativeNames = null;
        else {
            alternativeNames = new String[1];
            alternativeNames[0] = alternativeName;
        }
        this.maxExposedSurfaceArea = maxExposedSurfaceArea;
    }
    public String nameThreeLetters() {return nameThreeLetters;}
    public String nameOneLetter() {return nameOneLetter;}
    public String[] alternateNames() {return alternativeNames;}
    public AtomType hType()  {return hType;}
    public AtomType nType()  {return nType;}
    public AtomType caType() {return caType;}
    public AtomType cType()  {return cType;}
    public AtomType oType()  {return oType;}
    public AtomType cbType() {return cbType;}
    public boolean standard() {return standard;}

    public static ResidueType type(String name) {
	if (name.equalsIgnoreCase("CREATOR"))
        return ResidueType.CREATOR;
    if ((name.length() != 3) & (name.length() != 1))
	    return DMY;
	for (ResidueType type:ResidueType.values()) {
	    if ((name.equals(type.nameThreeLetters)) |
            (name.equals(type.nameOneLetter))) return type;
        if (type.alternativeNames != null) {       // chen 25/3/12
            for (String an: type.alternativeNames)
                if (name.equals(an)) return type;
        }
	}
	return DMY;
    }
    public static ResidueType type(char name) {
	for (ResidueType type:ResidueType.values())
	    if (name == type.nameOneLetter.charAt(0)) return type;
	return DMY;
    }
    
    public static ResidueType type(AtomType atomType) {
	for (ResidueType type:ResidueType.values()) {
	    if (atomType.equals(type.hType))  return type;
	    if (atomType.equals(type.nType))  return type;
	    if (atomType.equals(type.caType)) return type;
	    if (atomType.equals(type.cType))  return type;
	    if (atomType.equals(type.oType))  return type;
	    if (atomType.equals(type.cbType)) return type;
	}
	return DMY;
    }

    public static String typeOneLetter(String name) {
        if (name.equalsIgnoreCase("CREATOR"))
        return ResidueType.CREATOR.nameThreeLetters;
    if ((name.length() != 3))
            //throw new RuntimeException("Weird Residue type "+name+" of length "+name.length()+".");
        return null;
        for (ResidueType type:ResidueType.values()) {
            if (name.equals(type.nameThreeLetters)) return type.nameOneLetter;
        }
        return DMY.nameOneLetter;
    }

	public static boolean isPolar(ResidueType type) {
		if (type == SER) return true;
		if (type == THR) return true;
		if (type == ASN) return true;
		if (type == GLN) return true;
		if (type == CYS) return true;
		return false;
	}

	public static boolean isPositive(ResidueType type) {
		if (type == ARG) return true;
		if (type == LYS) return true;
		if (type == HIS) return true;
		return false;
	}

	public static boolean isNegative(ResidueType type) {
		if (type == ASP) return true;
		if (type == GLU) return true;
		return false;
	}

	public static boolean isAliphatic(ResidueType type) {
		if (type == GLY) return true;
		if (type == ALA) return true;
		if (type == ILE) return true;
		if (type == LEU) return true;
		if (type == MET) return true;
		if (type == PRO) return true;
		if (type == VAL) return true;

		return false;
	}
	public static boolean isAromatic(ResidueType type) {
		if (type == PHE) return true;
		if (type == TYR) return true;
		if (type == TRP) return true;
		return false;
	}

	public static ResidueType residueType(AtomType atomType) {
		return ResidueType.type(atomType.toString().charAt(0));
	}
	public static boolean isPolar(AtomType atomType) {
		return isPolar(residueType(atomType));
	}
	public static boolean isNegative(AtomType atomType) {
		return isNegative(residueType(atomType));
	}
	public static boolean isPositive(AtomType atomType) {
		return isPositive(residueType(atomType));
	}
	public static boolean isAromatic(AtomType atomType) {
		return isAromatic(residueType(atomType));
	}
	public static boolean isAliphatic(AtomType atomType) {
		return isAliphatic(residueType(atomType));
	}



}
			     

  
