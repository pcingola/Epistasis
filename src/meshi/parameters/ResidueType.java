package meshi.parameters;
public enum ResidueType {
    ALA("ALA","A",AtomType.AH,AtomType.AN,AtomType.ACA,AtomType.AC,AtomType.AO,AtomType.ACB,true),     //0
	CYS("CYS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //1
	ASP("ASP","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true), //2
	GLU("GLU","E",AtomType.EH,AtomType.EN,AtomType.ECA,AtomType.EC,AtomType.EO,AtomType.ECB,true), //3
	PHE("PHE","F",AtomType.FH,AtomType.FN,AtomType.FCA,AtomType.FC,AtomType.FO,AtomType.FCB,true), //4
	GLY("GLY","G",AtomType.GH,AtomType.GN,AtomType.GCA,AtomType.GC,AtomType.GO,null        ,true), //5
	HIS("HIS","H",AtomType.HH,AtomType.HN,AtomType.HCA,AtomType.HC,AtomType.HO,AtomType.HCB,true), //6
	ILE("ILE","I",AtomType.IH,AtomType.IN,AtomType.ICA,AtomType.IC,AtomType.IO,AtomType.ICB,true), //7
	LYS("LYS","K",AtomType.KH,AtomType.KN,AtomType.KCA,AtomType.KC,AtomType.KO,AtomType.KCB,true), //8
	LEU("LEU","L",AtomType.LH,AtomType.LN,AtomType.LCA,AtomType.LC,AtomType.LO,AtomType.LCB,true), //9
	MET("MET","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true), //10
	ASN("ASN","N",AtomType.NH,AtomType.NN,AtomType.NCA,AtomType.NC,AtomType.NO,AtomType.NCB,true), //11
	PRO("PRO","P",null,       AtomType.PN,AtomType.PCA,AtomType.PC,AtomType.PO,AtomType.PCB,true), //12
	GLN("GLN","Q",AtomType.QH,AtomType.QN,AtomType.QCA,AtomType.QC,AtomType.QO,AtomType.QCB,true),
	ARG("ARG","R",AtomType.RH,AtomType.RN,AtomType.RCA,AtomType.RC,AtomType.RO,AtomType.RCB,true),
	SER("SER","S",AtomType.SH,AtomType.SN,AtomType.SCA,AtomType.SC,AtomType.SO,AtomType.SCB,true),
	THR("THR","T",AtomType.TH,AtomType.TN,AtomType.TCA,AtomType.TC,AtomType.TO,AtomType.TCB,true),
	VAL("VAL","V",AtomType.VH,AtomType.VN,AtomType.VCA,AtomType.VC,AtomType.VO,AtomType.VCB,true),
	TRP("TRP","W",AtomType.WH,AtomType.WN,AtomType.WCA,AtomType.WC,AtomType.WO,AtomType.WCB,true),
	TYR("TYR","Y",AtomType.YH,AtomType.YN,AtomType.YCA,AtomType.YC,AtomType.YO,AtomType.YCB,true),
	    // Modified residues often found in PDB files
	MSE("MSE","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true), //selenium methionine
	YMC("YMC","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true), //??
	SLN("SLN","M",AtomType.MH,AtomType.MN,AtomType.MCA,AtomType.MC,AtomType.MO,AtomType.MCB,true), //selenium methionine
 	CSE("CSE","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //selenium cysteine
 	CSW("CSW","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //CYSTEINE-S-DIOXIDE
 	CSZ("CSZ","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //selenium cysteine
 	OCS("OCS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //CYSTEINESULFONIC ACID
 	SMC("SMC","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //CYSTEINE RESIDUE METHYLATED IN S-POSITION
 	SYG("SYG","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //GLUTAMYL-S-CYSTEINE
 	CME("CME","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //S,S-(2-HYDROXYETHYL)THIOCYSTEINE
 	CCS("CCS","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //????
 	CSX("CSX","C",AtomType.CH,AtomType.CN,AtomType.CCA,AtomType.CC,AtomType.CO,AtomType.CCB,true), //????
        SEP("SEP","S",AtomType.SH,AtomType.SN,AtomType.SCA,AtomType.SC,AtomType.SO,AtomType.SCB,true), //PHOSPHOSERINE
	ASQ("ASQ","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true), //PHOSPHOASPARTATE
	PHD("PHD","D",AtomType.DH,AtomType.DN,AtomType.DCA,AtomType.DC,AtomType.DO,AtomType.DCB,true), //PHOSPHORYLATION
	PCA("PCA","E",AtomType.EH,AtomType.EN,AtomType.ECA,AtomType.EC,AtomType.EO,AtomType.ECB,true), //PYROGLUTAMIC ACID
	ARO("ARO","R",AtomType.RH,AtomType.RN,AtomType.RCA,AtomType.RC,AtomType.RO,AtomType.RCB,true), //???
	HHS("HHS","H",AtomType.HH,AtomType.HN,AtomType.HCA,AtomType.HC,AtomType.HO,AtomType.HCB,true), //???
	DMY(" - ","-",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false), /* a dummy residue used 
													   as a place-holder 
													   or as a search key.*/
    CREATOR("XXX","X",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false),
    // for helium cluster minimization
    HEL("HEL","H",AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,AtomType.XXX,false);

    private final String nameThreeLetters;
    private final String nameOneLetter;
    private final AtomType caType,hType,nType,cType,oType,cbType;
    private final boolean standard;
    private ResidueType(String nameThreeLetters, String nameOneLetter, 
			AtomType hType, AtomType nType, AtomType caType, 
			AtomType cType, AtomType oType, AtomType cbType, boolean standard) {
	this.nameThreeLetters = nameThreeLetters;
	this.nameOneLetter = nameOneLetter;
	this.hType  = hType;
	this.nType  = nType;
	this.caType = caType;
	this.cType  = cType;
	this.oType  = oType;
	this.cbType = cbType;
	this.standard = standard;
    }
    public String nameThreeLetters() {return nameThreeLetters;}
    public String nameOneLetter() {return nameOneLetter;}
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
	    throw new RuntimeException("Weird Residue type "+name+" of length "+name.length()+".");
	for (ResidueType type:ResidueType.values()) {
	    if ((name.equals(type.nameThreeLetters)) | (name.equals(type.nameOneLetter))) return type;
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

}
			     

  
