package meshi.parameters;

public enum AtomType {
    //Termini 0-2
    TRN(BBatom.N, Element.N), // Charged Nitrogen at tne N-terminus 
    TRC(BBatom.C), // Carboxyl Carbon at the C-terminus
    TRO(BBatom.O), // Charged Oxygen at tne C-terminus 
    //Ala 3-8
    AH(BBatom.H, Element.H), AN(BBatom.N, Element.N), ACA(BBatom.CA), AC(BBatom.C), AO(BBatom.O, Element.O), ACB(BBatom.CB),
    //Cys 9-15
    CH(BBatom.H, Element.H), CN(BBatom.N, Element.N), CCA(BBatom.CA), CC(BBatom.C), CO(BBatom.O, Element.O),CCB(BBatom.CB),
    CSG(Element.S),
    //Asp 16-23
    DH(BBatom.H, Element.H), DN(BBatom.N, Element.N), DCA(BBatom.CA), DC(BBatom.C), DO(BBatom.O, Element.O),DCB(BBatom.CB), 
    DCG, 
    DOD(Element.O),
    //Glu 24-32
    EH(BBatom.H, Element.H), EN(BBatom.N, Element.N), ECA(BBatom.CA), EC(BBatom.C), EO(BBatom.O, Element.O),ECB(BBatom.CB),
    ECG,
    ECD,
    EOE(Element.O),
    //Phe 33-42
    FH(BBatom.H, Element.H), FN(BBatom.N, Element.N), FCA(BBatom.CA), FC(BBatom.C), FO(BBatom.O, Element.O),FCB(BBatom.CB),
    FCG,
    FCD,
    FCE,
    FCZ,
    //Gly 43-47
    GH(BBatom.H, Element.H), GN(BBatom.N, Element.N), GCA(BBatom.CA), GC(BBatom.C), GO(BBatom.O, Element.O),
    //His 48-60
    HH(BBatom.H, Element.H), HN(BBatom.N, Element.N), HCA(BBatom.CA), HC(BBatom.C), HO(BBatom.O, Element.O),HCB(BBatom.CB),
    HCG,
    HCD,
    HND(Element.N),
    HHD(Element.H),
    HCE,
    HNE(Element.N),
    HHE(Element.H),
    //Ile 61-69
    IH(BBatom.H, Element.H), IN(BBatom.N, Element.N), ICA(BBatom.CA), IC(BBatom.C), IO(BBatom.O, Element.O),ICB(BBatom.CB),
    ICG1,
    ICG2,
    ICD,
    //Lys 70-79
    KH(BBatom.H, Element.H), KN(BBatom.N, Element.N), KCA(BBatom.CA), KC(BBatom.C), KO(BBatom.O, Element.O),KCB(BBatom.CB),
    KCG,
    KCD,
    KCE,
    KNZ(Element.N),
    //Leu 80-88
    LH(BBatom.H, Element.H), LN(BBatom.N, Element.N), LCA(BBatom.CA), LC(BBatom.C), LO(BBatom.O, Element.O),LCB(BBatom.CB),
    LCG,
    LCD1,
    LCD2,
    //Met 89-97
    MH(BBatom.H, Element.H), MN(BBatom.N, Element.N), MCA(BBatom.CA), MC(BBatom.C), MO(BBatom.O, Element.O),MCB(BBatom.CB),
    MCG,
    MSD(Element.S),
    MCE,
    //Asn 98-108
    NH(BBatom.H, Element.H), NN(BBatom.N, Element.N), NCA(BBatom.CA), NC(BBatom.C), NO(BBatom.O, Element.O),NCB(BBatom.CB),
    NCG,
    NOD(Element.O),
    NND(Element.N),
    NHD1(Element.H),
    NHD2(Element.H),
    //Pro 109-115
                             PN(BBatom.N, Element.N), PCA(BBatom.CA), PC(BBatom.C), PO(BBatom.O, Element.O),PCB(BBatom.CB),
    PCG,
    PCD,

    //Gln 116-127
    QH(BBatom.H, Element.H), QN(BBatom.N, Element.N), QCA(BBatom.CA), QC(BBatom.C), QO(BBatom.O, Element.O),QCB(BBatom.CB),
    QCG,
    QCD,
    QOE(Element.O),
    QNE(Element.N),
    QHE1(Element.H),
    QHE2(Element.H),

    //Arg 128-139
    RH(BBatom.H, Element.H), RN(BBatom.N, Element.N), RCA(BBatom.CA), RC(BBatom.C), RO(BBatom.O, Element.O),RCB(BBatom.CB),
    RCG,
    RCD,
    RNE(Element.N),
    RHE(Element.H),
    RCZ,
    RNH(Element.N),
    //Ser
    SH(BBatom.H, Element.H), SN(BBatom.N, Element.N), SCA(BBatom.CA), SC(BBatom.C), SO(BBatom.O, Element.O),SCB(BBatom.CB),
    SOG(Element.O),
    //Thr
    TH(BBatom.H, Element.H), TN(BBatom.N, Element.N), TCA(BBatom.CA), TC(BBatom.C), TO(BBatom.O, Element.O),TCB(BBatom.CB),
    TCG,
    TOG(Element.O),
    //Val
    VH(BBatom.H, Element.H), VN(BBatom.N, Element.N), VCA(BBatom.CA), VC(BBatom.C), VO(BBatom.O, Element.O),VCB(BBatom.CB),
    VCG1,
    VCG2,
    //Trp
    WH(BBatom.H, Element.H), WN(BBatom.N, Element.N), WCA(BBatom.CA), WC(BBatom.C), WO(BBatom.O, Element.O),WCB(BBatom.CB),
    WCG,
    WCD1,
    WCD2,
    WCE2,
    WCE3,
    WNE(Element.N),
    WHE(Element.H),
    WCZ2,
    WCZ3,
    WCH2,
    //Tyr
    YH(BBatom.H, Element.H), YN(BBatom.N, Element.N), YCA(BBatom.CA), YC(BBatom.C), YO(BBatom.O, Element.O),YCB(BBatom.CB),
    YCG,
    YCD,
    YCE,
    YCZ,
    YOH(Element.O),
    // generic
    XXX,
    // for helium cluster minimization
    HE;

    private final BBatom bbAtom;
    private final Element element;
    private AtomType(BBatom bbAtom, Element element) {
	this.bbAtom = bbAtom;
	this.element = element;
    }
    private AtomType(Element element) {
	this.bbAtom = BBatom.NOT_BACKBONE;
	this.element = element;
    }
    private AtomType(BBatom bbAtom) {
	this.bbAtom = bbAtom;
	this.element = Element.C;
    }
    private AtomType() {
	this.bbAtom = BBatom.NOT_BACKBONE;
	this.element = Element.C;
    }

    public final boolean backbone() {return bbAtom != BBatom.NOT_BACKBONE;}
    public final boolean backboneH() {return bbAtom == BBatom.H;}
    public final boolean backboneN() {return bbAtom == BBatom.N;}
    public final boolean backboneCA() {return bbAtom == BBatom.CA;}
    public final boolean backboneC() {return bbAtom == BBatom.C;}
    public final boolean backboneO() {return bbAtom == BBatom.O;}
    public final boolean isCarbon() {return element == Element.C;}
    public final boolean isOxygen() {return element == Element.O;}
    public final boolean isNitrogen() {return element == Element.N;}
    public final boolean isSulfur() {return element == Element.S;}
    public final boolean isHydrogen() {return element == Element.H;}
    public final BBatom bbAtom() {return bbAtom;}
    public static AtomType type(String residueName,String atomName) {
	if (atomName.equals("OXT")) return TRO;
	ResidueType residueType = ResidueType.type(residueName);
	String name = residueType.nameOneLetter()+atomName;
	for(AtomType type:AtomType.values()) {
	    if (name.equals(type.toString())) return type;
	    else if ((name.length() >= 4) && 
		     name.substring(0,3).equals(type.toString())) return type;
	}
	if (residueType == ResidueType.ASN) {
	    if (atomName.equals("1HD2")) return AtomType.NHD1;
	    if (atomName.equals("HD21")) return AtomType.NHD1;
	    if (atomName.equals("2HD2")) return AtomType.NHD2;
	    if (atomName.equals("HD22")) return AtomType.NHD2;
	}
	if (residueType == ResidueType.GLN) {
            if (atomName.equals("1HE2")) return AtomType.QHE1;
            if (atomName.equals("HE21")) return AtomType.QHE1;
            if (atomName.equals("2HE2")) return AtomType.QHE2;
            if (atomName.equals("HE22")) return AtomType.QHE2;
	}
	return XXX;
    }
    public static AtomType type(String name) {
	for(AtomType type:AtomType.values()) {
	    if (name.equals(type.toString())) return type;
	}
	return XXX;
    }
     public static AtomType type(int ordinal) {
	for(AtomType type:AtomType.values())
	    if (ordinal == type.ordinal()) return type;
	return XXX;
    }
	 
	 
}
    //
