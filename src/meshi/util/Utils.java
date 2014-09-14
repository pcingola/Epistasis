package meshi.util;
import meshi.util.file.*;
import meshi.util.filters.*;
import meshi.util.dssp.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.PDB.*;
import meshi.energy.*;
import meshi.optimizers.*;
import meshi.parameters.*;
import meshi.sequences.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.ca.*;
import meshi.molecularElements.extendedAtoms.*;
import meshi.energy.simpleEnergyTerms.angle.*;
import meshi.energy.simpleEnergyTerms.plane.*;
import meshi.energy.simpleEnergyTerms.outOfPlane.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.*;
import meshi.applications.prediction.*;
import java.util.*;
import java.io.*;
/**
 * Where we store useful static methods that do not make sense anywhere else.
 * 
 **/

public class Utils implements MeshiPotential, KeyWords{
    public static final ExceptionHandler defaultExceptionHandler = new DefaultExceptionHandler();
    public static final ExceptionHandler doNothing = new DoNothingExceptionHandler();
 

    //----------------------------------- Commands -------------------------------------------------
    /**
     * A standard method to create the commands list.
     **/
    
    public static CommandList init(String[] args, String name) {
	return init(args, 1, ("Usage: java -XmxNNNm "+name+"  <commands file name>\n\n"+
			      "NNN is the size of the expected memory requirement in MegaBytes."));
    }
    
    public static CommandList init(String[] args, int numberOfParameters, String errorMessage) {
	if ((numberOfParameters > 0) && (args.length != numberOfParameters)) throw new RuntimeException(errorMessage);
	CommandList commands = new CommandList(args, 
					       new CommandsException(errorMessage));
	if (commands.keyExists(SEED)) 
	    MeshiProgram.initRandom(commands.firstWord(SEED).secondWordInt());
	return commands;
    }

    public static CommandList init(String[] args, int numberOfParameters, int seed, String errorMessage) {
	if ((numberOfParameters > 0) && (args.length != numberOfParameters)) throw new RuntimeException(errorMessage);
	CommandList commands = new CommandList(args,
					       new CommandsException(errorMessage));
    MeshiProgram.initRandom(seed);
	return commands;
    }

 
	
    //----------------------------------- Alignments ----------------------------------------------
    /**
     * Gets a residue alignment as a parameters and returns a new alignment that includes only residues with original atoms.
     **/
    public static ResidueAlignment getOriginalAtomsAlignment(ResidueAlignment alignment, OriginalAtoms originalAtoms, 
		    					     int row) {
	ResidueAlignment out = new ResidueAlignment();
	for (Iterator columns = alignment.iterator(); columns.hasNext();) {
	    ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
	    Residue residue = (Residue) column.cell(row).obj;
	    if (!residue.dummy()) {
		Atom atom = residue.ca();
		if (originalAtoms.accept(atom)) out.add(column);
	    }
	}
	return out;
    }
    
    public static void assignSecondaryStructure(Protein protein, CommandList commands, Key key) {
	assignSecondaryStructure(protein.chain(), commands, key);
    }

    public static void assignSecondaryStructure(Chain chain, CommandList commands, Key key) {
	String fileName = commands.firstWord(key).secondWord();
	SequenceList      ssList      = new SequenceList(fileName);
	SequenceAlignment ssAlignment = new SequenceAlignment(ssList);
	for (Iterator columns = ssAlignment.iterator(); columns.hasNext();) {
	    AlignmentColumn       column  = (AlignmentColumn) columns.next();
	    SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
	    SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
	    cell0.addAttribute(cell1);
	}
	ResidueSequence sequence = chain.sequence();
	SequenceAlignment sequenceAlignment = SequenceAlignment.identityAlignment(sequence, ssList.get(0));
	for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext();) {
	    SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentCell   from   = (SequenceAlignmentCell)column.cell(1);
	    SequenceAlignmentCell   from1  = (SequenceAlignmentCell) from.getAttribute(MeshiAttribute.SEQUENCE_ALIGNMENT_COLUMN_ATTRIBUTE);
	    char                    ss     = from1.getChar();
	    SequenceAlignmentCell   to     = (SequenceAlignmentCell)column.cell(0);
	    Residue                residue = (Residue) to.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    if (residue != null) {
		if (ss == 'H') residue.setSecondaryStructure(SecondaryStructure.HELIX);
		else if (ss == 'E') residue.setSecondaryStructure(SecondaryStructure.SHEET);
		else if (ss == 'C') residue.setSecondaryStructure(SecondaryStructure.COIL);
		else if (ss == 'A') residue.setSecondaryStructure(SecondaryStructure.ALL);
		else throw new RuntimeException("Unrecognized secondary structure "+ss);
	    }	    
	}
    }

    public static void assignSecondaryStructure(Protein protein, SequenceAlignment targetSecondaryStructure) {
	Iterator columns = targetSecondaryStructure.iterator();
	SequenceAlignmentColumn column;

	// Assign secondary structure.
	boolean first = true; //for the first dummy residue
	for (Residue residue:protein.chain()){
	    if (first) first = false;
	    else {
		try {
		    column = (SequenceAlignmentColumn) columns.next();
		} catch (Exception ex) {throw new RuntimeException("Alignment: \n"+targetSecondaryStructure+"\n"+
							       "Does not fit protein: "+protein+"\n"+ex);}
		String resName = column.cell(0).obj.toString();
		if (!resName.equals(residue.type().nameOneLetter()))
		    {throw new RuntimeException("\n"+
						"Alignment: \n"+targetSecondaryStructure+"\n"+
						"Does not fit protein:\n"+protein+"\n"+
						"resName = "+resName+"\n"+
						"residue.type() = "+residue.type()+"\n"+
					    "residue.type().nameOneLetter() = "+residue.type().nameOneLetter());}
		char ss = ((SequenceAlignmentCell)column.cell(1)).getChar();
		if (ss == 'H') residue.setSecondaryStructure(SecondaryStructure.HELIX);
		else if (ss == 'E') residue.setSecondaryStructure(SecondaryStructure.SHEET);
		else if (ss == 'C') residue.setSecondaryStructure(SecondaryStructure.COIL);
		else if (ss == 'A') residue.setSecondaryStructure(SecondaryStructure.ALL);
		else throw new RuntimeException("Unrecognized secondary structure "+ss);
	    }
	}
    }

    public static ResidueAlignment alignProteinsByAlignmentFile(CommandList commands, Key key,
								Chain chain0, Chain chain1) {
        Command      alignmentCommand        = commands.firstWord(key);
        SequenceList sequenceList            = new SequenceList(alignmentCommand.secondWord());
	return new ResidueAlignment(chain0, chain1, sequenceList);
    }

    //----------------------------------- DistanceMatrix -------------------------------------------
    /**
     * A standard method to build distance matrices.
     **/
    public static DistanceMatrix getDistanceMatrix(AtomList atomList, CommandList commands) {
	double rmax = commands.firstWord(R_MAX).secondWordDouble();
	double buffer = commands.firstWord(BUFFER_SIZE).secondWordDouble();
	double edge = commands.firstWord(GRID_EDGE).secondWordDouble();
	return new DistanceMatrix(atomList.molecularSystem(), rmax, buffer, edge, DistanceMatrix.DEFAULT_BONDED_LIST_DEPTH);
    }

    //----------------------------------- energy -------------------------------------------------
    public static double[] getEnergyValues(TotalEnergy totalEnergy) {
	double[] out = new double[totalEnergy.energyTerms().size()];
	boolean[] on = new boolean[totalEnergy.energyTerms().size()];
	int i = 0;
	for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
	    AbstractEnergy term = (AbstractEnergy) terms.next();
	    on[i] = term.isOn();
	    i++;
	}
	i = 0;
	for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
	    AbstractEnergy term = (AbstractEnergy) terms.next();
	    totalEnergy.off();
	    term.on();
	    out[i] = totalEnergy.evaluate();
	    i++;
	}
	i = 0;
	totalEnergy.off();
	for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
	    AbstractEnergy term = (AbstractEnergy) terms.next();
	    if (on[i]) term.on();
	    i++;
	}
	return out;
    }
    public static String[] getEnergyTermsNames(TotalEnergy totalEnergy) {
	String[] out = new String[totalEnergy.energyTerms().size()];
	int i = 0;
	for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
	    AbstractEnergy term = (AbstractEnergy) terms.next();
	    out[i] = term.comment();
	    i++;
	}
	return out;
    }

    public static TotalEnergy relax(AtomList atoms, Protein protein, EnergyCreator[] energyCreators, CommandList commands) {
	return relax(atoms, protein, energyCreators, commands, RELAX);
    }
    public static TotalEnergy oldRelax(AtomList atoms, Protein protein, EnergyCreator[] energyCreators, CommandList commands) {
	return oldRelax(atoms, protein, energyCreators, commands, RELAX);
    }

    public static TotalEnergy oldRelax(AtomList atoms, Protein protein, EnergyCreator[] energyCreators, CommandList commands, Key key) {
	OldDistanceMatrix distanceMatrix = new OldDistanceMatrix(atoms, 5.5, 2.0, DistanceMatrix.DEFAULT_EDGE(5.5,2.0), 4); 

	TotalEnergy energy = new TotalEnergy(protein, distanceMatrix, energyCreators, commands);
	Optimizer optimizer =  Utils.getLBFGS(energy, commands,key);
	try {
	    optimizer.run();
	}
	catch (OptimizerException ex) {
	    energy.test();
	    throw new RuntimeException(ex);
	}
	return energy;
    }
    public static TotalEnergy relax(AtomList atoms, Protein protein, EnergyCreator[] energyCreators, CommandList commands, Key key) {
	DistanceMatrix distanceMatrix = new DistanceMatrix(atoms.molecularSystem(), 5.5, 2.0, DistanceMatrix.DEFAULT_EDGE(5.5,2.0), 4); 

	TotalEnergy energy = new TotalEnergy(protein, distanceMatrix, energyCreators, commands);
	Optimizer optimizer =  Utils.getLBFGS(energy, commands,key);
	try {
	    optimizer.run();
	}
	catch (OptimizerException ex) {
	    energy.test();
	    throw new RuntimeException(ex);
	}
	return energy;
    }
  
    //----------------------------------- io and files ----------------------------------------------------
    public static StringList getStructureNames(CommandList commands){
	String fileName = commands.firstWord(STRUCTURE_NAMES).secondWord();
	return new StringList(new MeshiLineReader(fileName));
    }

    public static MeshiWriter  getOutputFile(CommandList commands){
	String fileName = commands.firstWord(OUTPUT_FILE_NAME).secondWord();
	try {
	    return new MeshiWriter(fileName);
	} catch (Exception ex) {throw new RuntimeException(ex);}
    }
    /**
     * Opens a file for writing (warped by MeshiWriter object) with a name and location specified
     * by the commands and a .<seed>.<runNumber>.pdb suffix.
     **/
    public static MeshiWriter newPdbWriter(CommandList commands, Key pathKey, Key nameKey, int runNumber) {
	try {
	    String suffix = "."+MeshiProgram.seed()+"."+runNumber+".pdb";
	    return new MeshiWriter(commands, pathKey, nameKey, suffix);
	}
	catch(Exception ex) {
	    throw new RuntimeException("A problem in opening output file with key words "+
				       pathKey+" and "+nameKey+"\n"+ex);
	}
    }

    //----------------------------------- analysis ------------------------------------------------
   public static void  RmsGdtEnergy(Object output, ResidueSequence  model, ResidueSequence  reference, TotalEnergy totalEnergy,String tag) {
	ResidueAlignment modelRefeenceAlignment = null;
	int              referenceSize          = -1;
	if (reference != null) {
	    modelRefeenceAlignment = new ResidueAlignment(model, reference); 
	    referenceSize          = reference.size();
	}
	if (output instanceof StringList )  RmsGdtEnergy((StringList ) output,modelRefeenceAlignment, totalEnergy, tag, referenceSize);
	else 
	    if (output instanceof MeshiWriter)  RmsGdtEnergy((MeshiWriter) output,modelRefeenceAlignment, totalEnergy, tag, referenceSize);
	    else throw new RuntimeException("Wrong first parameter for RmsGdtEnergy "+output);
    }

    public static void  RmsGdtEnergy(Object output, Protein model, Protein reference, TotalEnergy totalEnergy,String tag) {
	ResidueAlignment modelRefeenceAlignment = null;
	int              referenceSize          = -1;
	if (reference != null) {
	    modelRefeenceAlignment = new ResidueAlignment(reference.chain(), reference.name(),model.chain(), model.name());
	    referenceSize          = reference.atoms().CAFilter().size();
	}
	if (output instanceof StringList )  RmsGdtEnergy((StringList ) output,modelRefeenceAlignment, totalEnergy, tag, referenceSize);
	else 
	    if (output instanceof MeshiWriter)  RmsGdtEnergy((MeshiWriter) output,modelRefeenceAlignment, totalEnergy, tag, referenceSize);
	    else throw new RuntimeException("Wrong first parameter for RmsGdtEnergy "+output);
    }
	

    public static void  RmsGdtEnergy(MeshiWriter output, ResidueAlignment modelRefeenceAlignment, TotalEnergy totalEnergy,String tag, int refLength) {
	StringList stringList = new StringList();
	RmsGdtEnergy(stringList, modelRefeenceAlignment, totalEnergy, tag, refLength);
	stringList.print(output);
	output.flush();
    }

    public static void  RmsGdtEnergy(StringList output, ResidueAlignment modelRefeenceAlignment, TotalEnergy totalEnergy,String tag, int refLength) {
	double rms        = -1;
	double rmsOrig    = -1;
	double[] gdt      = {-1.0, -1.0, -1.0, -1.0, -1.0};
	double[] gdtOrig  = {-1.0, -1.0, -1.0, -1.0, -1.0};
	String header = "H_"+tag;
	String values = "V_"+tag;
	if (modelRefeenceAlignment != null) {
	    ResidueAlignment origRefeenceAlignment  = (ResidueAlignment) modelRefeenceAlignment.filter(OriginalAtom.filter);
	    rms        = Rms.rms(modelRefeenceAlignment);
	    if (origRefeenceAlignment.size() > 3)
		rmsOrig    = Rms.rms(origRefeenceAlignment);
	    gdt        = Rms.gdt(modelRefeenceAlignment, refLength); 
	    if (origRefeenceAlignment.size() > 3)
		gdtOrig    = Rms.gdt(origRefeenceAlignment, refLength); 
	    header += String.format("%6s %6s %5s %5s  %5s  %5s  %5s  %5s ","RMS_O","RMS","GDT_O","GDT","GDT1","GDT2","GDT3","GDT4");
	    values += String.format(" %6.3f %6.3f %5.3f %5.3f  %5.3f  %5.3f  %5.3f  %5.3f ",
				      rmsOrig, rms, gdtOrig[0], gdt[0], gdt[1], gdt[2], gdt[3] , gdt[4]);
	}
	if (totalEnergy != null) {
	    double energy = totalEnergy.evaluate();
	    header += String.format(" %12s","energy");
	    values += String.format(" %12.3f",energy);

	    boolean[] on = new boolean[totalEnergy.energyTerms().size()];
	    int i = 0;
	    for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
		AbstractEnergy term = (AbstractEnergy) terms.next();
		on[i] = term.isOn();
		i++;
	    }

	    i = 0;
	    for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
		AbstractEnergy term = (AbstractEnergy) terms.next();
		if (on[i]){
		    int length = term.comment().length();
		    if (length > 12) length = 12;
		    header += String.format(" %12s",term.comment().substring(0,length));
		    totalEnergy.off();
		    term.on();
		    values += String.format(" %12.3f",totalEnergy.evaluate());
		}
		i++;
	    }
	    
	    i = 0;
	    for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
		AbstractEnergy term = (AbstractEnergy) terms.next();
		if (on[i]) term.on();
		i++;
	    }
	}
	output.add(header);
	output.add(values);
    }
	

    //----------------------------------- Misc. ----------------------------------------------------
    /**
     * Static version of equals. Useful when the compared objects may be null.
     * This is NOT the right place for highly sofisticate algorithms. 
     **/
    public static boolean equals(Object obj1, Object obj2) {
	if (obj1 == null)
	    return (obj2 == null);
	else
	    return (obj1.equals(obj2));
    }
    
    /** 
     * Converts a string to an int.
     **/
    public static int toInt(String s) {
	return Integer.valueOf(s.trim()).intValue();
    }
    /** 
     * Converts a string to a double.
     **/
    public static double toDouble(String s) {
	return Double.valueOf(s.trim()).doubleValue();
    }

    // Used for loop building
    public static int[] getSeqOfProt(Protein prot, int fromRes, int toRes) {
        int[] output = new int[toRes-fromRes+1];
        for (int c=fromRes ; c<=toRes ; c++)
                if (prot.residue(c) == null)
                        throw new RuntimeException("\nCurrent setting asked for type of residue " + c + "which doesn't exist in the protein\n");
                else
		    output[c-fromRes] = prot.residue(c).type().ordinal();
        return output;
    }

   //----------------------------------- MolecularSystem ----------------------------------------------------
    public static double radius(MolecularSystem molecularSystem) {
	double sum2;
  	double dx, dy, dz, d2;
	double cmx, cmy, cmz; // center of mass x, y and z
	sum2 = 0;
	cmx = cmy = cmz = 0.0;
	int nWithCoordinates = 0;
	for (AtomCore atom:molecularSystem) {
	    if (! atom.status().nowhere()) {
	    cmx += atom.x();
	    cmy += atom.y();
	    cmz += atom.z();
	    nWithCoordinates++;
	  }
	}
	cmx /= nWithCoordinates;
	cmy /= nWithCoordinates;
	cmz /= nWithCoordinates;
	for (AtomCore atom:molecularSystem) {
	    if (! atom.status().nowhere()) {
		dx = cmx - atom.x();
		dy = cmy - atom.y();
		dz = cmz - atom.z();
		d2 = dx*dx + dy*dy + dz*dz;
		sum2 += d2;
	    }
	}
	return Math.sqrt(sum2 / nWithCoordinates);
    }
    
    public static void testMolecularSystemIntegrity(MolecularSystem ms, String comment) {
	int size = ms.size();
	System.out.println(comment+"   Molecular System of size = "+size);
	int normal = 0, frozen = 0, nowhere = 0, hidden = 0, clashes = 0;
	for (int i = 0; i < size; i++) {
	    AtomCore ai = ms.get(i);
	    if (ai.status() == AtomStatus.NORMAL) normal++;
	    else if (ai.status() == AtomStatus.FROZEN) frozen++;
	    else if (ai.status() == AtomStatus.NOWHERE) nowhere++;
	    else if (ai.status() == AtomStatus.HIDDEN) hidden++;
	    else throw new RuntimeException("MolecularSystemError:\n"+"weird AtomCore instance "+ai);
	    for (int j = i+1; j < size; j++) {
		AtomCore aj = ms.get(j);
		if (ai.atom.name.equals(aj.atom.name) && 
		    (ai.atom.residue().number() == aj.atom.residue().number()) && 
		    (ai.atom.residue().chain() == aj.atom.residue().chain()))
		    throw new RuntimeException("MolecularSystemError:\n"+"Two apparently identical atoms\n"+ai.atom+"\n"+aj.atom);
		if ((!ai.status().nowhere()) && (ai.x() == aj.x()) &&  (ai.y() == aj.y()) &&  (ai.z() == aj.z()))
		    throw new RuntimeException("MolecularSystemError:\n"+"Two clashing atoms\n"+ai.atom+"\n"+aj.atom);
	    }
	}
	System.out.println("normal  = "+normal+"\n"+
			   "frozen  = "+frozen+"\n"+
			   "nowhere = "+nowhere+"\n"+
			   "hidden  = "+hidden+"\n");
    }
	
    public static void hideAllBut(AtomList atomList) {
	MolecularSystem ms = atomList.molecularSystem();
	for (AtomCore atomCore:ms) {
	    if ((!atomCore.atom.nowhere()) & (!atomList.contains(atomCore.atom))) atomCore.atom.hide();
	}
    }
    //----------------------------------- Optimization ----------------------------------------------------
    public static LBFGS getLBFGS(TotalEnergy energy, CommandList commands) {
 	int maxSteps = commands.firstWord(MAX_STEPS).secondWordInt();
	double tolerance = commands.firstWord(TOLERANCE).secondWordDouble();
	int reportEvery = commands.firstWord(REPORT_EVERY).secondWordInt();
	return new LBFGS(energy,tolerance,maxSteps,reportEvery);
    }

    public static LBFGS getLBFGS(TotalEnergy energy, CommandList commands, Key key) {
	CommandList commands1 = commands.firstWordFilter(key);
	int maxSteps = commands1.secondWord(MAX_STEPS).thirdWordInt();
	double tolerance = commands1.secondWord(TOLERANCE).thirdWordDouble();
	int reportEvery = commands1.secondWord(REPORT_EVERY).thirdWordInt();
	return new LBFGS(energy,tolerance,maxSteps,reportEvery);
    }
   



    //----------------------------------- Proteins ----------------------------------------------------
    public static boolean addAtoms(Protein model, CommandList commands) { 
	return addAtoms(model,commands,100);
    }
    public static boolean addAtoms(Protein model, CommandList commands, int maxNumberOfRandomCoordinatesPerResidue) {
	EnergyCreator[] energyCreatorsBondedTermsOnly = {  
	    new BondCreator(),
	    new AngleCreator(),
	    new PlaneCreator(),
	    new OutOfPlaneCreator(),
	    new RamachandranSidechainEnergyCreator()
	};
	boolean  nonDummyFound = false;
	
	for (Iterator residues = model.residues().iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    if (residue.dummy() & nonDummyFound) {
		System.out.println("Cannot evaluate model "+model.name()+" "+residue+" is dummy");
		return false; //there is a missing residue
	    }
	    if (!residue.dummy()) nonDummyFound = true;
	}
	
	model.atoms().freeze();
	for (Iterator atoms = model.atoms().iterator(); atoms.hasNext();) {  
	    Atom atom         = (Atom) atoms.next();	  
            Residue residue   = atom.residue();
	    MeshiAttribute la = residue.getAttribute(MeshiAttribute.LOOP_RESIDUE);
	    /*if ((! atom.nowhere()) && (la == null))
                 atom.addAttribute(MeshiAttribute.ORIGINAL_ATOM);
        */
	} 
	
	for (Iterator residues = model.residues().iterator();residues.hasNext();){ 
	    Residue residue = (Residue) residues.next(); 
	    if (!residue.dummy()) {  
		boolean OK = Utils.assignRandomCoordinates(residue, commands, maxNumberOfRandomCoordinatesPerResidue); 
                if (! OK ) return false;
            }
	}
	Utils.relax(model.atoms(),model, energyCreatorsBondedTermsOnly, commands); 
	model.atoms().defrost();
        return true;
    }
    
    public static void printCaspFormat(ModelData modelData, CommandList commands,int modelNum) {
	String fileName = (String) modelData.get(0).value();
	Protein model   = new Protein(fileName, new PdbLineATOM(),ResidueExtendedAtomsCreator.creator);
	String target   = commands.firstWord(TARGET_NAME).secondWord();  /* with four digits T0302 T0277 etc. */
	String author   = commands.firstWord(CASP_GROUP).secondWord();  /* xxxxx-xxxxx-xxxxx */
	String parent   = commands.firstWord(TEMPLATE_NAME).secondWord();  /* 1qwe_A or 1qwe  */
	String method   = commands.firstWord(METHOD).secondWord();  /* Some description */
	try {
	    MeshiWriter mw = new MeshiWriter(target+"."+modelNum+".pdb");
	    mw.println("PFRMAT TS ");
	    mw.println("TARGET " + target);
	    mw.println("AUTHOR " + author);
	    mw.println("METHOD " + method);
	    mw.println("MODEL " + modelNum);
	    mw.println("PARENT " + parent);
	    for (Pair pair:modelData) mw.println("REMARK "+pair);
	    AtomList atoms = model.atoms();
	    for (int c=0 ; c<atoms.size() ; c++) 
    		if (!atoms.atomAt(c).nowhere()) mw.println(atoms.atomAt(c) + "                ");
	    mw.println("TER ");
	    mw.println("END ");
	    mw.close();
	}
	catch (Exception ex) {throw new RuntimeException(ex);}
    }

    public static int numberOfNonDummyResidues(Protein protein) {
	int number = 0;
	for (Iterator residues = protein.residues().iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    if (! residue.dummy()) number++;
	}
	return  number;
    }
    public static int numberOfAtomsWithCoordinates(Protein protein) {
	int number = 0;
	for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    if (! atom.nowhere()) number++;
	}
	return  number;
    }
    
    /** 
     * Reads a protein structure from the file pointed at by the command list.
     **/
    public static String getProteinNameFromPdbFileName(String fullFileName) {
	String fileName = (new File(fullFileName)).getName();	StringTokenizer words  = new StringTokenizer(fileName,".");
	return words.nextToken();
    }

    public static ResidueList getFragment(Protein protein, int start, int end) {
	Chain chain = protein.chain();
	ResidueList out = new ResidueList();
	for (int i = start; i <= end; i++)
	    out.add( chain.get(i));
	return out;
    }

    public static void addHydrogens(Protein protein, CommandList commands) {
	BondParametersList  bondParametersList = Utils.getBondParameters(commands);
	AngleParametersList angleParametersList = Utils.getAngleParameters(commands);
	for (Iterator iter = protein.residues().iterator(); iter.hasNext();) {
	    Residue residue = (Residue) iter.next();
	    if ((! residue.dummy()) & (residue instanceof ResidueExtendedAtoms))
		((ResidueExtendedAtoms) residue).addHydrogens(bondParametersList,angleParametersList);
	}	
    }

    /**
     * Creates a Protein object in a new MolecularSystem. 
     **/
    public static Protein getProtein(CommandList commands, 
				     Key key, ResidueCreator creator, ExceptionHandler exceptionHandler) {
	MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
	new MolecularSystem();
	
	try {
	    Protein out = new Protein(commands.firstWord(key).secondWord(),
				      new PdbLineATOM(),
				      creator); 
	    addHydrogens(out,commands);

	    MolecularSystem.setCurrentMolecularSystem(saveMS);
	    
	    return out;
	} catch (Exception ex) {exceptionHandler.handle(ex);}

	MolecularSystem.setCurrentMolecularSystem(saveMS);
	return null;
    }

    public static Protein getProtein(CommandList commands, String fileName, ResidueCreator creator, ExceptionHandler exceptionHandler) {
	MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
	new MolecularSystem();
	
	try {
	    Protein out = new Protein(fileName,
				      new PdbLineATOM(),
				      creator); 
	    addHydrogens(out,commands);

	    MolecularSystem.setCurrentMolecularSystem(saveMS);
	    
	    return out;
	} catch (Exception ex) {exceptionHandler.handle(ex);}

	MolecularSystem.setCurrentMolecularSystem(saveMS);
	return null;
    }

    /** 
     * Reads a protein structure from the file pointed at by the command list. The protein is created in a new MolecularSystem.
     **/
    public static Protein getProtein(CommandList commands, 
				     Key key1, Key key2, ResidueCreator creator, ExceptionHandler exceptionHandler) {
	MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
	new MolecularSystem();
	try {
	    CommandList commands1 = commands.firstWordFilter(key1);
	    Protein out = new Protein(commands1.secondWord(key2).thirdWord(),
			       new PdbLineATOM(),
			       creator); 
	    addHydrogens(out,commands);
	    MolecularSystem.setCurrentMolecularSystem(saveMS);
	    return out;
	} catch (Exception ex) {exceptionHandler.handle(ex);}
	MolecularSystem.setCurrentMolecularSystem(saveMS);
	return null;
    }

    public static Protein getReferenceProtein(CommandList commands) {
	if (commands.keyExists(SUPERIMPOSE)) 
	    return getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator,Utils.defaultExceptionHandler);
	else return null;
    }
	
    public static void assignRandomCaCoordinates(Chain chain) {
	// Assign arbitrary (extended) coordinate
	Coordinates prev = new Coordinates(0.0, 0.0, 0.0);
	boolean first = true;//the dummy residue
	for(Residue residue:chain) {
	    if (first) first = false;
	    else {
		if (! residue.dummy()) {
		    residue.ca().setXYZ(prev);
		    double dx = MeshiProgram.randomNumberGenerator().nextDouble()*3.8;
		    double tmp = Math.sqrt(3.8*3.8-dx*dx);
		    double dy = MeshiProgram.randomNumberGenerator().nextDouble()*tmp;
		    double dz = Math.sqrt(3.8*3.8-dx*dx-dy*dy);
		    residue.ca().addToX(dx);
		    residue.ca().addToY(dy);
		    residue.ca().addToZ(dz);
		    prev = new Coordinates(residue.ca());
		}
        }
	}
    }

    public static boolean assignRandomCoordinates(Residue residue, CommandList commands, int maxMissingAtoms) {
	List<Atom> nowhereAtoms      = new ArrayList<Atom>();
	List<Atom> nowhereHeavyAtoms = new ArrayList<Atom>();
	for (Iterator atoms = residue.atoms().iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
  	    if (atom.nowhere()) nowhereAtoms.add(atom);
	    if (atom.nowhere() & (!atom.type().isHydrogen())) nowhereHeavyAtoms.add(atom);
	}
	if (nowhereAtoms.size() == 0) {
	    System.out.println("assignRandomCaCoordinates to "+residue+" nothing to do");
	    return true;
	}
	if (nowhereHeavyAtoms.size() > maxMissingAtoms) {
	    System.out.println("assignRandomCaCoordinates to "+residue+" missig atoms "+nowhereHeavyAtoms.size()+" maximum "+ maxMissingAtoms+" do nothing");
	    return false;
	}
	for (Atom atom : nowhereAtoms) {
	    Atom neighbor = getCovalentNeighbor(atom);
	    if (neighbor == null) {
		System.out.println("\n\n Problem in Util.assignRandomCoordinates(Residue residue, CommandList commands, int maxMissingAtoms)");
		System.out.println("no CovalentNeighbor to "+atom);
		residue.atoms().print();
		throw new RuntimeException("This is weird"); 
	    }
	    atom.setXYZ(new Coordinates(new Coordinates(neighbor),3));
	    System.out.println("assignRandomCaCoordinates to "+atom);
	}
	// DistanceMatrix distanceMatrix = new DistanceMatrix(residue.atoms(), 5.5, 2.0, 4);
// 	EnergyCreator[] energyCreators = {  
// 	    new BondCreator(),
// 	    new AngleCreator(),
// 	    new PlaneCreator(),
// 	    new OutOfPlaneCreator(),
// 	    new ExcludedVolCreatorOLD(),
// 	    new RamachandranSidechainEnergyCreator()
// 	};
// 	TotalEnergy energy = new TotalEnergy(new Protein(residue), 
// 					     distanceMatrix, energyCreators, commands);
// 	Optimizer optimizer =  getLBFGS(energy, commands,RELAX);
// 	try {
// 	    optimizer.run();
// 	}
// 	catch (OptimizerException ex) {
// 	    energy.test();
// 	    throw new RuntimeException(ex);
// 	}
	return true;
    }

    public static Atom getCovalentNeighbor(Atom atom) {
	for (Iterator bonded = atom.bonded().iterator(); bonded.hasNext();) {
	    Atom bondedAtom = (Atom) bonded.next();
	    if (!bondedAtom.nowhere()) return bondedAtom;
	    else if (bondedAtom.number() < atom.number()) 
		return getCovalentNeighbor(bondedAtom);
	}
	return getCovalentNeighborUp(atom);
    }
    /**
     * May be used only by getCovalentNeighbor(Atom atom)
     **/
    private static Atom getCovalentNeighborUp(Atom atom) {
	for (Iterator bonded = atom.bonded().iterator(); bonded.hasNext();) {
	    Atom bondedAtom = (Atom) bonded.next();
	    if (!bondedAtom.nowhere()) return bondedAtom;
	    else if (bondedAtom.number() > atom.number()) 
		return getCovalentNeighborUp(bondedAtom);
	}
	return null;
    }
    

    public static void colorByEnergy(AtomList atomList) {
	double min = 10000; 
	double max = -10000;
	for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    double energy = atom.energy();
	    if (energy < min) min = energy;
	    if (energy > max) max = energy;
	}
	double range = max - min;
	for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
	    Atom atom = (Atom) atoms.next();
	    atom.setTemperatureFactor(99.9*((atom.energy()-min)/range));
	}
    }


    public static void assignBackboneCoordinates(AtomList modelAtoms, AtomList tempAtoms){
	Atom[] modelArray = modelAtoms.toArrayOfAtoms();
	Atom[] tempArray = tempAtoms.toArrayOfAtoms();
	Comparator comparator = new AtomComparator();
	Arrays.sort(modelArray,comparator);
	Arrays.sort(tempArray,comparator);
	int i=0;
	Atom oxt = null;
	for (Atom modelAtom:modelArray) {
	    if (!modelAtom.name().equals("OXT")) {
		    while((i < tempArray.length) &&
			  ((tempArray[i].residue().number() != modelAtom.residue().number()) |
			   (! tempArray[i].name().equals(modelAtom.name())))) {
			i++;
		    }
		    if (i >= tempArray.length) throw new RuntimeException("Can't find coordinates for "+modelAtom);
		    if ((!modelAtom.type().backboneCA()) || 
			(modelAtom.distanceFrom(tempArray[i]) < 0.8))
			modelAtom.setXYZ(new Coordinates(tempArray[i]));
		    i++;
	    }
	    else oxt = modelAtom;
	}
	Atom ca = oxt.residue().ca();
	double x = ca.x()+1;
	oxt.setXYZ(new Coordinates(x, ca.y(), ca.z())); 
    }

    public static void AssignDSSP(Protein protein, CommandList commands, Key key) {
        String dsspFileName = commands.firstWord(key).secondWord();
        AssignDSSP(protein, dsspFileName);    
    }

 public static void AssignDSSP(Protein protein, String dsspFileName) {

    DSSP dssp = new DSSP(dsspFileName);
	SequenceList OneDimView = dssp.sequenceList();
	Sequence residueSequence = OneDimView.getResidueSequence();
	SequenceAlignment OneDimAlignment = new SequenceAlignment(OneDimView);

	// Hang each column on its first cell
	Iterator columns = OneDimAlignment.iterator();
	Iterator resIter = residueSequence.iterator();
	while (columns.hasNext()) {
	    SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentColumn res = (SequenceAlignmentColumn) resIter.next();
	    ((SequenceAlignmentCell) res.cell(0)).addAttribute(column);
	}
	//
	  
	Sequence fullSequence = protein.sequence();
	SequenceAlignment alignment = SequenceAlignment.identityAlignment(fullSequence, residueSequence);
	columns = alignment.iterator();
	while (columns.hasNext()) {
	    SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
	    SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
	    Residue res = protein.chain().get(cell0.number);
//        Residue res = protein.residues().get(cell0.number);
        if (!res.dummy()) {
		SequenceAlignmentColumn OneDcolumn = (SequenceAlignmentColumn) cell1.getAttribute();
		if (OneDcolumn == null) {
		    res.setSecondaryStructure('C');
		    res.setAccessibility('A');
		}
		else {
		    SequenceAlignmentCell OneDcell1 = (SequenceAlignmentCell) OneDcolumn.cell(1);
		    SequenceAlignmentCell OneDcell2 = (SequenceAlignmentCell) OneDcolumn.cell(2);
		    char ss = OneDcell1.getChar();
		    char acc = OneDcell2.getChar();
		    res.setSecondaryStructure(ss);
		    res.setAccessibility(ss);
		}
	    }
	}
    }

    public static AtomList getNeighbors(ResidueList residueList, Protein protein, double clashDistance){
	List<Atom> all = new ArrayList<Atom>();
	for (Iterator residues = residueList.iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    for (Iterator atoms = residue.atoms().iterator(); atoms.hasNext();) {
		Atom atom = (Atom) atoms.next();
		if ((!atom.nowhere()) && (!all.contains(atom))) all.add(atom);
		List<Atom> tempList  = getNeighbors(atom, protein, clashDistance);
		for (Atom candidate:tempList) 
		    if ((!candidate.nowhere()) &&(!all.contains(candidate))) all.add(candidate);
	    }
	}
	Object[] tempArray = all.toArray();
	Arrays.sort(tempArray);

	AtomList out = new AtomList();
	for (Object atom:tempArray) {
	    if (weirdAtom((Atom)atom)) throw new RuntimeException("This is weird "+atom+" "+((Atom)atom).nowhere());
	    out.add((Atom)atom);
	}
	return out;
    }

    public static boolean weirdAtom(Atom atom) {
	if (atom.nowhere()) return true; 
	if (atom.x() < -1000) return true;
	if (atom.y() < -1000) return true;
	if (atom.z() < -1000) return true;
	if (atom.x() >  1000) return true;
	if (atom.y() >  1000) return true;
	if (atom.z() >  1000) return true;
	return false;
    }

    public static List<Atom> getNeighbors(Atom atom, Protein protein, double clashDistance){
	ArrayList<Atom> out = new ArrayList<Atom>();
	for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
	    Atom other = (Atom) atoms.next();
	    if ((!atom.nowhere())  && (atom != other) && (!other.nowhere()))
		if (atom.distanceFrom(other) < clashDistance*3) 
		    out.add(other);
	}
	return out;
    }

    public static int getClashes(Atom atom, List<Atom>neighbors, double clashDistance) {
	int out = 0;
	for (Atom neighbor: neighbors) {
	    if (atom.distanceFrom(neighbor) < clashDistance) out++;
	    if (atom.distanceFrom(neighbor) < 2) out+=20;
	}
	return out;
    }


    //----------------------------------- Sequences ----------------------------------------------------
    
    public static ResidueSequence getResidueSequence(CommandList commands, Key key) {
	String fileName=commands.firstWord(key).secondWord();
	FastaList fastaList = new FastaList(fileName);
	return new ResidueSequence(fastaList);
    }


    //----------------------------------- Lists ----------------------------------------------------
    public static void print(Iterable list) {
	for (Object obj:list) 
	    System.out.println(obj);
    }

    /**
     * Print the list elements in raws.
     * The number of elements in a raw is determined by the parameter
     **/
    public static void print(ArrayList list, int rawLength, String format) {
        if (rawLength < 1) throw new RuntimeException("Raw Length cannot be shorter than one.");
	int i = 0;
	for (Object obj:list) {
	    System.out.printf(format,obj);
	    if (i%rawLength == 0) System.out.println();
	    i++;
        }
        System.out.println();
    }

    public static ArrayList filter(ArrayList list, Filter filter,ArrayList out) {
	if (out.size() != 0) throw new RuntimeException("out needs to be empty");
	for (Object obj:list)
	    if (filter.accept(obj)) out.add(obj);
	return out;
    }
    public static DistanceList filter(DistanceList list, Filter filter,DistanceList out) {
	if (out.size() != 0) throw new RuntimeException("out needs to be empty");
	for (Distance dis:list)
	    if (filter.accept(dis)) out.add(dis);
	return out;
    }

  //-------------------------------------------------------------------------------------------------------------------------------------------
    public static void testNonBondedList(DistanceList nonBondedList){
      FreeDistance freeDistance;
         for (Distance distance :nonBondedList){
             if (distance.mode() == DistanceMode.INFINITE)
                                 throw new RuntimeException("Infinite distance in NonBondedList"+distance);
             freeDistance = new FreeDistance(distance.atom1(), distance.atom2());
             if (freeDistance.distance() != distance.distance())
                                  throw new RuntimeException("Unupdated distance in NonBondedList:"+distance+"\nReal distance is "+freeDistance);
         }
}
    //------------------------------------------------------------------------------------------------------------------------------------------   
    
     public static AtomList duplicateInAnewMolecularSystem(AtomList atomList) {
       MolecularSystem saveMS = MolecularSystem.currentMolecularSystem();
       new MolecularSystem();

       AtomList out = new AtomList();

       for (Atom atom:atomList){
           out.add(new Atom(new PdbLine(atom.toString())));
       }
       MolecularSystem.setCurrentMolecularSystem(saveMS);
       return out;
     }

    public static void moveHydrogensCloseToHeavy(AtomList atoms, double radius) {
	for (Iterator atomsIter = atoms.iterator(); atomsIter.hasNext();) {
	    Atom atom = (Atom) atomsIter.next();
	    if (atom.type().isHydrogen()) {
		Atom heavy = atom.bonded().get(0);
		atom.randomize(radius, heavy);
	    }
	}
    }

    public static BondParametersList getBondParameters(CommandList commands) {
	Command command = commands.firstWord(PARAMETERS_DIRECTORY); 
	String parametersDirectory = command.secondWord(); 
	return new BondParametersList(parametersDirectory+
				      "/"+BOND_PARAMETERS);
    }
    public static AngleParametersList getAngleParameters(CommandList commands) {
	Command command = commands.firstWord(PARAMETERS_DIRECTORY); 
	String parametersDirectory = command.secondWord(); 
	return new AngleParametersList(parametersDirectory+
				      "/"+ANGLE_PARAMETERS);
    }

    
    private static class DefaultExceptionHandler implements ExceptionHandler {
	public void handle(Exception ex) {
	    System.out.println("An Exception has occured and handled by the  DefaultExceptionHandler.\n"+ex);
	    ex.printStackTrace();
	    throw new RuntimeException(ex);
	}
    }
 
   private static class DoNothingExceptionHandler implements ExceptionHandler {
	public void handle(Exception ex) {
	    System.out.println("An Exception "+ex+" has occured.\n"+"The doNothing ExceptionHandler did nothing.");
	}
    }


}
