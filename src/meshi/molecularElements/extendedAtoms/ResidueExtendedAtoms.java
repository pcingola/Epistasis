package meshi.molecularElements.extendedAtoms;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;
import meshi.geometry.putH.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;
import java.util.*;
public class ResidueExtendedAtoms extends Residue {
    public static final ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();

    public final Atom H, N, CA, CB, C, O, OXT;    


    protected ResidueExtendedAtoms(String name) {
	super(name);
	H = N = CA = CB = C = O = OXT = null;
    }
    /**
     * Do not use this constructor to instantiate a creator object.
     **/
    public ResidueExtendedAtoms(ResidueType type, AtomList atomList, ResidueIdentifier id, ResidueMode mode) {
	super(type.nameThreeLetters(),type, id, new AtomList(),mode);
	Atom old;
	double tf;

	old = find(atomList,BBatom.CA);
	if (old == null) throw new RuntimeException("CA undefined");
	CA = new Atom("CA",this,type.caType(),new Coordinates(old),old.temperatureFactor());

	old = find(atomList,BBatom.N); 
	if (old == null) tf = new Double(0);
	else tf = old.temperatureFactor();
	N = new Atom("N",this, type.nType(), new Coordinates(old),tf);
	if (mode == ResidueMode.NTER) { 
	    N.setType(AtomType.TRN);
	    H = null;
	}
	else {
	    if (type == ResidueType.PRO) H = null;
	    else {
		old = find(atomList,BBatom.H); 
		if (old == null) tf = new Double(0);
		else tf = old.temperatureFactor();
		H = new Atom("H",this,type.hType(), new Coordinates(old),tf);
	    }
	}
	
	if (type != ResidueType.GLY) {
	    old = find(atomList,BBatom.CB); 
	    if (old == null) tf = new Double(0);
	    else tf = old.temperatureFactor();
	    CB = new Atom("CB",this,type.cbType(), new Coordinates(old),tf);
	}
	else CB = null;

	old = find(atomList,BBatom.C);
	if (old == null) tf = new Double(0);
	else tf = old.temperatureFactor();
        C = new Atom("C",this,type.cType(), new Coordinates(old),tf);

	old = find(atomList,BBatom.O); 
	if (old == null) tf = new Double(0);
	else tf = old.temperatureFactor();
	O = new Atom("O",this,type.oType(), new Coordinates(old),tf);

	if (mode != ResidueMode.CTER) { 
	    OXT = null;
	}
	else {
	    C.setType(AtomType.TRC);
	    O.setType(AtomType.TRO);
	    old = find(atomList,"OXT"); 
	    if (old == null) tf = new Double(0);
	    else tf = old.temperatureFactor();
	    OXT = new Atom("OXT",this, AtomType.TRO, new Coordinates(old),tf);
	}
	addAtomsAndBonds();
    }

    private void addAtomsAndBonds() {
	if (H != null) atoms.add(H);
	atoms.add(N);
	atoms.add(CA);
	if (CB != null) atoms.add(CB);
	atoms.add(C);
	atoms.add(O);
	if (OXT != null) atoms.add(OXT);
	if (H != null) bonds.add(H.bond(N));
	bonds.add(N.bond(CA));
	bonds.add(CA.bond(C));
	bonds.add(C.bond(O));
	if (CB != null) bonds.add(CA.bond(CB));
	if (OXT != null) bonds.add(C.bond(OXT));
    }
	

    public static Atom getAtom(String name, AtomType type, AtomList atomList, Residue residue) {
	Atom old = find(atomList,name);
	double tf;
	if (old == null) tf = new Double(0);
	else tf = old.temperatureFactor();
	return new Atom(name,residue, type, new Coordinates(old),tf);
    }

    public Atom head() {return  C;}
    public Atom tail() {return  N;}

 
    private static String nameConverter(String name) {
	if (name.equals("1HD2")) return "HD21";
	if (name.equals("2HD2")) return "HD22";
	if (name.equals("1HE2")) return "HE21";
	if (name.equals("2HE2")) return "HE22";
	return name;
    }
    
    public final Atom ca() { return CA;}
    public final Atom cb() { return CB;}
    public final Atom c() { return C;}
    public final Atom o() { return O;}
    public final Atom n() { return N;}
    public final Atom h() { return H;}

    public final Atom amideN() {return N;}

    // ***this is Guy's code***

    //a)
    /*      AA replacing methods        */
     
    public Ala toAla()  {
	Arg out ;
	ResidueMode mode = getMode();
	return new Ala(atoms,ID(),mode);
    }

    public Cys toCys()  {
	ResidueMode mode = getMode();
	return new Cys(atoms,ID(),mode);
    }

    public Asp toAsp()  {
	ResidueMode mode = getMode();
	return new Asp(atoms,ID(),mode);
    }

    public Glu toGlu()  {
	ResidueMode mode = getMode();
	return new Glu(atoms,ID(),mode);
    }

    public His toHis()  {
	ResidueMode mode = getMode();
	return new His(atoms,ID(),mode);
    }

    public Ile toIle()  {
	ResidueMode mode = getMode();
	return new Ile(atoms,ID(),mode);
    }

    public Lys toLys()  {
	ResidueMode mode = getMode();
	return new Lys(atoms,ID(),mode);
    }

    public Leu toLeu()  {
	ResidueMode mode = getMode();
	return new Leu(atoms,ID(),mode);
    }

     public Met toMet()  {
	ResidueMode mode = getMode();
	return new Met(atoms,ID(),mode);
    }

   public Asn toAsn()  {
	ResidueMode mode = getMode();
	return new Asn(atoms,ID(),mode);
    }

    public Pro toPro() {
	ResidueMode mode = getMode();
	return new Pro(atoms,ID(),mode);
    }

    public Gln toGln()  {
	ResidueMode mode = getMode();
	return new Gln(atoms,ID(),mode);
    }

    public Arg toArg()  {
	ResidueMode mode = getMode();
	return new Arg(atoms,ID(),mode);
   }

    public Ser toSer() {
	ResidueMode mode = getMode();
	return new Ser(atoms,ID(),mode);
     }

    public Thr toThr() {
	ResidueMode mode = getMode();
	return new Thr(atoms,ID(),mode);
    }

    public Val toVal() {
	ResidueMode mode = getMode();
	return new Val(atoms,ID(),mode);
    }

    public Trp toTrp()  {
	ResidueMode mode = getMode();
	return new Trp(atoms,ID(),mode);
    }

    public Tyr toTyr() {
	ResidueMode mode = getMode();
	return new Tyr(atoms,ID(),mode);
    }

    public void addHydrogens(BondParametersList bondParametersList, 
			     AngleParametersList angleParametersList) {
	if ((H != null) && 
	    H.nowhere() && 
	    (!N.nowhere()) &&
	    (!CA.nowhere())) {
	    try {
		PutHposLog log = PutHpos.pos(H,bondParametersList, angleParametersList);
		if (log.numberOfLegandsWithCoordinatesOfHeavy == 1) {/* The orientation of the hydrogen is random (as the 
									coordinates of the carbonyl atom of the previous residue
									are not set.*/
		    H.resetCoordinates();
		}
	    } catch (NotEnoughBoundAtomsException ex) {H.resetCoordinates();} /* Ignore the exception the class Protein
										   handles missing atoms anyway.*/
	}
    }
}
    //c)
    /*    private static void setBackboneAtoms(ResidueExtendedAtoms residue, ResidueExtendedAtoms in){

        setNewCoordinates(residue.CA,  in.CA.x(),  in.CA.y(),in.CA.z()) ;
        setNewCoordinates(residue.C,  in.C.x(),  in.C.y(),in.C.z())  ;
        setNewCoordinates(residue.O,  in.O.x(),  in.O.y(),in.O.z())  ;
        setNewCoordinates(residue.N,  in.N.x(),  in.N.y(),in.N.z())  ;

        if (in.H!=null && !(residue instanceof Pro)){
	    setNewCoordinates(residue.H,  in.H.x(),  in.H.y(),in.H.z())  ;
	}
            

    }
    //d)
        private Atom[] getSortedAtoms()
    {

	int atomsCounter=0;
	int copyItr=0;
	AtomList aList=this.atoms();
	if (aList==null) return null;

	for (int j=0 ; j<aList.size(); j++)
	    if  (((aList.atomAt(j)).name()).length()>1  &&  !((aList.atomAt(j)).name()).equals("OXT")  && !(((aList.atomAt(j)).name()).charAt(0)=='H')) atomsCounter++;


	Atom[] sortedList= new Atom[atomsCounter];
	for (int j=0 ; j<aList.size(); j++)
	    {//copy the pointers of the relevant atoms to a new list
		if  (((aList.atomAt(j)).name()).length()>1   &&  !((aList.atomAt(j)).name()).equals("OXT")    && !(((aList.atomAt(j)).name()).charAt(0)=='H'))
		    {
			sortedList[copyItr]=aList.atomAt(j);
			copyItr++;
		    }
	    }
	sortList(sortedList);
	return sortedList;

    }  //end method getSortedAtoms

    //e)
    private static int mapMatchingAtoms(Atom[] oldList,Atom[] newList){
	int newDist=1, oldDist=1;
    	int oldPtr=0,newPtr=0;

    	while (newPtr<newList.length && oldPtr<oldList.length )
	    {
	
		if ( (newPtr>0 && changedDistance(newList,newPtr)) && oldDist>=newDist)  newDist++;
		if (oldPtr>0 && changedDistance(oldList,oldPtr))  oldDist++;

		if (oldDist>=newDist )
		    {
        		newList[newPtr].setX(oldList[oldPtr].x());
        		newList[newPtr].setY(oldList[oldPtr].y());
        		newList[newPtr].setZ(oldList[oldPtr].z());
       		 	newPtr++;
		    }

		oldPtr++;
	    } //end while
    	
    	return newPtr;
			
    }
	
    //f)
    private static Atom[] listAddN(Atom[] newList, ResidueExtendedAtoms residue) {
			
        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.N;
        
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    private static Atom[] listAddO(Atom[] newList, ResidueExtendedAtoms residue) {

        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.O;
		
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    private static Atom[] listAddC(Atom[] newList, ResidueExtendedAtoms residue) {

        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.C;
		
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    //g)
    private static void setNewCoordinates(Atom atom, double X, double Y,double Z){
	
	atom.setX(X);
	atom.setY(Y);
	atom.setZ(Z);

    }
		
		
		
    private static Atom	getNextAtom(int newPtr, Atom[] newList,ResidueExtendedAtoms residue){
		
	double centerX,centerY,centerZ;       
	boolean condition3;		
	Atom chosenAtom,result=null;
	int modelAtom;
   
        if (newPtr>3 &&  !changedDistance(newList,newPtr) ) modelAtom=newPtr-2;

        else modelAtom=newPtr-1;
         
        centerX=newList[modelAtom].x() ;
        centerY=  newList[modelAtom].y() ;
        centerZ= newList[modelAtom].z() ;
        double minDistance=1.45,maxDistance=1.55;
        int itr=0;
	double DistFromAlpha=0;;
	for (int i=0;i<1000;i++){

	    double newDistFromAlpha;
						          
     	    do{

		chosenAtom=createNewAtom(centerX,centerY,centerZ,residue,maxDistance);//new Atom(centerX,centerY,centerZ,1.55,"locationGen",residue,1);
          
		itr++;
		if (itr==3000) { minDistance-=0.1; itr=0;}
		if (itr==1000) { maxDistance+=0.1;}
    		condition3=(minDist(chosenAtom.x(),chosenAtom.y(),chosenAtom.z(),newList,newPtr)< minDistance);

	    } while (condition3);
          
	    newDistFromAlpha=distance(chosenAtom.x(),chosenAtom.y(),chosenAtom.z(),newList[3].x(),newList[3].y(),newList[3].z());
	    if (newDistFromAlpha>DistFromAlpha) {DistFromAlpha=newDistFromAlpha;result=chosenAtom;}

	}		  
          

	return result;
		
    }
    //h)
    private final static double BACKBONE_RELIABILITY=200.0; 
    
    //i)
    private void sortList(Atom[] aList)   //performs insertion sort
    {
	if (aList==null) return;
	for (int i=1; i<aList.length; i++)
	    insert(aList,i);

    }

    //j)

    private static void insert(Atom[] aList, int i)
    {
	Atom value=aList[i];
	while (i>0 && compareAtoms(aList[i-1],value)>0)
	    {
		aList[i]=aList[i-1];
		i=i-1;
	    }

	aList[i]=value;
    }

    //k)

    private static boolean changedDistance(Atom[] list, int ptr){

	String firstAtomName= list[ptr].name();
	String secondAtomName = list[ptr-1].name();
	if (firstAtomName.charAt(1)!=secondAtomName.charAt(1) ) return true;
	return false;

    }

    //l)

		
    private static Atom createNewAtom(double centerX,double centerY,double centerZ,ResidueExtendedAtoms residue,double radius){

	/*      double atomsMaxDistance=1.55;

	double xDirection= ((atomCreationRand.nextInt()%1000)+1);
	double yDirection= ((atomCreationRand.nextInt()%1000)+1);
	double zDirection= ((atomCreationRand.nextInt()%1000)+1);
	double sum= Math.sqrt(Math.pow(xDirection,2)+Math.pow(yDirection,2)+Math.pow(zDirection,2));
	xDirection/=sum;
	yDirection/=sum;
	zDirection/=sum;
	xDirection*= atomsMaxDistance;
	yDirection*= atomsMaxDistance;
	zDirection*= atomsMaxDistance;                                             */
    /*	return new Atom(centerX,centerY,centerZ,radius,"random",residue,1);
	//      return new Atom(centerX+xDirection,centerY+yDirection,centerZ+zDirection,"random",residue,1);


    }
    //n)

    private static double minDist(double x,double y,double z,Atom[] list,int place)
    {
	double min=10000;
	for (int i=0; i<place; i++)
	    {
		double dist=  distance(x,y,z,list[i].x(),list[i].y(),list[i].z());
		if (min>dist) min=dist;
	    }
	return min;
    }
    //o)
    private static double distance(double x1,double y1,double z1,double x2,double y2,double z2)
    {
        return Math.sqrt((x1-x2)*(x1-x2)+(y1-y1)*(y1-y2)+(z1-z2)*(z1-z2));
    }
    //p)
    private static int compareAtoms(Atom a1,Atom a2)
    {
	if ((a1==null) || (a2==null)) return 0;
	if ((a1.name()).equals(a2.name())) return 0;

	char c1,c2;
	char[] cArray={'A','B','G','D','E','Z','H'};
	char[] dArray={'1','2','3'};
	c1=(a1.name()).charAt(1);
	c2=(a2.name()).charAt(1);

	if (c1==c2)  //same distance from CA, compare positions( 1 or 2 )
	    {
		for   (int i=0; i<dArray.length ; i++)
		    {
			if ((a1.name()).charAt(2)==dArray[i]) return -1;
			if ((a2.name()).charAt(2)==dArray[i]) return 1;
		    }
	    }

	else
	    {
		for   (int i=0; i<cArray.length ; i++)
		    {
			if (c1==cArray[i]) return -1;
			if (c2==cArray[i]) return 1;

		    }
	    }



	return 0;
	}*///end method compareAtoms 

    //}
    
   
