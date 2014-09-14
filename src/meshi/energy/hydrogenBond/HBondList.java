package meshi.energy.hydrogenBond;

import java.util.*;
import java.lang.Integer;

import meshi.geometry.*;
import meshi.util.Updateable;
import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.*;

/**
 * @author amilev
 *   Note: if the 2 atoms of the HB are frozen - we don't add them to the list.
 **/
public class HBondList extends HBdistanceList implements Updateable {

    //--------------------------------- data fields -----------------------------
     public int insertionToList = 0;
    public int deletionFromList = 0;
    
    /*
     * The  PairsList updated every X steps.
     */
    private final int UPDATE_EVERY_X_STEPS = 50;
    ///*
    // * List of the candidate to hydrogen bonding given by the DistanceMatrix object
    // */
    //protected DistanceList nonBondedList;
    /*
     * List of all the new HB elements that were added to the  hBondList in the last update call.
     */
    protected DistanceList newhBondList;
    protected static HBdistanceList inputNewHBList = new HBdistanceList(new GoodResiduesForHB());
    /*
    * List of all the HB elements
    */
    //protected DistanceList  hBondList;

    /*
     * holds the relevant row numbers
     */
    private ArrayList<Integer> relevantRows;
    /*
     * get the parameters (epsilon, sigma) of each element
     */
    private HydrogenBondsParametersList parametersList;
    private DistanceMatrix distanceMatrix;

    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates =0;
    /*
     *  hBondList/this should be updated when countUpdates >= UPDATE_EVERY_X_STEPS
     */
    private int countUpdates = 50;
    /*
     * Filter for update
     */
    private GoodResiduesForHB goodResiduesForHB;
    private final IsAlive isAlive;

    /**
     * to find easily the last element that was in this list before the new update
     */
    private int oldSize = -1;

   public final int getOldSize() { return oldSize; }
    public final DistanceList newhBondList() { return newhBondList;}
    public static HBdistanceList inputNewHBList() { return inputNewHBList;}
    public final int countUpdates() {return countUpdates;}

    //-------------------------------- constructors --------------------------------

    /**
	 * @param distanceMatrix
	 */
	public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList) {
        super(new GoodResiduesForHB());
        goodResiduesForHB = new GoodResiduesForHB();
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        //hBondList = new DistanceList();
        newhBondList = new DistanceList(100);
        relevantRows = new ArrayList<Integer>();
        //nonBondedList = distanceMatrix.nonBondedList();
        //update(nonBondedList);
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
    }

    public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList,
                     DistanceList specialDistances) {
        super(new GoodResiduesForHB(specialDistances));
        goodResiduesForHB = new GoodResiduesForHB(specialDistances);
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        newhBondList = new DistanceList(100);
        relevantRows = new ArrayList<Integer>();
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
    }

    //--------------------------------------- methods ---------------------------------------

	/* (non-Javadoc)
	 * @see meshi.util.Updateable#update(int)
	 */
	public void update(int numberOfUpdates) {
		if (numberOfUpdates == this.numberOfUpdates+1) {
		    this.numberOfUpdates++;
            //System.out.println("HBondList: in update(int numberOfUpdates):");
		    update();
		}
		else if (numberOfUpdates != this.numberOfUpdates)
		    throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
	}


    private void update() {
        oldSize = -1;
        newhBondList.clear();
        if (countUpdates == UPDATE_EVERY_X_STEPS) {//TODO can delete this if - i leave it for now for possible future needs
            int prevrousSize = size() ;
            updateTwoLists(inputNewHBList);
            int secondSize = size() ;
            countUpdates = 0;
            cleanList_dm();
            oldSize = size() - inputNewHBList .size() ;
            //System.out.println("HBondList: prevouse size "+prevrousSize +" secondSize  "+secondSize +" current size "+size());
            if(size() > secondSize | secondSize < prevrousSize | size() < oldSize )
                throw new RuntimeException("HBondsList: problem at update - the list after adding elements is shorter or after cleaning is longer");
        }
        else if (countUpdates > UPDATE_EVERY_X_STEPS)
            throw new RuntimeException("Something weird with HbondList.update()\n"+
                                       "countUpdates = "+countUpdates+" UPDATE_EVERY_X_STEPS  = "+UPDATE_EVERY_X_STEPS);
        else {
            //System.out.println("HBondList: in update: distanceMatrix.newHydrogenBondsList():");
            //distanceMatrix.newHydrogenBondsList().print();
            countUpdates++;
           updateTwoLists(inputNewHBList);
           cleanList_dm();
           oldSize = size() - inputNewHBList .size() ;
           if(size() < oldSize)
                     throw new RuntimeException("HBondList: size < old size");
            }

    }

	/**
	 * @param atomPairList of new distances that where added after  HbondList was updated in the last time
	 */
    private void updateTwoLists(HBdistanceList atomPairList) {
        oldSize = -1;
        oldSize = size();
        for (Distance pair:atomPairList){
            if (!isAlive.accept(pair))
                  System.out.println("HBondList: there is a problem");
            //System.out.println("HBondList: in updateTwoLists: pair: "+pair);
            //if (goodResiduesForHB.accept(pair)) {     this term is checked on HBdistanceList
            //System.out.println("!!!!!! HBondList: in updateTwoLists: Goodpair !!!!! : "+pair);
                HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                distanceAttribute .set(pair.atom1() ,pair .atom2() );
                pair.addAttribute(distanceAttribute);
		if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("Weird pair in  updateTwoLists "+pair);
                add(pair);
                newhBondList.add(pair);
                insertionToList ++;
             //}
        }
    }

    /*
     * the update is done by first go over the heads atoms of the row in distance matrix
     * and just then (only if the head atom is Hydrogen or Oxygen) go over the rows themself.
     * (more afficient ?)
     */
    private void update_dm(){
        // System.out.println("HBondList: in update_dm: start");
        Atom headAtom;
        MatrixRow matrixRow;
        Iterator rows = distanceMatrix.rowIterator();
        //System.out.println("HBondList: dm  length = " +distanceMatrix .matrix.length );
        HB_AtomAttribute headAttribute;
        if (relevantRows.isEmpty()){
            // System.out.println("HBondList: relevantRows is empty");
            while (rows.hasNext() ){
              if ((matrixRow = (MatrixRow) rows.next()) != null){
            //while((matrixRow = (MatrixRow) rows.next()) != null){
                headAtom = matrixRow.atom.atom;
                //System.out.println("HBondList: headAtom = " +headAtom);
                headAttribute = (HB_AtomAttribute) headAtom.getAttribute(HB_AtomAttribute.key);
                if (headAttribute != null){//meens that it is H or O !
                    //System.out.println("HBondList: headAttribute != null");
                    relevantRows.add(headAtom.number());//add the row number
		    for (Distance pair:matrixRow) {
			if (pair.distance() < Distance.INFINITE_DISTANCE) {
			    if(goodResiduesForHB.accept(pair) ){
				if(!isAlive .accept(pair ))
				    throw new RuntimeException("need the second check !\n"+
							       pair+"\n"+
							       pair.atom1()+" "+pair.atom1().frozen()+
							       pair.atom2()+" "+pair.atom2().frozen());
				HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
				distanceAttribute.setParameters
				    ((HydrogenBondsParameters) parametersList.parameters(pair));
				distanceAttribute .set(pair.atom1() ,pair .atom2() );
				pair.addAttribute(distanceAttribute);
				//hBondList.add(pair);
				if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("Weird pair in update_dm  "+pair);
				add(pair);
				insertionToList ++;
			    }
			}
                    }//while
                }//if
              }//if
            }//while there is more rows
        }//if not empty
        else {

            int row;
            for (Integer current:relevantRows){
                row = current;
                matrixRow = distanceMatrix.rowNumber(row);
                headAtom = matrixRow.atom.atom;
		for (Distance pair:matrixRow) {
		    if (pair.distance() < Distance.INFINITE_DISTANCE) {
			if(goodResiduesForHB.accept(pair)){
			    if(!isAlive .accept(pair ))
                                throw new RuntimeException("need the secod check !");
			    HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
			    distanceAttribute.setParameters
				((HydrogenBondsParameters) parametersList.parameters(pair));
			    pair.addAttribute(distanceAttribute);
			    // hBondList.add(pair);
			    if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("Weird pair in update_dm #2  "+pair);			    
			    add(pair);
			    insertionToList ++;
			}
		    }
                }//while
            }//while there is more rows
        }//else
        // System.out.println("HBondList: in update_dm: end");
        // this.print();
    }


    public Iterator newWithinRmaxIterator() {
        return new WithinRmaxIterator(newhBondList);// WithinRmaxIterator is an internal class.
    }

    /**
     * save only the live elements in this list
     */
    public void cleanList_dm(){
	Distance[] newInternalArray = new Distance[size()];
	int currentIndex =0;
	for(Distance pair:this){
	    if (isAlive .accept(pair )){
		newInternalArray [currentIndex] = pair ;
		currentIndex++;
               }
	    else
		deletionFromList++;
	}
	clear();
	for (Distance pair:newInternalArray) {
	    //if ((pair.mode() == DistanceMode.INFINITE)) throw new RuntimeException("cleanList_dm  "+pair);
	    if (pair != null) add(pair);
	}
    }

       public Iterator withinRmaxIterator() {
           return new WithinRmaxIterator(this);// WithinRmaxIterator is an internal class.
       }


    //--------------------------- internal class IsWithInRmax ---------------------------

    static class IsWithInRmax implements Filter{
       private final double rMax;
        public IsWithInRmax(){
            super();
            rMax =   DistanceMatrix.rMax();
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            //double dis = rMax - distance.distance();
            return ((rMax >= distance.distance()) & (!distance.dead()));
        }
    } //--------------------------- internal class IsAlive ---------------------------

    static class IsAlive implements Filter{
        DistanceMatrix dm;
        boolean firstTimeWornning = true;
        public IsAlive(DistanceMatrix  matrix){
            super();
            dm = matrix ;
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            boolean frozen = distance.atom1() .frozen() & distance.atom2() .frozen() ;
            double dis = distance.distance();
            if (distance.dead()) return false;
            //boolean alive = distance .update(rMax2 ,rMaxPlusBuffer2 );
            if (frozen){
                if(firstTimeWornning )   //TODO check this - maybe frozen atoms should be looked at ?
                {
                    System.out.println("*** NOTE: frozen atoms does not get into HBondList !!! *****");
                    firstTimeWornning = false;
                }

                return false;
            }
            return dis != Distance.INFINITE_DISTANCE;
        }
    }
    


    //--------------------------- internal class WithinRmaxIterator ---------------------------
    
    private  class WithinRmaxIterator implements Iterator   {
        IsWithInRmax isWithInRmax= new IsWithInRmax();
        DistanceList list;
	int current;
	int listSize = 0;
        public WithinRmaxIterator(DistanceList list){
	    this.list = list;
	    listSize = list.size();
	    current = 0;
        }

        /*
         * @return the next element that itas 2 atoms are within Rmax or Null if there is no such element
         */
        public Object next() {
	    if (current >= listSize) return null;
	    Object obj = list.get(current);
	    current++;
	    if (isWithInRmax.accept(obj)) return (obj);
	    else return next();
        }
	public boolean hasNext() {return (current < listSize);}
	public void remove() {throw new RuntimeException("do not do that");}
    }    
    
}




