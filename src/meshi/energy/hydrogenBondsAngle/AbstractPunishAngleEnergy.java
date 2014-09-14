package meshi.energy.hydrogenBondsAngle;

import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HBondList;
import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.HydrogenBondsParameters;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.util.MeshiException;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:13:56
 * To change this template use File | Settings | File Templates.
 */
public abstract class AbstractPunishAngleEnergy extends NonBondedEnergyTerm {
    protected HBondList hBondList;
    //private DistanceList specialDisatnces = null;

    public AbstractPunishAngleEnergy(){}

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              weight,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement();
    }

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      double xMax)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              weight,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement(xMax);
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      double xMax,
                                      double maxAngle)
    {
        super(toArray(distanceMatrix,hBondList),//updatabel resources
              weight,
              distanceMatrix);
        this.hBondList = hBondList;
        setComment();
        setEnergyElement(xMax,maxAngle);
    }

    public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      DistanceList specialDis){
        this(distanceMatrix,hBondList,weight);
        //specialDisatnces = specialDis;
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      DistanceList specialDis,
                                      double xMAx){
        this(distanceMatrix,hBondList,weight,xMAx);
        //specialDisatnces = specialDis;
    }

     public AbstractPunishAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      DistanceList specialDis,
                                      double xMax,
                                      double maxAngle ){
        this(distanceMatrix,hBondList,weight,xMax,maxAngle);
        //specialDisatnces = specialDis;
    }

    /**
     * set the energyElement to point on the relavent instance
     */
    public abstract void  setEnergyElement();

    /**
        * set the energyElement to point on the relavent instance
        */
       public abstract void  setEnergyElement(double xMax);

    /**
            * set the energyElement to point on the relavent instance
            */
           public abstract void  setEnergyElement(double xMax,double maxAngle);

    /**
     * set the comment to be the relavent comment according to the instance
     */
    public abstract void  setComment();

    public double evaluate() {
        if (! on) return 0.0;
        double energy = 0;
        double element_energy;

        //go over all the relevent elements in hBondList
        Iterator hBondListIter = hBondList.withinRmaxIterator();           //TODO i think the "withinRmaxIterator()" is not needed
        Distance pair;

    /*    if(specialDisatnces != null){
            Iterator specDisIter = specialDisatnces.iterator() ;
            while((pair  = (Distance) specDisIter .next()) != null){
                HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                //distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                distanceAttribute .set(pair.atom1() ,pair .atom2() );
                pair.addAttribute(distanceAttribute);
                    energyElement .set(pair );
                    element_energy = energyElement.evaluate(50*weight); //TODO
                    energy += element_energy;
            }
        }*/
       /* Distance d1=null,d2=null,d3=null,d4=null;
        if(specialDisatnces != null){
           d1 = (Distance)specialDisatnces.fastElementAt(0);
           d2 = (Distance)specialDisatnces.fastElementAt(1);
           d3 =  (Distance)specialDisatnces.fastElementAt(2);
           d4 =  (Distance)specialDisatnces.fastElementAt(3);

        }*/
        while ((pair  = (Distance)  hBondListIter.next()) != null) {
	    if (pair.largeType == AtomType.YOH) continue; //Ugly patch by Chen
            energyElement.set(pair);
            /*boolean isSpecial = false;
	      if(specialDisatnces != null && !specialDisatnces.isEmpty() ){
	      Distance speDis;
	      Iterator specDisIter = specialDisatnces.iterator() ;
	      
	      while(!isSpecial&  (speDis  = (Distance) specDisIter .next()) != null){
	      isSpecial = pair.amiEquals(speDis );
	      }
	      }*/
	    // if(d1 != null && d2 !=null && (pair.amiEquals(d1) | pair.amiEquals(d2) | pair.amiEquals(d3) | pair.amiEquals(d4) ) ){
            /*if(isSpecial)   {
                     element_energy = energyElement.evaluate(weight);
            }
            else*/
	    try {
		element_energy = energyElement.evaluate();
	    }
	    catch (RuntimeException ex) {
		System.out.println("\n Failed to evaluate "+pair+"\n");
		throw ex;
	    }
            energy += element_energy;
            if (element_energy < 0) //energy should always be positive (if the angle is bad angle) or zero if it is good angle or pair of hidrogen-oxygen that are not connected
                throw new MeshiException(comment +": energy ("+energy +") should always be positive or zero: "+pair );

          //TODO add freeElement
        }
        return energy;
    }

    public void evaluateAtoms() {
        if (! on) return;
        Iterator hBondListIter = hBondList.withinRmaxIterator();
        Distance pair;
        while ((pair  = (Distance)  hBondListIter.next()) != null) {
            energyElement.set(pair);
            energyElement.evaluateAtoms();
            //TODO add freeElement
        }

    }

    public void test(TotalEnergy totalEnergy, Atom atom){
        System.out.println("Start Test "+comment);
        if (! on) System.out.println(""+this +" is off");
        Iterator hBondListIter = hBondList.withinRmaxIterator();
        Distance pair;
        while ((pair  = (Distance)  hBondListIter.next()) != null) {
        energyElement.set(pair);
        energyElement.test(totalEnergy,atom);
            //TODO add freeElement
    }
        System.out.println("End Test "+comment);
    }
}
