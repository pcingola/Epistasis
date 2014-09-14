/*
 * Created on 26/01/2005
 * as part of meshi.1.5
 * 
 */
package meshi.energy.hydrogenBondsAngle;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

/**
 * @author amilev
  */
public class HBondsPunishOHNAngleCreator extends EnergyCreator {
	
	HBondsPunishOHNAngleEnergy hBondsPunishOHNAngleEnergy;
	HydrogenBondsCreator hydrogenBondsCreator;
	
	
	/**
	 * @param hydrogenBondsCreator
	 **/
	public HBondsPunishOHNAngleCreator(HydrogenBondsCreator hydrogenBondsCreator) {
		super(HYDROGEN_BONDS_ANGLES);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}
	
	/**
	 * @param weight
	 * @param hydrogenBondsCreator
	 **/
	public HBondsPunishOHNAngleCreator (double weight,HydrogenBondsCreator hydrogenBondsCreator)
	{
		super(weight);
		this.hydrogenBondsCreator = hydrogenBondsCreator;
	}
	
	/**
	 * @param protein
	 * @param distanceMatrix
     * @param commands - NOT IN USE !
	 **/
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
		      CommandList commands)
	{
		HydrogenBondsEnergy HBE = hydrogenBondsCreator.getHydrogenBondsEnergy();
        double xMax;
                double maxAngle;
                try{
                    xMax  = commands.firstWordFilter(HYDROGEN_BONDS_ANGLES ).secondWord(MAX_DISTANCE ).thirdWordDouble() ;
                }
                catch (RuntimeException r){
                       xMax = -1;
                }
                try{
                    maxAngle  = commands.firstWordFilter(HYDROGEN_BONDS_ANGLES ).secondWord(MAX_ANGLE ).thirdWordDouble() ;
                }
                catch (RuntimeException r){
                       maxAngle = -1;
                }

          if (xMax == -1){
        return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                HBE.hBondList(),
                weight(),
                 hydrogenBondsCreator.getSpecialDis() ) ;
               }
          else {
              if (maxAngle == -1){
                   return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        weight(),
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax) ;
              }
              else
              return new HBondsPunishOHNAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        weight(),
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax,
                        maxAngle);

          }
    }

}
