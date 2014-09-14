package meshi.energy.hydrogenBondsAngle;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;
import meshi.util.CommandList;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:47:17
 * To change this template use File | Settings | File Templates.
 */
public class HbondsPunishHOCAngleCreator extends EnergyCreator {

    HbondsPunishHOCAngleEnergy hBondsPunishHOCAngleEnergy;

    HydrogenBondsCreator hydrogenBondsCreator;


    /**
     * @param hydrogenBondsCreator
     **/
    public HbondsPunishHOCAngleCreator(HydrogenBondsCreator hydrogenBondsCreator) {
        super(HYDROGEN_BONDS_ANGLES);
        this.hydrogenBondsCreator = hydrogenBondsCreator;
    }

    /**
     * @param weight
     * @param hydrogenBondsCreator
     **/
    public HbondsPunishHOCAngleCreator (double weight,HydrogenBondsCreator hydrogenBondsCreator)
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
                return new HbondsPunishHOCAngleEnergy(distanceMatrix,
                        HBE.hBondList(),
                        weight(),
                        hydrogenBondsCreator.getSpecialDis() ) ;
          }
          else {
              if (maxAngle == -1){
                return new HbondsPunishHOCAngleEnergy(distanceMatrix,
                          HBE.hBondList(),
                          weight(),
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax);
              }
              else
              return new HbondsPunishHOCAngleEnergy(distanceMatrix,
                          HBE.hBondList(),
                          weight(),
                        hydrogenBondsCreator.getSpecialDis(),
                        xMax,
                        maxAngle);

          }
    }
}
