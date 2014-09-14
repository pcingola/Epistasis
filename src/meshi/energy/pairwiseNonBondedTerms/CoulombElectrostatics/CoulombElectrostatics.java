package meshi.energy.pairwiseNonBondedTerms.CoulombElectrostatics;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import java.util.*;

/**
 *
 *  An implementation of Coulomb's law.
 *  The fundamental equation of electrostatics is Coulomb's law, 
 *  Describes the force between two point charges :
 * 
 *  Electrostatics(distance_AB) = q1*q2/(dielectricConstant*distance_AB)
 *  dElectrostatics/ddistance_AB = -1*q1*q2/(dielectricConstant*distance_AB^2)
 * 
 *  Where:
 *    1. q1 and q2 are the electrical charges of atoms A and B respectively.
 *    2. A and B are non-covalently-bonded atoms (typically at list 4 bonds apart)
 *    3. distance_AB is the distance between atoms A and B.
 *    4. dielectricConstant is a constant number (selected by the user).
 *
 *
 *  Constructors:
 *     CoulombElectrostatics()
 *     CoulombElectrostatics(DistanceMatrix distanceMatrix,ChargeParametersList  parametersList,
 *			                double weight, double dielectricConstant)
 * 
 *   Object methods
 *      public double evaluate()
 *      public void evaluateAtoms()
 *      public void test(TotalEnergy totalEnergy, boolean testAll, boolean verbose)      
 *
 *
 * 
 **/

public class CoulombElectrostatics extends NonBondedEnergyTerm{
    /**
     * default Constructor
     *
     **/
    public CoulombElectrostatics(){super();}

    /**
     * constructor
     * @param distanceMatrix
     * @param parametersList
     * @param weight
     * @param dielectricConstant
     **/
    public CoulombElectrostatics(DistanceMatrix distanceMatrix, 
				 double weight, double dielectricConstant, double[] charges)  {
        super(toArray(distanceMatrix), weight,distanceMatrix);
        comment = "Electrostatics";
        energyElement = new CoulombElectrostaticEnergyElement(distanceMatrix, 
                                                              weight,dielectricConstant, charges);
    }
}

	
 
