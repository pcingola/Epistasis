/**
 * Lennard-Jones potential
 * A classic pair potential with an attractive van der Waals potential tail 
 * and a stiff repulsive core. It is simply (12-6 form): 
 * 
 * E(distanceAb) = 4*epsilonAb{(sigmaAb/distanceAb)^12-(sigmaAb/distanceAb)^6}
 *                = 4*epsilonAb*sigma^6*{(sigma^6*distanceAb^-12) - distanceAb^-6}
 * dE/ddistanceAb = 4*epsilonAb*sigmaAb^6*{(-12*sigmaAb^6)*(distanceAb^-13)+6*(distanceAb^-7)}
 * 
 * Where:
 *    1. a and b are non-covalently-bonded atoms (typically at list 4 bonds from one another)
 *    2. distanceAb is the distance between atoms a and b.
 *    3. epsilonAb and sigmaAb are constants depending on the nature of the atoms.
 *
 * Typically only homo-atom parameters are supplied by the forcefield and hetero-atom parameters 
 * are calculated by the Lorentz-Berthelot mixing rule:
 *
 * epsilonAb = sqrt(epsilonAa * epsilonBb)
 * sigmaAb = (sigmaAa + sigmaBb) /2
 **/
package meshi.energy.pairwiseNonBondedTerms.LennardJones;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.filters.Filter;
import java.util.*;
public class LennardJones extends NonBondedEnergyTerm{
    public LennardJones(){super();}

    public LennardJones(DistanceMatrix distanceMatrix, 
			double weight, double[][][] parameters){
	super(toArray(distanceMatrix), weight,distanceMatrix);
	comment = "LennardJones";
	energyElement = new LennardJonesEnergyElement(distanceMatrix,weight,parameters);
    }
}
