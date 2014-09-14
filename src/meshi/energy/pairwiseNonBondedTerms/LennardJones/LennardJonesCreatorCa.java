package meshi.energy.pairwiseNonBondedTerms.LennardJones;
import meshi.energy.pairwiseNonBondedTerms.*;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.util.*;
import meshi.util.string.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;   
import meshi.geometry.*;
import meshi.util.*;
import meshi.parameters.*;

public class LennardJonesCreatorCa extends LennardJonesCreator{
    private static double[][][] parameters = null;
    public LennardJonesCreatorCa(double weight) {
  	super(weight);
	parametersFileName = LENNARD_JONES_PARAMETERS_CA;
    }
    public LennardJonesCreatorCa() {
  	super();
	parametersFileName = LENNARD_JONES_PARAMETERS_CA;
    }
}
