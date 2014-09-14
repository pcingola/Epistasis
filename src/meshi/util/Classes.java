package meshi.util;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.pairwiseNonBondedTerms.LennardJones.*;

public interface Classes {
    Class BOND_PARAMETERS_CLASS = (new BondParameters()).getClass();    
    Class BOND_ELEMENT_CLASS    = (new BondEnergyElement()).getClass();
    //Class LENNARD_JONES_PARAMETERS_CLASS = (new LennardJonesParameters()).getClass();    
    //Class LENNARD_JONES_ELEMENT_CLASS    = (new LennardJonesEnergyElement()).getClass();
}
