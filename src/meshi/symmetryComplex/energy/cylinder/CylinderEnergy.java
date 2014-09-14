package meshi.symmetryComplex.energy.cylinder;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;

/**
 * Cylinder energy term. The energy of every atom is w * f * 
 * ([radius term] + [height term]), where:<BR>
 * <UL>
 * <LI>w is the weight.
 * <LI>f is the force constant. Currently f=1 but this can be changed
 * <LI>[radius term] is zero when the atom is between the two (infinite)
 * coaxial cylinders, or the square of the distance to the nearer one.
 * <LI>[height term] is zero when the atom is between the cylinders' (infinite)
 * top and bottom surfaces, or the square of the distance to the nearer one.
 * </UL>
 * 
 * @version 0.1
 * @author Oren Wolfshtat
 */
public class CylinderEnergy extends SimpleEnergyTerm {

	public CylinderEnergy(AtomList atomList,
						  DistanceMatrix distanceMatrix,
						  ParametersList parametersList,
						  double weight,
						  Filter filter) {

		super(toArray(distanceMatrix), parametersList, weight);
		comment = "Cylinder";
		createElementsList(filter == null ?
							atomList :
							atomList.filter(filter));
	}

	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		return new CylinderEnergyElement(
				(Atom) baseElement, (CylinderParameters) parameters, weight);

//		return null;
	}

}