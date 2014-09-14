package meshi.symmetryComplex.energy.edmEnergy;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.util.filters.Filter;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;

public class EDMEnergy extends SimpleEnergyTerm {

    public EDMEnergy( Protein protein, EDMEnergyParametersList eepl,
                      double weight, String comment, Filter filter ) {
        super( toArray(), eepl, weight );

        this.comment = comment;

        createElementsList( filter == null ?
                            protein.atoms() :
                            protein.atoms().filter(filter));
    }

    public EnergyElement createElement(
            Object baseElement, Parameters parameters) {
        return new EDMEnergyElement(
                (Atom) baseElement, (EDMEnergyParameters) parameters, weight );
    }

}
