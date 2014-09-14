package meshi.symmetryComplex.energy;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;

import java.util.ArrayList;

public class SymmetryEnergy extends SimpleEnergyTerm {

    public SymmetryEnergy(SymmetricComplex symmetricComplex, DistanceMatrix distanceMatrix,
                          ParametersList parametersList, double weight) {
        super(toArray(symmetricComplex, distanceMatrix), parametersList, weight);  // distanceMatrix must be the second updatable resource!!! because it must be updated after SymmetryComplex
        comment = "Symmetry";
        elementsList = new ArrayList();
        EnergyElement newElement = createElement(symmetricComplex, null);
        elementsList.add(newElement); //toDo to check elementList
    }

    public EnergyElement createElement(Object baseElement, Parameters parameters) {

        return new SymmetryEnergyElement((SymmetricComplex)baseElement, weight);
    }

     public void test(TotalEnergy totalEnergy, Atom atom) {
         System.out.println("No energy test for the SymmetryEnergy. This term is used to update only and does not entend any energy calculations.");
     }
}