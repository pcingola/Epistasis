package meshi.symmetryComplex.energy;

import meshi.energy.EnergyElement;
import meshi.symmetryComplex.molecularImageElements.SymmetricComplex;

/**
 * This energy element always contributes zero to the total energy,
 * but in every iteration it updates image atoms/proteins... TODO
 */
public class SymmetryEnergyElement extends EnergyElement {

    protected SymmetricComplex symmetricComplex;
    double weight;

    public SymmetryEnergyElement(SymmetricComplex symmetricComplex, double weight) {
        this.symmetricComplex = symmetricComplex;
        this.weight = weight;
    }

    // TODO: update?
    public double evaluate() {return 0.0;}

    protected void setAtoms() {
        atoms = symmetricComplex.atoms();
    }

    public String toString() {
        return("SymmetryEnergyElement, symmetricComplex: "+
            symmetricComplex+", weight: "+weight);
    }
}