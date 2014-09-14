package meshi.applications.prediction.beautify.unWarpEnergy;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;
import meshi.util.UpdateableException;

public class UnWarpEnergy extends AbstractEnergy implements Updateable {
    private static UnWarpEnergyElementsList elements;

    public UnWarpEnergy(double weight, Protein protein) {
        super(toArray(elements = new UnWarpEnergyElementsList(protein, weight)), weight);
        comment = "UnWarpEnergy";
    }

    public void update() {
    }

    public void update(int i) {
        try {
            elements.update(i);
        } catch (UpdateableException ex) {
            throw new RuntimeException("Update problem in UnWarpEnergy:\n" + ex);
        }
    }

    public double evaluate() {
        if (!on) return 0;
        double e = 0;
        for (UnWarpEnergyElement element : elements)
            e += element.evaluate();
        return e;
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (!on) System.out.println("" + this + " is off");
        for (UnWarpEnergyElement element : elements) {
            if (!element.frozen()) {
                element.test(totalEnergy, atom);
            }
        }
    }

    public void evaluateAtoms() {
        if (on) {
            for (UnWarpEnergyElement element : elements)
                element.evaluateAtoms();
        }
    }
}

