package meshi.applications.prediction.beautify.unWarpEnergy;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;
import meshi.util.UpdateableException;

import java.util.ArrayList;

public class UnWarpEnergyElementsList extends ArrayList<UnWarpEnergyElement> implements Updateable {
    private int numberOfUpdates;

    public UnWarpEnergyElementsList(Protein protein, double weight) {
        super();
        Atom first = protein.atoms().atomAt(0);
        int number = first.residue().number();
        for (int i = 1; i < protein.atoms().size(); i++) {
            Atom second = protein.atoms().atomAt(i);
            if (number + 1 != second.residue().number()) {
                throw new RuntimeException("Very weird " + i + " " + number + " " + first + " " + second);
            } else {
                number++;
                add(new UnWarpEnergyElement(first, second, i * 3.8, weight));
            }
        }
    }

    /*------------------------------------------ update --------------------------------------------*/
    /**
     * Updates the distances.
     */
    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new UpdateableException();
    }

    public void update() {
        for (UnWarpEnergyElement element : this) element.update();
    }
}
