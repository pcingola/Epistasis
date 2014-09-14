package meshi.applications.prediction.beautify.warpEnergy;

import meshi.util.Updateable;
import meshi.util.UpdateableException;

import java.util.ArrayList;

public class WarpEnergyElementsList extends ArrayList<WarpEnergyElement> implements Updateable {
    private int numberOfUpdates;

    public WarpEnergyElementsList() {
        super();
        numberOfUpdates = 0;
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
        for (WarpEnergyElement element : this) element.update();
    }
}
