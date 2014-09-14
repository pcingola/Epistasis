package meshi.applications.prediction.beautify.unWarpEnergy;

import meshi.energy.EnergyElement;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;


public class UnWarpEnergyElement extends EnergyElement {
    private FreeDistance distance;
    private Atom first, second;
    private double weight;
    private boolean on = true;
    private double targetDistance;

    public UnWarpEnergyElement(Atom first, Atom second, double targetDistance, double weight) {
        this.first = first;
        this.second = second;
        setAtoms();
        this.weight = weight;
        this.targetDistance = targetDistance;
	distance = null;
    }

    public void update() {
	if (distance == null) distance = new FreeDistance(first, second);
        if (on) distance.update();
    }

    public double evaluate() {
        if (!on) return 0;
        double targetMinusDis = targetDistance - distance.distance();
        double targetMinusDisPlus1 = targetMinusDis + 1;
        double energy = targetMinusDis * targetMinusDis * weight / targetMinusDisPlus1;
        double dEnergyDdistance = -weight * (targetMinusDis / targetMinusDisPlus1) * (2 - targetMinusDis / targetMinusDisPlus1);
        first.addToFx(-dEnergyDdistance * distance.dDistanceDx());
        first.addToFy(-dEnergyDdistance * distance.dDistanceDy());
        first.addToFz(-dEnergyDdistance * distance.dDistanceDz());
        second.addToFx(dEnergyDdistance * distance.dDistanceDx());
        second.addToFy(dEnergyDdistance * distance.dDistanceDy());
        second.addToFz(dEnergyDdistance * distance.dDistanceDz());
        return energy;
    }

    protected void setAtoms() {
        atoms = new AtomList();
        atoms.add(first);
        atoms.add(second);
    }

    public void off() {
        on = false;
    }
}
