package meshi.applications.prediction.beautify.warpEnergy;

import meshi.applications.prediction.beautify.BeautifyAttribute;
import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.util.Updateable;
import meshi.util.UpdateableException;

import java.util.Iterator;

public class WarpEnergy extends AbstractEnergy implements Updateable {
    private static WarpEnergyElementsList elements;

    public WarpEnergy(double weight) {
        super(toArray(elements = new WarpEnergyElementsList()), weight);
        comment = "WarpEnergy";
    }

    public void update() {
    }

    public void update(int i) {
        try {
            elements.update(i);
        } catch (UpdateableException ex) {
            throw new RuntimeException("Update problem in WarpEnergy:\n" + ex);
        }
    }

    public void ignoreLoopResidues() {
        for (WarpEnergyElement element : elements) {
            if (BeautifyAttribute.isLoop(element.residue())) element.off();
        }
    }

    public void ignoreUnhappyResidues(double threshold) {
        for (WarpEnergyElement element : elements)
            element.offIfUnhappy(threshold);
    }

    public double evaluate() {
        if (!on) return 0;
        double e = 0;
        for (WarpEnergyElement element : elements)
            e += element.evaluate();
        return e;
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (!on) System.out.println("" + this + " is off");
        for (WarpEnergyElement element : elements) {
            if (!element.frozen()) {
                element.test(totalEnergy, atom);
            }
        }
    }

    public void evaluateAtoms() {
        if (on) {
            for (WarpEnergyElement element : elements)
                element.evaluateAtoms();
        }
    }

    public void add(ResidueAlignmentColumn column) {
        Residue resBeautify = (Residue) column.cell0().obj;
        Residue resShotgun = (Residue) column.cell1().obj;
        if ((!resShotgun.dummy()) &&
                (!resBeautify.dummy())) {
            Atom caBeautify = resBeautify.ca();
            Atom caShotgun = resShotgun.ca();
            elements.add(new WarpEnergyElement(caShotgun, caBeautify, weight));
        }
    }

    public void add(ResidueAlignment residueAlignment) {
        for (Iterator columns = residueAlignment.iterator(); columns.hasNext();)
            add((ResidueAlignmentColumn) columns.next());
    }
}

