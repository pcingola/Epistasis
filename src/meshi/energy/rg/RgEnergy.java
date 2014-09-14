package meshi.energy.rg;

import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;

public class RgEnergy extends CooperativeEnergyTerm {

	public RgEnergy() {
	}

	public RgEnergy(AtomList atomList, DistanceMatrix dm, double weight) {
		super(toArray(), atomList, dm, null, weight);
		comment = "RG";
	}

	@Override
	public double evaluate() {
		if (!on) return 0;
		double energy = 0;
		double Rg = 0;
		double dRg = 0;
		double cmx = 0, cmy = 0, cmz = 0; // center of mass x, y and z
		Atom atom;
		Iterator iter = atomList.iterator();
		while ((atom = (Atom) iter.next()) != null) {
			cmx += atom.x();
			cmy += atom.y();
			cmz += atom.z();
		}
		cmx /= atomList.size();
		cmy /= atomList.size();
		cmz /= atomList.size();
		iter = atomList.iterator();
		while ((atom = (Atom) iter.next()) != null)
			Rg += (cmx - atom.x()) * (cmx - atom.x()) + (cmy - atom.y()) * (cmy - atom.y()) + (cmz - atom.z()) * (cmz - atom.z());
		Rg /= atomList.size();

		energy = Math.sqrt(Rg);
		dRg = 1 / energy / atomList.size(); // Multiplication by 2 was made to save time
		energy *= weight;
		dRg *= weight;

		iter = atomList.iterator();
		while ((atom = (Atom) iter.next()) != null) {
			if (!atom.frozen()) {
				atom.addToFx(-dRg * (atom.x() - cmx)); // Negating so that it is force
				atom.addToFy(-dRg * (atom.y() - cmy)); // Negating so that it is force
				atom.addToFz(-dRg * (atom.z() - cmz)); // Negating so that it is force
			}
		}
		return energy;
	}

	@Override
	public void evaluateAtoms() {
	}

	@Override
	public void update(int numberOfUpdates) {
	}
}
