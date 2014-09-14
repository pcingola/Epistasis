package meshi.symmetryComplex.energy.cylinder;

import meshi.energy.EnergyElement;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.geometry.Coordinates;

/**
 * Responsible for energy and forces calculations of the Cylinder
 * energy function on a single atom.<BR>
 * The atom's cartesian coordinates are transformed to 
 * <A href="http://en.wikipedia.org/wiki/Cylindrical_coordinate_system">Cylindrical Coordinates</A>
 * and forces are transformed back to cartesian coordinates.
 * 
 * @version 0.1
 * @author Oren Wolfshtat
 */
public class CylinderEnergyElement extends EnergyElement {

    /**
     * The atom which is kept in the cylinder.
     */
    protected final Atom atom;

    protected final double weight;

    /**
     * The cylinder's dimensions.
     */
    protected final double innerR, outerR, height;

    /**
     * The atom's cylindrical coordinates.
     */
    protected double atomR, atomTheta, atomHeight;

    /**
     * Force constants. Currently both equal to 1 for all atoms but can be
     * changed to atom-specific.
     */
    protected final double force , force2;

    /**
     * Constructs an energy elements for the given atom.
     */
    public CylinderEnergyElement(
            Atom atom, CylinderParameters parameters, double weight) {

        this.atom = atom;
        this.height = parameters.height;
        this.innerR = parameters.innerR;
        this.outerR = parameters.outerR;
        this.weight = weight;

        force = 1 * weight;      // replace 1 with a force constant.
        force2 = 2 * 1 * 1 * weight; // replace 1 with a force constant.

        updateCylindricalCoors();

        setAtoms();
        updateFrozen();
    }

    /**
     * Updates the cylindrical coordinates according to the cartesian coordinates.
     */
    private void updateCylindricalCoors() {
        double atomX = atom.x(), atomY = atom.y();
        if (Math.abs(atomX) < VERY_SMALL & Math.abs(atomY) < VERY_SMALL)
            atomX = VERY_SMALL;    //toDo what does it mean
        atomR = Math.hypot(atomX, atomY);
        atomTheta = Math.atan2(atomY, atomX);
        atomHeight = atom.z();
    }

    public void setAtoms() {
        atoms = new AtomList(atom);
    }

    /**
     * Calculates the cylinder energy and updates forces on the atom.
     */
    public double evaluate() {
        if (frozen()) //TODO Check when to call super.updateFrozen()
            return 0.0;
        updateCylindricalCoors();
        double energy = 0;

        // Calculate the atom's position relative to the cylinder.
        double disFromOuterR = atomR - outerR;
        double disFromInnerR = atomR - innerR;
        double disFromTop = atomHeight - height;
        double disFromBottom = atomHeight + height;

        double forceOnR = 0, forceOnHeight = 0;

        if (disFromOuterR > 0) { // The atom's outside the outer radius.
            double disFromOuterR2 = disFromOuterR * disFromOuterR;
            energy += disFromOuterR2 * force;
            forceOnR = -1 * disFromOuterR * force2;
        }
        if (disFromInnerR < 0) { // The atom's inside the inner radius.
            double disFromInnerR2 = disFromInnerR * disFromInnerR;
            energy += disFromInnerR2 * force;
            forceOnR = -1 * disFromInnerR * force2;
        }
        if (disFromTop > 0) { // The atom's above the cylinder.
            double disFromTop2 = disFromTop * disFromTop;
            energy += disFromTop2 * force;
            forceOnHeight = -1 * disFromTop * force2;
        }
        if (disFromBottom < 0) { // The atom's below the cylinder.
            double disFromBottom2 = disFromBottom * disFromBottom;
            energy += disFromBottom2 * force;
            forceOnHeight = -1 * disFromBottom * force2;
        }

        updateForces(forceOnR, forceOnHeight);
        return energy;
    }

    /**
     * Updates the forces on the atom. The parameters are forces in cylidrical
     * coordinates and are transformed to cartesian coordinates.
     */
    private void updateForces(double forceOnR, double forceOnHeight) {
        double forceOnX = forceOnR * Math.cos(atomTheta);
        double forceOnY = forceOnR * Math.sin(atomTheta);
        atom.addToFx(forceOnX);
        atom.addToFy(forceOnY);
        atom.addToFz(forceOnHeight);
    }

    /**
     * For loops building. Should be changed.
     */
    public static boolean isPointInCylinder( Coordinates coordinates, double outerR, double innerR, double height) {
        double atomR = Math.hypot(coordinates.x(), coordinates.y());
        double atomHeight = coordinates.z();

        return (atomR <= outerR) & (atomR >= innerR) & (Math.abs(atomHeight) <= height);
    }
}
