package meshi.applications.prediction.homology;

import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*; 
import meshi.molecularElements.atoms.*;
import meshi.util.filters.Filter;

import java.util.Iterator;

public class NonFrozensRadiusFilter implements Filter {

    AtomList /*atoms,*/nonFrozens;
    double radius;
    int bondDepth;

    public NonFrozensRadiusFilter(AtomList atoms,double radius,int bondDepth){

        /*this.atoms = atoms;*/
        this.radius = radius;
        this.bondDepth = bondDepth;

        nonFrozens = new AtomList();
        Atom atom;
        Iterator atomListIter = atoms.iterator();
        while ((atom = (Atom)atomListIter.next())!=null){
            if (!atom.frozen()) nonFrozens.add(atom);
        }
    }

    public boolean accept(Object obj) {
        Atom atom = (Atom)obj;
        /*if (!atoms.contains(atom))
        throw new RuntimeException("atom "+atom+" is not in the filter list");
*/
        if (!atom.frozen()) return true;

        if (isConnectedToNonFrozen(atom,bondDepth)) return true;

        if (inRadiusFromNonFrozen(atom)) return true;

        return false;
    }

    private boolean isConnectedToNonFrozen(Atom atom, int depth) {

        if (!atom.frozen()) return true;
        if (depth == 0) return false;

        Atom bonded;
        Iterator bondeds = atom.bonded().iterator();
        while((bonded=(Atom)bondeds.next())!=null){
            if (isConnectedToNonFrozen(bonded,depth-1)) return true;
        }

        return false;
    }

    private boolean inRadiusFromNonFrozen(Atom atom) {

        Atom nonFrozen;
        Iterator nonFrozensIter = nonFrozens.iterator();
        while ((nonFrozen = (Atom)nonFrozensIter.next())!=null){
            if (atom.distanceFrom(nonFrozen) <= radius) return true;
        }
        return false;
    }
}
