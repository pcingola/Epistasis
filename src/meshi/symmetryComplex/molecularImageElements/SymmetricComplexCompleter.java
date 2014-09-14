package meshi.symmetryComplex.molecularImageElements;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.symmetryComplex.transformations.Transformation;
import meshi.symmetryComplex.topologyMap.Topology;
import meshi.symmetryComplex.proteinsGJ.CoordinatesAssigner;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.geometry.Coordinates;

import java.io.IOException;


public class SymmetricComplexCompleter extends SymmetricComplex implements KeyWords {
    //CommandList commands;
    Topology topology;
    public SymmetricComplexCompleter (AtomList initialAtoms, ResidueCreator creator,
                                      Transformation[] transformations,
                                      CommandList commands, Topology topology)throws IOException{
        super(initialAtoms, creator, transformations);
        this.topology = topology;
        assignCoordinates(commands, topology);
        test();
    }

    private void assignCoordinates(CommandList commands, Topology topology) throws IOException {
    int numberOfResiduesWithoutCoordinates;
    int report = 0;
    CoordinatesAssigner coordinatesAssigner;
    int maxNumberOfClashes;
    double clashDistance;
    int nTries;

    clashDistance = commands.firstWord(CLASH_DISTANCE).secondWordDouble();
    System.out.println("clashDistance = "+clashDistance);
    maxNumberOfClashes = commands.firstWord(MAX_CLASHES).secondWordInt();
    System.out.println("maxNumberOfClashes = "+maxNumberOfClashes);
    nTries = commands.firstWord(N_TRIES).secondWordInt();
    System.out.println("nTries = "+nTries);

    ResidueList caResidues = getSource().filter(new ResidueList.NonDummyFilter());// nonDummyResidues();
    if (getNumberOfResiduesWithCoordinates(caResidues) < 1)
                      throw new RuntimeException("Can't to assign residues to empty complex!");

    coordinatesAssigner = new CoordinatesAssigner(this, commands, caResidues, clashDistance, maxNumberOfClashes, nTries, topology);

    System.out.println("Starting to add missing residues to "+this);
    while ((numberOfResiduesWithoutCoordinates = getNumberOfResiduesWithoutCoordinates(caResidues)) > 0) {
        if (report % 100 == 0) {
        System.out.println("\n"+numberOfResiduesWithoutCoordinates+" residues do not have coordinates.\n");
        }
        report++;
        coordinatesAssigner.step();
        if (report > caResidues.size()*100)    {
            throw new RuntimeException("Cann't find coordinates for all residues." +
                    "\n"+numberOfResiduesWithoutCoordinates+" residues do not have coordinates. \n");
        }
    }
    }

    private int getNumberOfResiduesWithoutCoordinates(ResidueList residues) {
    int count = 0;
        for (Residue residue : residues)
                  if ((!residue.dummy()) && (residue.ca().nowhere()))
                      count++;
    return count;
    }

    private int getNumberOfResiduesWithCoordinates(ResidueList residues) {
    int count = 0;
        for (Residue residue : residues)
                if ((!residue.dummy()) && (!residue.ca().nowhere())) count++;
    return count;
    }

    private void  test(){
            for (Atom atom :atoms){
                 if (atom.nowhere())
                    throw new RuntimeException("Nowhere atom is found in the normal part of SymmetryComplex in chain "+atom.chain()+".\n"+atom+atom.residue());
                if (atom.core.status().image() && ((atom.core.x()== Coordinates.NOWHERE_CONST ) || (atom.core.y()== Coordinates.NOWHERE_CONST ) ||
                                (atom.core.z()== Coordinates.NOWHERE_CONST )) )
                    throw new RuntimeException("Nowhere atom in image part of SymmetryComplex in chain "+atom.chain()+".\n"+atom+atom.residue());
            }
    }


}
