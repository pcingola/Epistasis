package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain;

import meshi.energy.*;
import meshi.parameters.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.geometry.DistanceMatrix;

import java.util.*;

/** A Ramachandran plot and sidechain torsions optimization energy term.
 * Statistical analysis of residue backbone and sidechain torsions in a
 * large database of residue observations has been smoothed using
 * polynomial spline interpolation.
 * For a given residue the energy value approximates the percentage
 * of finding its current backbone and sidechain torsion angles.
 * 
 * @author El-ad David Amir
 *
 */
public class RamachandranSidechainEnergy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {
    
    protected static double[] sumPerResidueType       = new double[ResidueType.values().length];
    protected static double[] sm2PerResidueType       = new double[ResidueType.values().length];
    protected static int[]    numberOfResiduesPerType = new int[ResidueType.values().length];
    public static double[] sumPerResidueType()       {return sumPerResidueType;}
    public static double[] sm2PerResidueType()       {return sm2PerResidueType;}
    public static int[]    numberOfResiduesPerType() {return numberOfResiduesPerType;}
    public static double[] residueEnergies;
    protected ResidueTorsionsList residueTorsionsList;
    public ResidueTorsionsList residueTorsionsList() {return residueTorsionsList;}

    public RamachandranSidechainEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			RamachandranSidechainParametersList rspl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), rspl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
		for (int i = 0; i < ResidueType.values().length;i++)
		    numberOfResiduesPerType[i] = 0;
		for (ResidueTorsions rt:residueTorsionsList) 
		    numberOfResiduesPerType[rt.getResidueType().ordinal()]++;
		this.residueTorsionsList = residueTorsionsList;
	}

    public double evaluate() {
	for (int i = 0; i < ResidueType.values().length;i++) {
	    sumPerResidueType[i] = 0;
	    sm2PerResidueType[i] = 0;
	}
    for (ResidueTorsions r:residueTorsionsList )
              r.resetEnergy();

    return super.evaluate();
    }
	
    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		RamachandranSidechainParameters rsp =
			(RamachandranSidechainParameters) parameters;
		
		switch( NUM_SIDECHAIN_TORSIONS[resTorsions.getResidueType().ordinal()] ) {
		case 0:
			return new RamachandranSidechainEnergyElementChi0(
					resTorsions, rsp, weight );
		case 1:
			return new RamachandranSidechainEnergyElementChi1(
					resTorsions, rsp, weight );
		case 2:
			return new RamachandranSidechainEnergyElementChi2(
					resTorsions, rsp, weight );
		case 3:
			return new RamachandranSidechainEnergyElementChi3(
					resTorsions, rsp, weight );
		case 4:
			return new RamachandranSidechainEnergyElementChi4(
					resTorsions, rsp, weight );
		default:
			throw new RuntimeException( "unidentified residue type" );
		}
	}
	
}
