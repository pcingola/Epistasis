package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.geometry.*;
import meshi.parameters.ResidueType;

public class CompositePropensityEnergy extends SimpleEnergyTerm implements CompositeTorsionsDefinitions {
    protected static double[] sumPerResidueType       = new double[ResidueType.values().length];
    protected static double[] sm2PerResidueType       = new double[ResidueType.values().length];
    protected static int[]    numberOfResiduesPerType = new int[ResidueType.values().length];
    public static double[] sumPerResidueType()       {return sumPerResidueType;}
    public static double[] sm2PerResidueType()       {return sm2PerResidueType;}
    public static int[]    numberOfResiduesPerType() {return numberOfResiduesPerType;}
    protected ResidueTorsionsList residueTorsionsList;
    public ResidueTorsionsList residueTorsionsList() {return residueTorsionsList;}

    public CompositePropensityEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			CompositePropensityParametersList cppl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), cppl, weight );
        this.residueTorsionsList = residueTorsionsList;
        for (int i = 0; i < ResidueType.values().length;i++)
            numberOfResiduesPerType[i] = 0;

        this.comment = comment;
		createElementsList( residueTorsionsList );
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
		CompositePropensityParameters cpp =
			(CompositePropensityParameters) parameters;
		
		return new CompositePropensityEnergyElement(
					resTorsions, cpp, weight );
    }
	
}
