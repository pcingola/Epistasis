package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.*;
import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.util.*;

public class AtomicPairwisePMFSummaCreator extends
		EnergyCreator implements KeyWords {

	private static AtomicPairwisePMFSummaParameters parameters = null;

	public AtomicPairwisePMFSummaCreator() {
		super( ATOMIC_PAIRWISE_PMF_SUMMA_ENERGY );
	}

	public AtomicPairwisePMFSummaCreator( double weight ) {
		super( weight );
	}

	public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix,
			CommandList commands ) {
		/* load parameters if this is the first time the function is called */
		if( parameters == null )
			parameters = new AtomicPairwisePMFSummaParameters( commands );

		return term = new AtomicPairwisePMFSumma( distanceMatrix, weight(), parameters );
	}

}
