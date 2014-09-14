package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.*;
import meshi.molecularElements.*;
import meshi.geometry.*;
import meshi.util.*;

public class CooperativeAtomicPairwisePMFSummaCreator extends EnergyCreator implements KeyWords {

        private static CooperativeAtomicPairwisePMFSummaParameters parameters = null;
        protected      AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator;
        public CooperativeAtomicPairwisePMFSummaCreator(AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator ) {
            super(COOPERATIVE_ATOMIC_PAIRWISE_PMF_SUMMA_ENERGY );
            this.atomicPairwisePMFSummaCreator = atomicPairwisePMFSummaCreator;
        }

	public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
	    /* load parameters if this is the first time the function is called */
	    if(( parameters == null ) | (weight() == 0))
		parameters = new CooperativeAtomicPairwisePMFSummaParameters( commands );
	    
	    return new CooperativeAtomicPairwisePMFSumma( distanceMatrix, weight(), parameters, 
							  (AtomicPairwisePMFSumma) atomicPairwisePMFSummaCreator.term());
	}

}