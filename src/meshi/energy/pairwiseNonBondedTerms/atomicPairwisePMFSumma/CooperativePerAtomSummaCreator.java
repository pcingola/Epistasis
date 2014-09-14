package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 01/06/2009
 * Time: 13:15:40
 * To change this template use File | Settings | File Templates.
 */
public class CooperativePerAtomSummaCreator  extends EnergyCreator implements KeyWords {
    private static CooperativePerAtomSummaParameters parameters = null;
    protected      AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator;
    public CooperativePerAtomSummaCreator(AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator ) {
            super(COOPERATIVE_PERATOM_SUMMA_ENERGY );
            this.atomicPairwisePMFSummaCreator = atomicPairwisePMFSummaCreator;
        }

    public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
        /* load parameters if this is the first time the function is called */
        if(( parameters == null ) | (weight() == 0))
        parameters = new CooperativePerAtomSummaParameters( commands );

        return new CooperativePerAtomSumma( distanceMatrix, weight(), parameters,
                              (AtomicPairwisePMFSumma) atomicPairwisePMFSummaCreator.term());
    }


}
