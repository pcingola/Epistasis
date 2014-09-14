package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.energy.EnergyCreator;
import meshi.energy.AbstractEnergy;
import meshi.util.KeyWords;
import meshi.util.CommandList;
import meshi.molecularElements.Protein;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 04/06/2009
 * Time: 11:47:03
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZSummaCreator extends EnergyCreator implements KeyWords {
    private static CooperativeZParameters parameters = null;
    protected      AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator;

    public CooperativeZSummaCreator(AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator ) {
            super(COOPERATIVE_Z_SUMMA_ENERGY );
            this.atomicPairwisePMFSummaCreator = atomicPairwisePMFSummaCreator;
        }

    public AbstractEnergy createEnergyTerm( Protein protein, DistanceMatrix distanceMatrix, CommandList commands ) {
        /* load parameters if this is the first time the function is called */
        if(( parameters == null ) | (weight() == 0))
        parameters = new CooperativeZParameters( commands );

        return new CooperativeZSumma( distanceMatrix, weight(), parameters,
                              (AtomicPairwisePMFSumma) atomicPairwisePMFSummaCreator.term());
    }


}
