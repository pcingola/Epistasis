package meshi.applications.prediction.beautify.warpEnergy;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentColumn;
import meshi.util.CommandList;

import java.util.Iterator;


public class WarpEnergyCreator extends EnergyCreator {
    public final Protein template;
    public final double warpThreshold;

    public WarpEnergyCreator() {
        super(WARP_ENERGY);
        template = null;
        warpThreshold = -1;
    }

    public WarpEnergyCreator(Protein template, double warpThreshold) {
        super(WARP_ENERGY);
        this.template = template;
        this.warpThreshold = warpThreshold;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands) {
        term = new WarpEnergy(weight());
        if (template != null) {
            ResidueAlignment residueAlignment = new ResidueAlignment(protein.chains().get(0), protein.name(),
                    template.chains().get(0), template.name());
            //residueAlignment.print();
            for (Iterator columns = residueAlignment.iterator(); columns.hasNext();) {
                ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
                Residue templateResidue = (Residue) column.cell(1).obj;
                Residue proteinResidue = (Residue) column.cell(0).obj;
                if (!templateResidue.dummy()) {
                    if (templateResidue.ca().distanceFrom(proteinResidue.ca()) < warpThreshold) {
                        ((WarpEnergy) term).add(column);
                    }
                }
            }
        }
        return term;
    }
}
   
