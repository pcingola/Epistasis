package meshi.molecularElements;
import meshi.molecularElements.atoms.*; 
import meshi.parameters.*;
public interface ResidueCreator {
    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode);
}
