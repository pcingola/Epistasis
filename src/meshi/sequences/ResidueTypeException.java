package meshi.sequences;
public class ResidueTypeException extends RuntimeException {
    public ResidueTypeException(Character weirdChar, Sequence sequence) {
	super("Weird character "+weirdChar+" in \n"+sequence);
    }
    public ResidueTypeException(Character weirdChar) {
	super("Weird character "+weirdChar);
    }
}
