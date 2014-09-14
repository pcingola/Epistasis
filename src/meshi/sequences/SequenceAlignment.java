package meshi.sequences;
import meshi.util.filters.*;
import meshi.util.*;
import meshi.util.string.*;
import meshi.sequences.aligner.*;
import java.util.*;

public class SequenceAlignment extends ArrayList<SequenceAlignmentColumn>{
    private Double score = null;
    public final StringList comments;
   public SequenceAlignment() {
	super();
	comments = new StringList();
    }

    public SequenceAlignment(ResidueAlignment residueAlignment) {
	this();
	for (Iterator resColumns = residueAlignment.iterator(); resColumns.hasNext();)
	    add(new SequenceAlignmentColumn((ResidueAlignmentColumn) resColumns.next()));
	for (String s:residueAlignment.comments)
	    comments.add(s);
    }
    public SequenceAlignment(Sequence sequence1, Sequence sequence2) {
	this();
	SequenceList sequenceList = new SequenceList();
	sequenceList.add(sequence1);
	sequenceList.add(sequence2);
	SequenceAlignment temp = new SequenceAlignment(sequenceList);
	for (SequenceAlignmentColumn sac:temp)
	    add(sac);
	for (String s:temp.comments)
	    comments.add(s);
    }
			
    public SequenceAlignment(SequenceList sequenceList) { 
	this();
	int numberOfSequences = sequenceList.size();
	boolean done = false;

	int length = (sequenceList.get(0)).size();
	for (int i = 0; i < numberOfSequences; i++) {
	    Sequence sequence = sequenceList.get(i);
	    if (sequence.size() != length) {
		System.out.println("A problem with sequences:");
		sequenceList.print();
		throw new RuntimeException("All sequences must have the same Length");
	    }
	    comments.add(sequence.comment());
	}
	    
	for (int iPosition = 0; iPosition < length; iPosition++) {
	    SequenceAlignmentColumn column = new SequenceAlignmentColumn(numberOfSequences);
	    for (int iSequence = 0; iSequence < numberOfSequences; iSequence++) {
		Sequence sequence = sequenceList.get(iSequence);
		column.add(iSequence, sequence.cell(iPosition));
	    }
	    add(column);
	}
    } 
    




    public String toString() {
	SequenceList sequenceList = new SequenceList(this);
	String out = "";
	for (Iterator sequences = sequenceList.iterator(); sequences.hasNext();)
	    out+=""+sequences.next();
	return out;
    }

    public boolean isExactMach() {
	for (Iterator columns = iterator(); columns.hasNext();)
	    if (! ((SequenceAlignmentColumn) columns.next()).isExactMach()) return false;
	return true;
    }
                   
    public boolean isExactMachWithGaps() {
	for (Iterator columns = iterator(); columns.hasNext();)
	    if (! ((SequenceAlignmentColumn) columns.next()).isExactMachWithGaps()) return false;
	return true;
    }

    public static SequenceAlignment identityAlignment(Sequence sequence1, Sequence sequence2) {
	DpMatrix matrix = new DpMatrix(sequence1, sequence2, new Identity(-0.2));
	SequenceAlignment out = matrix.backTrack();
	return out;
    }
	
    public double score() {
	if (score == null) throw new RuntimeException("Alignment without score");
	return score.doubleValue();
    }

    public void setScore(double score) {
	this.score = new Double(score);
    }

    public void print() {
	for (SequenceAlignmentColumn sac:this)
	    System.out.println(sac);
    }

    /**
     * Fetches a column with a specific cell number in a specific row. This column now is
     * a handle to the corresponding elements of the other proteins.
     **/
    public AlignmentColumn getColumn(int row, int number) { 
        for (Iterator columns = iterator(); columns.hasNext();) { 
            AlignmentColumn column = (AlignmentColumn) columns.next(); 
            if (column.cell(row).number == number) return column; 
        } 
        return null; 
    } 



}
