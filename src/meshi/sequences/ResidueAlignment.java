package meshi.sequences;
import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.util.string.*;
import meshi.util.*;
import java.util.*;

public class ResidueAlignment extends ArrayList<ResidueAlignmentColumn> {
    private ResidueAlignmentColumn lastColumn;
    public final StringList comments;
     /**
     * Empty alignment.
     **/
    public ResidueAlignment() {
	super();
	lastColumn = null;
	comments = new StringList();
    }

    /**
     * A trevial alignment of two protein objects which are assumed to be two models of the same protein.
     **/
    public ResidueAlignment(Chain chain1, String name1, Chain chain2, String name2) {
	this();
	comments.add(name1);
	comments.add(name2);
	Sequence sequence1 = chain1.sequence();
	Sequence sequence2 = chain2.sequence();

	SequenceAlignment sa = SequenceAlignment.identityAlignment(sequence1,sequence2);
	for (Iterator columns = sa.iterator(); columns.hasNext();) {
	    AlignmentColumn column   = (AlignmentColumn) columns.next();
	    AlignmentCell   cell1    = column.cell(0);
	    AlignmentCell   cell2    = column.cell(1);
	    Residue         residue1 = (Residue) cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    Residue         residue2 = (Residue) cell2.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    if ((residue1 != null) & (residue2 != null)) {
		ResidueAlignmentColumn newColumn = new ResidueAlignmentColumn(residue1,residue2);
		add(newColumn);
	    }
	}
    }
		

    public ResidueAlignment(ResidueSequence  residueSequence0, ResidueSequence  residueSequence1) {
	this(SequenceAlignment.identityAlignment(residueSequence0,residueSequence1));
    }

    public ResidueAlignment(Chain chain0, Chain chain1, SequenceList sequenceList) {
	this();
	AlignmentColumn column;
	AlignmentCell   from;
	AlignmentCell   to;
	Residue                 residue, residue0, residue1;
	AlignmentCell cell0, cell1;

	Sequence          sequence0          = chain0.sequence();
	Sequence          sequence1          = chain1.sequence();
	Sequence          sequence0FromSL    = sequenceList.get(0);
	Sequence          sequence1FromSL    = sequenceList.get(1);
	SequenceAlignment sequence0Alignment = SequenceAlignment.identityAlignment(sequence0, sequence0FromSL);
	SequenceAlignment sequence1Alignment = SequenceAlignment.identityAlignment(sequence1, sequence1FromSL);

	for (Iterator columns = sequence0Alignment.iterator(); columns.hasNext();) {
	    column  = (AlignmentColumn) columns.next();
	    from    = column.cell(0);
	    to      = column.cell(1);
	    residue = (Residue)from.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    if (residue != null) to.addAttribute(residue);
	}
	for (Iterator columns = sequence1Alignment.iterator(); columns.hasNext();) {
	    column  = (SequenceAlignmentColumn) columns.next();
	    from    = column.cell(0);
	    to      = column.cell(1);
	    residue = (Residue)from.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    if (residue != null) to.addAttribute(residue);
	}
	
	SequenceAlignment sequenceAlignment = new SequenceAlignment(sequenceList);
	for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext();) {
	    column   = (SequenceAlignmentColumn) columns.next();
            cell0    = column.cell(0);
	    cell1    = column.cell(1);
	    residue0 = (Residue) cell0.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    residue1 = (Residue) cell1.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
	    if ((residue0 != null) & (residue1 != null)& (! SequenceAlignmentCell.wildcardOrGap((SequenceAlignmentColumn)column))) {
		add(new ResidueAlignmentColumn(residue0, residue1));
	    }
	}
    }



    public ResidueAlignment(SequenceAlignment sequenceAlignment) {
	this();
	int key = MeshiAttribute.RESIDUE_ATTRIBUTE;
	for (Iterator columns = sequenceAlignment.iterator(); columns.hasNext();) {
	    SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
	    SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
	    Residue residue0 = (Residue) cell0.getAttribute(key);
	    Residue residue1 = (Residue) cell1.getAttribute(key);
	    if ((residue0!= null) & (residue1 != null) & (! SequenceAlignmentCell.wildcardOrGap(column))) {
		ResidueAlignmentCell residueCell0 = new ResidueAlignmentCell(residue0);
		ResidueAlignmentCell residueCell1 = new ResidueAlignmentCell(residue1);
		ResidueAlignmentColumn residueColumn = new ResidueAlignmentColumn(residueCell0,residueCell1);
		add(residueColumn);
	    }
	}
    }
    


    public boolean add(ResidueAlignmentColumn column) {
	if (lastColumn == null) {
	    lastColumn = column; 
	    return super.add(column);
	}
	if (!column.cell0().gap())
		if (column.residue0().number() <= lastColumn.residue0().number()) return false;
	if (!column.cell1().gap())
		if (column.residue1().number() <= lastColumn.residue1().number()) return false;
	if ((!column.cell0().gap()) & 
	    (!column.cell1().gap())) lastColumn = column;
	return super.add(column);
    }

   public ResidueAlignment insert(ResidueAlignmentColumn column) {
	ResidueAlignment out = new ResidueAlignment();
	for (String c:comments)
	    out.comments.add(c);
	boolean added = false;
	for (Iterator columns = iterator(); columns.hasNext();) {
	    ResidueAlignmentColumn current = (ResidueAlignmentColumn) columns.next();
	    if  (! added) {
		if ((current.residue0().number() < column.residue0().number()) &
		    (current.residue1().number() < column.residue1().number())){
		    if (! out.add(current)) throw new RuntimeException("weird situation 1");
		}
		else if ((current.residue0().number() > column.residue0().number()) &
			 (current.residue1().number() > column.residue1().number())) {
		    if (! out.add(column)) throw new RuntimeException("weird situation 2");
		    if (! out.add(current)) throw new RuntimeException("weird situation 3");
		    added = true;
		}
		else {
			return (null);
		}
	    }
	    else {
		if ((current.residue0().number() <= column.residue0().number()) |
		    (current.residue1().number() <= column.residue1().number())){
		    throw new RuntimeException("weird situation 4");
		}
		if (! out.add(current)) throw new RuntimeException("weird situation 5");
	    }
	}
	if (! added) if (! out.add(column)) throw new RuntimeException("weird situation 6");
	return out;
    }
		
		    
	    
	    
    public AtomList getCaList(int row) {
	return (new AtomAlignment(this,new CaFilter())).atomList(row);
    }

    public String toString() {
	return (new SequenceAlignment(this)).toString();
    }
       

    public ResidueAlignment filter(Filter filter) {
	ResidueAlignment out = new ResidueAlignment();
	for(ResidueAlignmentColumn rc:this)
	    if (filter.accept(rc)) out.add(rc);
	return out;
    }
}
