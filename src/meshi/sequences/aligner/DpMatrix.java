package meshi.sequences.aligner;
import meshi.sequences.*;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 29/07/2005
 * Time: 15:47:31
 * To change this template use File | Settings | File Templates.
 */
public class DpMatrix {
    Sequence sequence1, sequence2; // One dimentional sequenceAlignments
    public final CellScorer cellScorer;
    private Cell [][] matrix;
    


    public DpMatrix(Sequence sequence1, Sequence sequence2, CellScorer scrr){
	this.sequence1 = sequence1;
	this.sequence2 = sequence2;
	cellScorer = scrr;
	matrix = new Cell[sequence1.size()+1] [sequence2.size()+1];
	for (int row=0; row < sequence1.size()+1;row++ ){
	    for (int col=0; col < sequence2.size()+1;col++ ){
		setCell(row, col, new Cell(row , col, this));
	    }
	}
    }
    public void setCell(int row, int column, Cell cell){
	matrix[row][column] = cell;
    }
    public Cell getCell(int rowNumber, int colNumber) {
		// TODO Auto-generated method stub
	return matrix[rowNumber][colNumber];
    }

    public SequenceAlignment backTrack() {
	SequenceAlignment inversAlignment = new SequenceAlignment();

	Cell cell = matrix[sequence1.size()][sequence2.size()];
	double score = cell.score;
	while (cell != matrix[0][0]) {
	    Cell back = cell.getBack();
	    Cell up = cell.upCell;
	    Cell left = cell.leftCell;
	    Cell upLeft = cell.upLeftCell;
	    int rowNumber = cell.rowNumber;
	    int colNumber = cell.colNumber;
	    SequenceAlignmentColumn column;
	    //System.out.println("xxxxxx0 "+cell);
	    if (back == up) {
		column = new SequenceAlignmentColumn((SequenceAlignmentCell) sequence1.cell(rowNumber-1), 
						     new SequenceAlignmentCell());
		cell = up;
		//System.out.println("xxxxxx01 "+column);
	    }
	    else if (back == left){
		column = new SequenceAlignmentColumn(new SequenceAlignmentCell(), 
						     (SequenceAlignmentCell) sequence2.cell(colNumber-1));
		cell = left;
		//System.out.println("xxxxxx02 "+column);
	    }
	    else {
		column = new SequenceAlignmentColumn((SequenceAlignmentCell) sequence1.cell(rowNumber-1), 
						     (SequenceAlignmentCell) sequence2.cell(colNumber-1));
		cell = upLeft;
		//System.out.println("xxxxxx03 "+column);
	    }
	    //System.out.println("xxxxxx "+column);
	    inversAlignment.add(column);
	}
	
	SequenceAlignment out = new SequenceAlignment();
	for (int i = inversAlignment.size() - 1; i >= 0; i--)
	    out.add(inversAlignment.get(i));
	out.comments.add(sequence1.comment());
	out.comments.add(sequence2.comment());
	out.setScore(score);
	return out;
    }
		
    public char rowChar(int index) {
	if (index == 0) return SequenceAlignmentCell.GAP_CHAR;
	return ((SequenceAlignmentCell) sequence1.get(index-1).cell(0)).getChar();
    }

    public char columnChar(int index) {
	if (index == 0) return SequenceAlignmentCell.GAP_CHAR;
	return ((SequenceAlignmentCell) sequence2.get(index-1).cell(0)).getChar();
    }
	
}


   

	
