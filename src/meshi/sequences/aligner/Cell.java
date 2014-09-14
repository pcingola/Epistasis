package meshi.sequences.aligner;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 29/07/2005
 * Time: 15:39:55
 * To change this template use File | Settings | File Templates.
 */
public class Cell {
    public final int rowNumber;
    public final int colNumber;
    public final Cell upCell, upLeftCell, leftCell;//the next cell will be either the cell in the right of this cell or the cell in bottom of this one
    protected Cell back;
    public final double score;
    public final DpMatrix matrix;

    Cell(int i, int j, DpMatrix mat){
    	rowNumber = i;
    	colNumber = j;
    	matrix = mat;
    	//the initation of the neighbour cells
    	if(rowNumber == 0) {
	    upCell = null;
    	}
    	else{
	    upCell = matrix.getCell(rowNumber - 1, colNumber);
    	}
 
   	if(rowNumber == 0 || colNumber == 0) {
	    upLeftCell = null;
    	}
    	else{
	    upLeftCell = matrix.getCell(rowNumber - 1, colNumber - 1);
    	}

    	if(colNumber == 0) {
	    leftCell = null;
    	}
    	else{
	    leftCell = matrix.getCell(rowNumber , colNumber - 1 );
    	}
    	//thes core for the cell
    	
	score = matrix.cellScorer.getScore(this);
    }   		
        
    public void setBack(Cell back) {
	this.back = back;
    }
    public Cell getBack() {
	return back;
    }
	
    public String toString() {
	return "cell "+rowNumber+" "+colNumber+" "+score;
    }
	
}
