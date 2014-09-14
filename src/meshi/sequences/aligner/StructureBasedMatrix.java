package meshi.sequences.aligner;
import meshi.sequences.*;

public class StructureBasedMatrix implements CellScorer {
    public final double gapPenalty;
    public final double[][] scores;
    
    public StructureBasedMatrix(double gapPenalty, double[][] scores) {
	this.gapPenalty = gapPenalty;
	this.scores = scores;
    }

    public double getScore(Cell cell) {
	DpMatrix matrix = cell.matrix;

	int row = cell.rowNumber;
	int column = cell.colNumber;
	double score = scores[row][column];

	Cell up = cell.upCell;
	Cell left = cell.leftCell;
	Cell upLeft = cell.upLeftCell;

	if ((up == null) & (left == null)) {
	    cell.setBack(null);
	    return 0;
	}	
	if (up == null) {
	    cell.setBack(left);
	    return 0; 
	}
	if (left == null) {
	    cell.setBack(up);
	    return 0;
	}

	double scoreUp = (matrix.sequence2.size() == column)?up.score:gapPenalty + up.score;
	double scoreLeft = (matrix.sequence1.size() == row)?left.score:gapPenalty + left.score;
	double scoreUpLeft = upLeft.score +score;
	
	double max = Math.max(scoreUp, scoreLeft);
	max = Math.max(max, scoreUpLeft);

	if (scoreUp == max) cell.setBack(up);
	else if (scoreLeft == max) cell.setBack(left);
	else if (scoreUpLeft == max) cell.setBack(upLeft);

	return max;
    }
}
