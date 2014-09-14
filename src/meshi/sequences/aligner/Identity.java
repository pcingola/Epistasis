package meshi.sequences.aligner;
import meshi.sequences.*;

public class Identity implements CellScorer {
    public final double gapPenalty;
    
    public Identity(double gapPenalty) {
	this.gapPenalty = gapPenalty;
    }

    public double getScore(Cell cell) {
	DpMatrix matrix = cell.matrix;

	int row = cell.rowNumber;
	char rowChar = matrix.rowChar(row);
	int column = cell.colNumber;
	char columnChar = matrix.columnChar(column);
	double score;
	if ((rowChar == SequenceAlignmentCell.GAP_CHAR) |
	    (columnChar == SequenceAlignmentCell.GAP_CHAR) |
	    (rowChar == SequenceAlignmentCell.WILDCARD_CHAR) |
	    (columnChar == SequenceAlignmentCell.WILDCARD_CHAR))
	    score = 0;
	else score = ((rowChar == columnChar)?1:0);

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
