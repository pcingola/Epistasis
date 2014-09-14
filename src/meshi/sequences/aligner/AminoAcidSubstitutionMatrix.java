package meshi.sequences.aligner;

import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 30/07/2005
 * Time: 10:25:13
 * To change this template use File | Settings | File Templates.
 */
//supports matrices in the format which each line includes the aa it represents and the matrix is full not only diaganol half
public class AminoAcidSubstitutionMatrix implements SubstitutionMatrix {
	Integer[][] scoringMatrix = new Integer[20][20];
	String amino_a = new String(); //will be our help to know the indices in the matrix

	public AminoAcidSubstitutionMatrix(String mat_name) {
		String str = new String();
		FileToString fts = new FileToString(mat_name);
		StringTokenizer aaTknzr = new StringTokenizer(fts.getText().trim(), "ARNDCQEGHILKMFPSTWYV", true); //we want the AA as delimiters
		// str=  tknzr.nextToken();  //the first token will be the first letter in the matrix
		//  amino_a=str;
		for (int i = 0; i < 20; i++) {
			str = aaTknzr.nextToken(); // will get the AA
			amino_a = amino_a.concat(str);
			str = aaTknzr.nextToken(); //will get the line of numbers after the AA
			//System.out.println(i+"0)"+str);
			// System.out.println(amino_a);
			StringTokenizer scoresTknzr = new StringTokenizer(str.trim()); //making a tokenizer which will go htrew the scores
			for (int j = 0; j < 20; j++) {
				scoringMatrix[i][j] = Integer.valueOf(scoresTknzr.nextToken());
			}

		}

	}

	@Override
	public int getScore(char sourceAA, char imageAA) {
		int i = amino_a.indexOf(sourceAA);
		int j = amino_a.indexOf(imageAA);
		return scoringMatrix[i][j].intValue();
	}

	public void print() {
		for (int i = 0; i < 20; i++) {//change to global
			System.out.println("");
			System.out.print(i + " " + amino_a.charAt(i));//change to substring
			for (int j = 0; j < 20; j++) {
				System.out.print("  " + scoringMatrix[i][j]);

			}

		}
	}

}
