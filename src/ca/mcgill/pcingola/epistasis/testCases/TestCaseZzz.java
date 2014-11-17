package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.GprSeq;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.Epistasis;
import ca.mcgill.pcingola.epistasis.gwas.InteractionLikelihood;
import ca.mcgill.pcingola.epistasis.phylotree.LikelihoodTreeAa;

/**
 * Test cases for logistic regression
 *
 * @author pcingola
 */
public class TestCaseZzz extends TestCase {

	public static boolean debug = false;
	public static boolean verbose = true || debug;
	public static final double MAX_LL_ERROR = 0.2;

	byte[] sequence(char aa1, char aa2, int n, int len) {
		byte s[] = new byte[len];
		byte a1 = GprSeq.aa2Code(aa1);
		byte a2 = GprSeq.aa2Code(aa2);

		for (int i = 0; i < n; i++)
			s[i] = a1;

		for (int i = n; i < len; i++)
			s[i] = a2;

		return s;
	}

	//	public void test_zzz() {
	//		String treeFile = "test/hg19.100way.nh";
	//		String multAlignFile = "test/msas.best.head.fa.gz";
	//		String idMapFile = "test/idMap_ensemblId_refseq_pdbId.best.txt";
	//		String qMatrixFile = "test/Qhat.txt";
	//		String aaFreqsFile = "test/aa.frequencies.txt";
	//		String q2MatrixFile = "test/Qhat2.txt";
	//		String aaFreqsContactFile = "test/aa.contact.frequencies.txt";
	//
	//		Epistasis ep = new Epistasis();
	//		ep.setTreeFile(treeFile);
	//		ep.setMultAlignFile(multAlignFile);
	//		ep.setIdMapFile(idMapFile);
	//		ep.setqMatrixFile(qMatrixFile);
	//		ep.setAaFreqsFile(aaFreqsFile);
	//		ep.setQ2MatrixFile(q2MatrixFile);
	//		ep.setAaFreqsContactFile(aaFreqsContactFile);
	//
	//		ep.load();
	//		InteractionLikelihood il = ep.newInteractionLikelihood();
	//		il.precalcExps();
	//
	//		LikelihoodTreeAa treeNull = il.getTreeNull();
	//		LikelihoodTreeAa treeAlt = il.getTreeAlt();
	//
	//		// Sequences are 100 bases long
	//		int len = 100;
	//		double llmax = Double.NEGATIVE_INFINITY;
	//
	//		int n = 50;
	//		llmax = 0;
	//		String best = "";
	//		Timer.showStdErr("n: " + n);
	//		for (char aa1s1 : GprSeq.AMINO_ACIDS) {
	//			for (char aa2s1 : GprSeq.AMINO_ACIDS) {
	//				double llmaxlocal = 0;
	//				Timer.showStdErr(aa1s1 + "\t" + aa2s1);
	//				for (char aa1s2 : GprSeq.AMINO_ACIDS) {
	//					for (char aa2s2 : GprSeq.AMINO_ACIDS) {
	//						byte seq1[] = sequence(aa1s1, aa2s1, n, len);
	//						byte seq2[] = sequence(aa1s2, aa2s2, n, len);
	//
	//						double llnull = il.likelihoodNullModel(treeNull, seq1, seq2);
	//						double llalt = il.likelihoodAltModel(treeAlt, seq1, seq2);
	//						double logLikRatio = -2.0 * (Math.log(llnull) - Math.log(llalt));
	//
	//						if (logLikRatio > llmax) {
	//							best = "GLOBAL\tn: " + n + "\t" + aa1s1 + "|" + aa2s1 + "\t" + aa1s2 + "|" + aa2s2 + "\tLL: " + logLikRatio + "\tLL_null: " + llnull + "\tLL_alt: " + llalt + "\t" + GprSeq.code2aa(seq1) + "\t" + GprSeq.code2aa(seq2);
	//							System.out.println(best);
	//							llmax = logLikRatio;
	//						}
	//
	//						//						if (logLikRatio > llmaxlocal) {
	//						//							System.out.println("LOCAL \tn: " + n + "\t" + aa1s1 + "|" + aa2s1 + "\t" + aa1s2 + "|" + aa2s2 + "\tLL: " + logLikRatio + "\tLL_null: " + llnull + "\tLL_alt: " + llalt + "\t" + GprSeq.code2aa(seq1) + "\t" + GprSeq.code2aa(seq2));
	//						//							llmaxlocal = logLikRatio;
	//						//						}
	//					}
	//				}
	//			}
	//		}
	//
	//		Timer.showStdErr(best);
	//
	//	}

	//D|C	V|C
	public void test_zzz() {
		String treeFile = "test/hg19.100way.nh";
		String multAlignFile = "test/msas.best.head.fa.gz";
		String idMapFile = "test/idMap_ensemblId_refseq_pdbId.best.txt";
		String qMatrixFile = "test/Qhat.txt";
		String aaFreqsFile = "test/aa.frequencies.txt";
		String q2MatrixFile = "test/Qhat2.txt";
		String aaFreqsContactFile = "test/aa.contact.frequencies.txt";

		Epistasis ep = new Epistasis();
		ep.setTreeFile(treeFile);
		ep.setMultAlignFile(multAlignFile);
		ep.setIdMapFile(idMapFile);
		ep.setqMatrixFile(qMatrixFile);
		ep.setAaFreqsFile(aaFreqsFile);
		ep.setQ2MatrixFile(q2MatrixFile);
		ep.setAaFreqsContactFile(aaFreqsContactFile);

		ep.load();
		InteractionLikelihood il = ep.newInteractionLikelihood();
		il.precalcExps();

		LikelihoodTreeAa treeNull = il.getTreeNull();
		LikelihoodTreeAa treeAlt = il.getTreeAlt();

		// Sequences are 100 bases long
		String best = "";
		int len = 100;
		double llmax = Double.NEGATIVE_INFINITY;
		//		char aa1s1 = 'D', aa2s1 = 'C';
		//		char aa1s2 = 'V', aa2s2 = 'C';
		char aa1s1 = 'A', aa2s1 = 'R';
		char aa1s2 = 'V', aa2s2 = 'D';

		for (int n = 0; n < len; n++) {
			byte seq1[] = sequence(aa1s1, aa2s1, n, len);
			byte seq2[] = sequence(aa1s2, aa2s2, n, len);

			double llnull = il.likelihoodNullModel(treeNull, seq1, seq2);
			double llalt = il.likelihoodAltModel(treeAlt, seq1, seq2);
			double logLikRatio = -2.0 * (Math.log(llnull) - Math.log(llalt));

			if (logLikRatio > llmax) {
				best = "GLOBAL\tn: " + n + "\t" + aa1s1 + "|" + aa2s1 + "\t" + aa1s2 + "|" + aa2s2 + "\tLL: " + logLikRatio + "\tLL_null: " + llnull + "\tLL_alt: " + llalt + "\t" + GprSeq.code2aa(seq1) + "\t" + GprSeq.code2aa(seq2);
				System.out.println(best);
				llmax = logLikRatio;
			}

		}

		Timer.showStdErr(best);

	}

}
