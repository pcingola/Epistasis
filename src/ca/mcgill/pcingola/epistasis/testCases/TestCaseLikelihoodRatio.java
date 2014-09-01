package ca.mcgill.pcingola.epistasis.testCases;

import junit.framework.TestCase;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.Epistasis;

/**
 * Test case likelihood ratio
 *
 * @author pcingola
 */
public class TestCaseLikelihoodRatio extends TestCase {

	public static boolean debug = true;
	public static boolean verbose = true;
	public static final double EPSILON = 1E-6;
	public static final double MAX_LL_ERROR = 0.2;

	public void test_01() {

		String treeFile = "test/hg19.100way.nh";
		String multAlignFile = "test/msas.best.head.fa.gz";
		String idMapFile = "test/idMap_ensemblId_refseq_pdbId.best.txt";
		String qMatrixFile = "test/Qhat.txt";
		String aaFreqsFile = "test/aa.frequencies.txt";
		String q2MatrixFile = "test/Qhat2.txt";
		String aaFreqsContactFile = "test/aa.contact.frequencies.txt";

		String tests = "test/likelihood.test.test_01.txt";

		Epistasis ep = new Epistasis();
		ep.setTreeFile(treeFile);
		ep.setMultAlignFile(multAlignFile);
		ep.setIdMapFile(idMapFile);
		ep.setqMatrixFile(qMatrixFile);
		ep.setAaFreqsFile(aaFreqsFile);
		ep.setQ2MatrixFile(q2MatrixFile);
		ep.setAaFreqsContactFile(aaFreqsContactFile);

		ep.load();
		ep.precalcExps();

		Timer timer = new Timer();
		int lineNum = 1;
		for (String line : Gpr.readFile(tests).split("\n")) {
			// Parse line
			if (verbose) System.out.println(lineNum + ":\t" + line);
			String field[] = line.split("\t");
			String msaId1 = field[0];
			int msaIdx1 = Gpr.parseIntSafe(field[1]);
			String msaId2 = field[2];
			int msaIdx2 = Gpr.parseIntSafe(field[3]);
			double ll = Gpr.parseDoubleSafe(field[4]);

			// Calculate
			String res = ep.likelihoodRatio(msaId1, msaIdx1, msaId2, msaIdx2, false);

			// Parse result
			String fieldsRes[] = res.split("\t");
			double llres = Gpr.parseDoubleSafe(fieldsRes[2]);

			if (debug) {
				System.out.println("\t" + res);
				System.out.println("\t\tResult: " + llres + "\t" + fieldsRes[2] + "\n");
			}

			assertEquals(ll, llres, MAX_LL_ERROR);
			lineNum++;
		}

		timer.end();
		System.out.println("Elapsed: " + timer);
	}
}
