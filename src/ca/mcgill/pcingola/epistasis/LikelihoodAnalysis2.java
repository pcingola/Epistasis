package ca.mcgill.pcingola.epistasis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

/**
 * Logistic regression log-likelihood analysis of 2 VCF entries + phenotype data
 *
 * @author pcingola
 */
public class LikelihoodAnalysis2 extends LikelihoodAnalysis {

	HashMap<String, byte[]> gtByKey;

	public static void main(String[] args) {
		Timer.showStdErr("Start");

		boolean debug = false;

		LikelihoodAnalysis2 zzz = new LikelihoodAnalysis2(args);

		if (debug) {
			zzz.setDebug(debug);
			zzz.setLogLikInfoField(VCF_INFO_LOG_LIKELIHOOD);
		}

		zzz.run(debug);

		Timer.showStdErr("End");
	}

	public LikelihoodAnalysis2(String args[]) {
		super(args);
	}

	@Override
	protected List<VcfEntry> run(VcfFileIterator vcf, boolean createList) {
		gtByKey = new HashMap<String, byte[]>();

		// Process and populate list of VCF entries (single thread, used for debugging and test cases)
		Timer.show("Reading VCF file");
		int count = 1;
		for (VcfEntry ve : vcf) {
			String key = ve.toStr();
			byte gt[] = ve.getGenotypesScores();
			gtByKey.put(key, gt);
			Gpr.showMark(count++, 1000);
		}

		Timer.show("Done VCF file: " + gtByKey.size() + " entries");

		return new ArrayList<VcfEntry>();
	}
}
