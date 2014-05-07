package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.util.Gpr;

public class MsaAaAnalysis {

	public static final boolean debug = false;

	MultipleSequenceAlignmentSet msas;

	public static void main(String[] args) {
		String msaFileName = Gpr.HOME + "snpEff/db/multiz100way/head.fa";
		MultipleSequenceAlignmentSet msas = new MultipleSequenceAlignmentSet(msaFileName, 100);
		msas.load();

		MsaAaAnalysis maa = new MsaAaAnalysis(msas);
		maa.run();
	}

	public MsaAaAnalysis(MultipleSequenceAlignmentSet msas) {
		this.msas = msas;
	}

	public void run() {
		//		// Pre-calculate skip on all msas 
		//		Timer.showStdErr("Calculating AA probablities.");
		//		msas.getMsas().stream().forEach(MsaAaAnalysis::probAa);
		//		Timer.showStdErr("Done.");
	}

}
