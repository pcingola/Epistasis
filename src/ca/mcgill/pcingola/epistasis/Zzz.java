package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

public class Zzz {

	public static void main(String[] args) {
		String dir = Gpr.HOME + "/snpEff/epistasis";
		String vcfFileName = dir + "/pos.vcf";
		VcfFileIterator vcf = new VcfFileIterator(vcfFileName);

		StringBuilder sb = new StringBuilder();
		int entryNum = 1;
		byte[] scoresPrev = null;
		String label = "", labelPrev = "";
		for (VcfEntry ve : vcf) {
			byte[] scores = ve.getGenotypesScores();

			label = ve.toStr();
			sb.append(show(label, scores));

			// Every two lines
			if (entryNum % 2 == 0) {
				byte[] epiMin = new byte[scores.length];
				byte[] epi = new byte[scores.length];
				for (int j = 0; j < scores.length; j++) {
					if ((scores[j] >= 0) && (scoresPrev[j] >= 0)) {
						epiMin[j] = (byte) (Math.min(scores[j], scoresPrev[j]));
						epi[j] = (byte) (scores[j] * scoresPrev[j]);
					}
				}

				sb.append(show("Epi", epi));
				sb.append(show("Epi min", epiMin));

				sb.append(show("Var " + labelPrev, variant(scoresPrev)));
				sb.append(show("Var " + label, variant(scores)));
				sb.append(show("Var Epi", variant(epi)));
				sb.append("\n");
				System.out.println("\n");
			}

			entryNum++;
			scoresPrev = scores;
			labelPrev = label;
		}

		Gpr.toFile(dir + "/pos.out", sb);
	}

	static byte[] variant(byte[] gt) {
		byte[] var = new byte[gt.length];
		for (int j = 0; j < gt.length; j++) {
			var[j] = (gt[j] > 0 ? 1 : gt[j]);
		}

		return var;
	}

	static String show(String label, byte[] gt) {
		int ac = 0, count = 0;

		StringBuilder gtSb = new StringBuilder();
		for (int j = 0; j < gt.length; j++) {
			if (gt[j] >= 0) {
				gtSb.append(gt[j]);
				ac += gt[j];
				count++;
			} else gtSb.append(' ');
		}

		double af = (100.0 * ac) / (2.0 * count);
		String strShort = String.format("%-30s\tac: %8d\taf: %.2f\tlen: %8d", label, ac, af, gt.length);
		String str = String.format("%-30s\tac: %8d\taf: %.2f\tlen: %8d\t%s", label, ac, af, gt.length, gtSb);
		System.out.println(strShort);
		return str + "\n";
	}
}
