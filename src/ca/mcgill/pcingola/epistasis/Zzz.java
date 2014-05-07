package ca.mcgill.pcingola.epistasis;

import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.NextProt;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;

public class Zzz extends SnpEff {

	public static void main(String[] args) {
		String genome = "testHg3771Chr1";
		String argsSnpEff[] = { "eff", "-v", "-nextProt", "-c", Gpr.HOME + "/snpEff/" + Config.DEFAULT_CONFIG_FILE, genome };

		Zzz zzz = new Zzz(argsSnpEff);

		zzz.setGenomeVer(genome);
		zzz.parseArgs(argsSnpEff);
		zzz.loadConfig();
		zzz.loadDb();

		CountByType count = new CountByType();
		for (Marker m : zzz.getConfig().getSnpEffectPredictor().getMarkers()) {
			if (m instanceof NextProt) {
				NextProt np = (NextProt) m;
				Gene g = (Gene) np.findParent(Gene.class);
				if (g != null && g.getId().equals("ENSG00000198125")) System.out.println(np.getChromosomeName() + "\t" + np.getStart() + "\t" + np.getEnd() + "\t" + np.getId() + "\t" + np.size());
				count.addScore(np.getId(), np.size());
			}
		}

		for (String type : count.keysSorted()) {
			double avg = 0;
			if (count.get(type) > 0) avg = count.getScore(type) / count.get(type);
			if (avg < 10) System.out.println(type + "\t" + count.get(type) + "\t" + count.getScore(type) + "\t" + avg);
		}
	}

	public Zzz(String[] args) {
		super(args);
	}
}
