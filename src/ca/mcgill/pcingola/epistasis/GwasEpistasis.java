package ca.mcgill.pcingola.epistasis;

import java.util.HashMap;

import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.SnpEff;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * Perform GWAS using epistasis data
 *
 * @author pcingola
 */
public class GwasEpistasis extends SnpEff {

	public static final double LOG_LIKELIHOOD_RATIO_THRESHOLD = -1000.0;

	int countOk, countErr;
	String genesLikeFile;
	String vcfFile;
	String genomeVer, pdbDir;
	HashMap<String, Transcript> trancriptById;

	public GwasEpistasis(String configFile, String genomeVer, String genesLikeFile, String vcfFile) {
		this.configFile = configFile;
		this.genomeVer = genomeVer;
		this.genesLikeFile = genesLikeFile;
		this.vcfFile = vcfFile;
	}

	public int getCountErr() {
		return countErr;
	}

	public int getCountOk() {
		return countOk;
	}

	/**
	 * Load all data
	 */
	public void initialize() {
		// Initialize SnpEff
		String argsSnpEff[] = { "eff", "-v", "-c", configFile, genomeVer };
		args = argsSnpEff;
		setGenomeVer(genomeVer);
		parseArgs(argsSnpEff);
		loadConfig();

		// Load SnpEff database
		if (genomeVer != null) loadDb();

		// Initialize trancriptById
		trancriptById = new HashMap<>();
		for (Gene g : config.getSnpEffectPredictor().getGenome().getGenes())
			for (Transcript tr : g) {
				String id = tr.getId();
				if (id.indexOf('.') > 0) {
					// When using RefSeq transcripts, we don't store sub-version number
					id = id.substring(0, id.indexOf('.'));
				}
				trancriptById.put(id, tr);
			}

	}

	/**
	 * Parse MSA_ID from gene likelihood file
	 *
	 * Format: trId _ chr : start - end [ aaIdx ]
	 *
	 * E.g.  : 	NM_004635_3:50679137-50679201[0]
	 * 			NM_006936_21:46226865-46226954[22]
	 * 
	 * @return A marker that contains the codon referenced by the ID
	 *			Note: The marker has ALL bases in the codon. 
	 *				  For instance, is the codon is split between two exons, the 
	 *				  marker will contain the intron

	 *
	 */
	protected Marker parseMsaId(String id, char aaExpected) {
		String idRep = id.replace(':', '_').replace('-', '_').replace('[', '_').replace(']', '_');
		String f[] = idRep.split("_");
		String trId = f[0] + "_" + f[1];
		String chr = f[2];
		int start = Gpr.parseIntSafe(f[3]);
		int end = Gpr.parseIntSafe(f[4]);
		int idx = Gpr.parseIntSafe(f[5]);

		//---
		// Calculate position witihn CDS
		//---

		// Find transcript and exon
		Transcript tr = trancriptById.get(trId);
		if (tr == null) return null;
		Exon ex = tr.findExon(start);
		if (ex == null) return null;

		// Calculate start position
		int startPos;
		int fr = 0;
		if (ex.getFrame() != 0) {
			if (ex.isStrandPlus()) {
				idx--;
				if (ex.getFrame() == 2) idx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			} else {
				idx--;
				if (ex.getFrame() == 2) idx++; // I don't know why UCSC numbers the AA differentlt when frame is 2
				fr = 3 - ex.getFrame(); // Offset based on frame
			}
		}

		// Find AA start position
		if (ex.isStrandPlus()) {
			int exStart = Math.max(start, tr.getCdsStart());
			startPos = exStart + (idx * 3 + fr);
		} else {
			int exEnd = Math.min(end, tr.getCdsStart());
			startPos = exEnd - (idx * 3 + fr);
		}

		// Get position within CDS
		int cdsBase = tr.baseNumberCds(startPos, false);
		int cds2pos[] = tr.baseNumberCds2Pos();
		if ((ex.isStrandPlus() && (startPos < ex.getStart())) //
				|| (ex.isStrandMinus() && (startPos > ex.getEnd()))) {
			// If the position is outside the exon, then we must jump to previous exon
			startPos = cds2pos[cdsBase - ex.getFrame()];
			cdsBase = tr.baseNumberCds(startPos, true);
		}

		//---
		// Sanity check: Make sure that AA matches between transcript model and MSA data from 'genes likelihood' file
		//---

		// Extract codon
		String cdsSeq = tr.cds();
		String codonStr = cdsSeq.substring(cdsBase, cdsBase + 3);
		String aa = genome.codonTable().aa(codonStr);

		if (aa.equals("" + aaExpected)) {
			countOk++;
			if (debug) Gpr.debug("OK: " + id + " : " + aa);
			else System.out.print('.'); // Show OK
		} else {
			countErr++;
			if (debug) Gpr.debug("Entry ID     : " + id //
					+ "\ntr ID        : " + trId + ", chr: " + chr + ", start: " + start + ", end: " + end + ", idx: " + idx + ", fr: " + fr//
					+ "\nTranscript : " + tr //
					+ "\nExon       : " + ex //
					+ "\nStart pos: " + startPos //
					+ "\nCodon    : " + codonStr + ", aa (real): " + aa + ", aa (exp): " + aaExpected //
			);
			else System.out.print('*'); // Show error
		}

		// Add a newline every now and then
		if (!debug && (countOk + countErr) % 100 == 0) System.out.println("");

		//---
		// Create marker
		// Important: The marker has ALL bases in the codon. 
		//            For instance, is the codon is split between two exons, the 
		//            marker will contain the intron
		//---
		int markerStart, markerEnd;
		if (tr.isStrandPlus()) {
			markerStart = cds2pos[cdsBase];
			markerEnd = cds2pos[cdsBase + 2];
		} else {
			markerStart = cds2pos[cdsBase + 2];
			markerEnd = cds2pos[cdsBase];
		}

		Chromosome chromo = genome.getChromosome(chr);
		Marker marker = new Marker(chromo, markerStart, markerEnd, ex.isStrandMinus(), id);
		return marker;
	}

	/**
	 * Perform GWAS analysis using epistatic data
	 */
	public void readGenesLogLikelihood() {

		//---
		// Read "genes likelihood" file
		//---
		Timer.showStdErr("Reading genes likelihood file '" + genesLikeFile + "'.");
		int count = 0, countKept = 0;
		LineFileIterator lfi = new LineFileIterator(genesLikeFile);
		for (String line : lfi) {
			if (line.isEmpty()) continue;

			// Parse line
			String f[] = line.split("\t");
			String msa1 = f[0];
			String msa2 = f[1];
			double logLikRatio = Gpr.parseDoubleSafe(f[2]);
			String seq1 = f[5];
			String seq2 = f[6];

			// Filter by log likelihood
			count++;
			if (logLikRatio < LOG_LIKELIHOOD_RATIO_THRESHOLD) continue;

			// Map to genomic coordinates
			// if (debug) Gpr.debug(lfi.getLineNum() + ": " + line);
			countKept++;

			parseMsaId(msa1, seq1.charAt(0));
			parseMsaId(msa2, seq2.charAt(0));

		}

		int tot = countErr + countOk;
		Timer.showStdErr("Genes likelihood file '" + genesLikeFile + "'." //
				+ "\n\tEntries kept: " + countKept + " / " + count + " [ " + (countKept * 100.0 / count) + "% ]" //
				+ "\n\tmapping. Err / OK : " + countErr + " / " + tot + " [ " + (countErr * 100.0 / tot) + "% ]" //
		);
	}

}
