package ca.mcgill.pcingola.epistasis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBFileReader;

/**
 * PDB distance analysis
 * 
 * References: http://biojava.org/wiki/BioJava:CookBook:PDB:read
 * 
 * @author pcingola
 */
public class PdbDistanceAnalysis {

	public static final int AA_MIN_SEPARATION = 25;
	public static final int MAX_AA_LEN = 10000;
	public static final boolean debug = true;
	public static final boolean verbose = true;

	String pdbDir;
	double distanceThreshold = 5.0;
	double maxResolution = 3.0;
	double sumDist[];
	int count[];
	int countTh[];
	IdMapper idMapper;

	public PdbDistanceAnalysis(String pdbDir, double distanceThreshold, IdMapper idMapper) {
		this.pdbDir = pdbDir;
		this.distanceThreshold = distanceThreshold;
		this.idMapper = idMapper;
		sumDist = new double[MAX_AA_LEN];
		count = new int[sumDist.length];
		countTh = new int[sumDist.length];
	}

	/**
	 * Get all AAs in a chain
	 * @param chain
	 * @return
	 */
	List<AminoAcid> aminoAcids(Chain chain) {
		ArrayList<AminoAcid> aas = new ArrayList<AminoAcid>();
		for (Group group : chain.getAtomGroups())
			if (group instanceof AminoAcid) aas.add((AminoAcid) group);
		return aas;
	}

	/**
	 * Distances within two chains of the same structure
	 * @param chain1
	 */
	void distance(Chain chain1, Chain chain2) {
		List<AminoAcid> aas1 = aminoAcids(chain1);
		List<AminoAcid> aas2 = aminoAcids(chain2);

		for (int i = 0; i < aas1.size(); i++) {
			int minj = 0;
			if (chain1.getChainID().equals(chain2.getChainID())) minj = i + AA_MIN_SEPARATION;

			for (int j = minj; j < aas2.size(); j++) {
				AminoAcid aa1 = aas1.get(i);
				AminoAcid aa2 = aas1.get(j);

				int aadist = Math.abs(i - j);
				if (aadist == 0) throw new RuntimeException("WTF!?\n\t" + aa1 + "\t" + chain1 + "\n\t" + aa2 + "\t" + chain2);

				double d = distanceMin(aa1, aa2);
				sumDist[aadist] += d;
				count[aadist]++;
				if (d <= distanceThreshold) {
					countTh[aadist]++;
					// Show geneName, transcript_ID, etc.
					if (debug) {
						Structure structure = chain1.getParent();
						String pdbId = structure.getPDBCode();
						String refSeqIds = IdMapper.refSeqIds(idMapper.getByPdbId(pdbId));

						if (verbose) System.out.println("pdbId: " + pdbId //
								+ "\trefSeqIDs: " + refSeqIds //
								+ "\tdistance: " + d //
								+ "\taa1: " + aa1.getChemComp().getId() + ", " + chain1.getChainID() + ", " + aa1.getResidueNumber() + ", " + (i + 1) //
								+ "\taa2: " + aa2.getChemComp().getId() + ", " + chain2.getChainID() + ", " + aa2.getResidueNumber() + ", " + (j + 1) //
						);
					}
				}

			}
		}
	}

	/**
	 * Distances within all chains in a structure 
	 * @param chain
	 */
	void distance(Structure structure) {
		// Distance
		for (Chain chain1 : structure.getChains())
			distance(chain1, chain1);
	}

	/**
	 * Minimum distance between all atoms in two amino acids
	 * @param aa1
	 * @param aa2
	 * @return
	 */
	double distanceMin(AminoAcid aa1, AminoAcid aa2) {
		double distMin = Double.POSITIVE_INFINITY;

		for (Atom atom1 : aa1.getAtoms())
			for (Atom atom2 : aa2.getAtoms()) {
				try {
					double dist = Calc.getDistance(atom1, atom2);
					distMin = Math.min(distMin, dist);
				} catch (StructureException e) {
					throw new RuntimeException(e);
				}
			}

		return distMin;
	}

	/**
	 * Run 
	 */
	public void run() {
		PDBFileReader pdbreader = new PDBFileReader();

		File dir = new File(pdbDir);
		int count = 0, countProcessed = 0;
		for (File pdbFile : dir.listFiles()) {
			try {
				String pdbFileName = pdbFile.getCanonicalPath();
				if (!pdbFileName.endsWith(".pdb")) continue;

				System.err.println("Reading: " + pdbFileName);
				Structure structure = pdbreader.getStructure(pdbFileName);

				// Within resolution limits? => Process
				double res = structure.getPDBHeader().getResolution();

				if (res <= maxResolution) {
					// Distance
					distance(structure);
					countProcessed++;
				} else {
					System.err.println("\tFiltered protein " + structure.getPDBCode() + ". Resolution: " + structure.getPDBHeader().getResolution());
				}

				count++;
			} catch (IOException e) {
				throw new RuntimeException(e);
			}

		}

		System.out.println("Done." //
				+ "\n\tTotal files     : " + count //
				+ "\n\tProcessed files : " + countProcessed //
				+ "\n\tFiltered files  : " + (count - countProcessed) //
		);

	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("pos.dist\tdist\tcount.threshold\tp.threshold\tcount\n");
		for (int i = 0; i < sumDist.length; i++) {
			double dist = 0;
			double pth = 0;
			if (count[i] > 0) {
				dist = sumDist[i] / count[i];
				pth = ((double) countTh[i]) / count[i];
			}

			sb.append(i + "\t" + dist + "\t" + countTh[i] + "\t" + pth + "\t" + count[i] + "\n");
		}

		return sb.toString();
	}
}
