package ca.mcgill.pcingola.epistasis;

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
	public static final boolean debug = false;
	public static final boolean verbose = false;

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
	List<DistanceResult> distance(Chain chain) {
		ArrayList<DistanceResult> results = new ArrayList<>();
		List<AminoAcid> aas = aminoAcids(chain);

		for (int i = 0; i < aas.size(); i++) {
			int minj = i + AA_MIN_SEPARATION;

			for (int j = minj; j < aas.size(); j++) {
				AminoAcid aa1 = aas.get(i);
				AminoAcid aa2 = aas.get(j);

				int aadist = Math.abs(i - j);
				double d = distanceMin(aa1, aa2);
				sumDist[aadist] += d;
				count[aadist]++;

				if (d <= distanceThreshold) {
					countTh[aadist]++;
					DistanceResult dres = new DistanceResult(aa1, aa2, d);
					results.add(dres);
					if (verbose) System.out.println(dres);
				}
			}
		}

		return results;
	}

	/**
	 * Distances within all chains in a structure
	 * @param chain
	 */
	List<DistanceResult> distance(Structure structure) {
		ArrayList<DistanceResult> results = new ArrayList<>();

		// Distance
		for (Chain chain1 : structure.getChains())
			results.addAll(distance(chain1));

		return results;
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
	public List<DistanceResult> run() {
		PDBFileReader pdbreader = new PDBFileReader();
		ArrayList<DistanceResult> results = new ArrayList<>();

		for (IdMapperEntry ime : idMapper.getEntries()) {
			try {
				String pdbFileName = pdbDir + "/" + ime.pdbId.toLowerCase() + ".pdb";

				if (verbose) System.err.println("Distance: " + pdbFileName);
				Structure psbStruct = pdbreader.getStructure(pdbFileName);

				// Within resolution limits? => Process
				double res = psbStruct.getPDBHeader().getResolution();
				if (res > maxResolution) continue;

				// Does it have associated transcripts?
				String pdbId = psbStruct.getPDBCode();
				if (idMapper.getByPdbId(pdbId) == null) continue;

				// Distance
				results.addAll(distance(psbStruct));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		return results;
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
