package ca.mcgill.pcingola.epistasis;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.io.PDBFileReader;
import org.biojava.bio.structure.io.util.PDBTemporaryStorageUtils.LinkRecord;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * PDB distance analysis
 * 
 * References: http://biojava.org/wiki/BioJava:CookBook:PDB:read
 * 
 * @author pcingola
 */
public class ZzzPdb {

	public static final boolean debug = false;
	public static final int AA_MIN_SEPARATION = 10;

	String pdbFileName;

	public static void main(String[] args) {
		String pdbFileName = Gpr.HOME + "/snpEff/db/pdb/3RGK.pdb";
		ZzzPdb zzz = new ZzzPdb(pdbFileName);
		zzz.run();
	}

	public ZzzPdb(String pdbFileName) {
		this.pdbFileName = pdbFileName;
	}

	/**
	 * Distances within a chain 
	 * @param chain
	 */
	void distance(Chain chain) {
		List<Group> groups = chain.getAtomGroups();

		for (int i = 0; i < groups.size(); i++) {
			double dmin = Double.MAX_VALUE;
			Group gmin = null;
			Group group1 = groups.get(i);

			for (int j = i + AA_MIN_SEPARATION; j < groups.size(); j++) {
				Group group2 = groups.get(j);

				double d = distanceMin(group1, group2);
				if (d < dmin) {
					dmin = d;
					gmin = group2;
				}

			}
			if (gmin != null) System.out.println(dmin + "\t" + group1.getType() + ":" + group1.getChemComp().getId() + "\t" + gmin.getType() + ":" + gmin.getChemComp().getId());
		}
	}

	/**
	 * Distances within all chains in a structure 
	 * @param chain
	 */
	void distance(Structure structure) {
		// Distance
		for (Chain chain : structure.getChains())
			distance(chain);
	}

	/**
	 * Minimum distance between all atoms in two amino acids
	 * @param aa1
	 * @param aa2
	 * @return
	 */
	double distanceMin(Group aa1, Group aa2) {
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
	 * Read link record
	 * @param line
	 * @return
	 */
	protected LinkRecord readLinkLine(String line) {
		String name1 = line.substring(12, 16).trim();
		String altLoc1 = line.substring(16, 17).trim();
		String resName1 = line.substring(17, 20).trim();
		String chainID1 = line.substring(21, 22).trim();
		String resSeq1 = line.substring(22, 26).trim();
		String iCode1 = line.substring(26, 27).trim();

		String name2 = line.substring(42, 46).trim();
		String altLoc2 = line.substring(46, 47).trim();
		String resName2 = line.substring(47, 50).trim();
		String chainID2 = line.substring(51, 52).trim();
		String resSeq2 = line.substring(52, 56).trim();
		String iCode2 = line.substring(56, 57).trim();

		String sym1 = line.substring(59, 65).trim();
		String sym2 = line.substring(66, 72).trim();

		return new LinkRecord(name1, altLoc1, resName1, chainID1, resSeq1, iCode1, name2, altLoc2, resName2, chainID2, resSeq2, iCode2, sym1, sym2);
	}

	/**
	 * Read LINK lines
	 * @return
	 */
	protected List<LinkRecord> readLinkLines(String pdbFileName) {
		try {
			return Files.lines(Paths.get(pdbFileName)) //
					.filter(line -> line.startsWith("LINK ")) //
					.map(line -> readLinkLine(line)) //
					.collect(Collectors.toList()) //
			;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * Run 
	 */
	public void run() {
		PDBFileReader pdbreader = new PDBFileReader();

		try {
			System.out.println("Reading: " + pdbFileName);
			Structure structure = pdbreader.getStructure(pdbFileName);

			List<LinkRecord> links = readLinkLines(pdbFileName);
			links.stream().forEach(System.out::println);

			for (Chain chain : structure.getChains()) {
				System.out.println("Chain: " + chain.getId());

				for (Group group : chain.getAtomGroups()) {
					System.out.println("\tGroup: " + group);

					for (Atom atom : group.getAtoms()) {
						System.out.println("\t\t" + atom);
					}
				}

				distance(chain);
			}

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		System.out.println("Done.");

	}

}
