package ca.mcgill.pcingola.epistasis.pdb;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.DBRef;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;

import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.pcingola.epistasis.IdMapper;
import ca.mcgill.pcingola.epistasis.IdMapperEntry;
import ca.mcgill.pcingola.epistasis.coordinates.MsaCoordinates;
import ca.mcgill.pcingola.epistasis.coordinates.PdbCoordinate;
import ca.mcgill.pcingola.epistasis.likelihood.InteractionLikelihood;
import ca.mcgill.pcingola.epistasis.msa.MultipleSequenceAlignmentSet;

public class PdbInteracionAnalysis {

	public static final String UNIPROT_DATABASE = "UNP";

	public static boolean debug = false;
	public static boolean verbose = false;

	int neighbours;
	double distanceThreshold;
	String pdbInteractionsFile;

	IdMapper idMapper;
	InteractionLikelihood interactionLikelihood;
	MultipleSequenceAlignmentSet msas;
	PdbGenomeMsas pdbGenomeMsas;

	Map<String, String> pdbIdChainToMolecule;
	Map<String, List<String>> pdbIdToChains;

	public PdbInteracionAnalysis(InteractionLikelihood interactionLikelihood, double distanceThreshold, int neighbours, String pdbInteractionsFile) {
		this.interactionLikelihood = interactionLikelihood;
		this.distanceThreshold = distanceThreshold;
		this.neighbours = neighbours;
		this.pdbInteractionsFile = pdbInteractionsFile;

		pdbGenomeMsas = interactionLikelihood.getPdbGenomeMsas();
		idMapper = interactionLikelihood.getIdMapper();
		msas = pdbGenomeMsas.getMsas();
	}

	/**
	 * Get AA sequence
	 */
	String aaSequence(Structure pdbStruct, String chain) {
		// AA sequence
		List<AminoAcid> aas = aminoAcids(pdbStruct, chain);
		StringBuilder sb = new StringBuilder();

		for (AminoAcid aa1 : aas)
			sb.append(aa1.getAminoType());

		return sb.toString();
	}

	void addChainToMolecule(String pdbId, String chain, String moleculeName) {
		pdbIdChainToMolecule.put(pdbId + ":" + chain, moleculeName);
	}

	/**
	 * Get list of amino acids in a chain
	 */
	List<AminoAcid> aminoAcids(Structure pdbStruct, String chainName) {
		for (Chain chain : pdbStruct.getChains()) {
			if (chain.getChainID().equals(chainName)) {
				ArrayList<AminoAcid> aas = new ArrayList<AminoAcid>();

				for (Group group : chain.getAtomGroups())
					if (group instanceof AminoAcid) aas.add((AminoAcid) group);

				return aas;
			}
		}

		return null;
	}

	Map<String, String> chainUniprotIds(Structure pdbStruct) {
		// Get db references
		Map<String, String> chain2uniproId = new HashMap<String, String>();
		for (DBRef dbref : pdbStruct.getDBRefs()) {
			if (debug) Gpr.debug("DBREF\tchain:" + dbref.getChainId() + "\tdb: " + dbref.getDatabase() + "\tID: " + dbref.getDbAccession());
			if (dbref.getDatabase().equals(UNIPROT_DATABASE)) chain2uniproId.put(dbref.getChainId(), dbref.getDbAccession());
		}
		return chain2uniproId;
	}

	/**
	 * Minimum distance between 2 amino acids (compare every atom)
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
	 * Analyze interacting sites in a pdb structure
	 */
	void findInteracting(String pdbId, Structure pdbStruct, String chainName1, String chainName2, String trId1, String trid2) {
		List<AminoAcid> aas1 = aminoAcids(pdbStruct, chainName1);
		List<AminoAcid> aas2 = aminoAcids(pdbStruct, chainName2);

		// AA sequence
		if (verbose) {
			Gpr.debug("AA seq 1: " + aaSequence(pdbStruct, chainName1));
			Gpr.debug("AA seq 2: " + aaSequence(pdbStruct, chainName2));
		}

		double sum = 0;
		int count = 0;
		for (AminoAcid aa1 : aas1) {
			MsaCoordinates msa1 = null;
			int aaIdx1 = aa1.getResidueNumber().getSeqNum() - 1;

			for (AminoAcid aa2 : aas2) {
				double dmin = distanceMin(aa1, aa2);
				sum += dmin;
				count++;
				if (dmin <= distanceThreshold) {
					int aaIdx2 = aa2.getResidueNumber().getSeqNum() - 1;

					PdbCoordinate pdb1 = new PdbCoordinate(pdbId, chainName1, aaIdx1, aa1.getAminoType());
					PdbCoordinate pdb2 = new PdbCoordinate(pdbId, chainName2, aaIdx2, aa2.getAminoType());
					if (verbose) Gpr.debug("Interacting\t" + dmin + "\t" + pdb1 + "\t" + pdb2);

					if (msa1 == null) msa1 = pdbGenomeMsas.mapToMsa(pdb1);
					MsaCoordinates msa2 = pdbGenomeMsas.mapToMsa(pdb2);

					if (msa1 == null) {
						Gpr.debug("Could not map PDB coordinates to MSA: " + pdb1);
						continue;
					}
					if (msa2 == null) {
						Gpr.debug("Could not map PDB coordinates to MSA: " + pdb2);
						continue;
					}

					// Calculate LL(MSA)
					String llstr = interactionLikelihood.logLikelihoodRatioStr(msa1.msaId, msa1.msaIdx, msa2.msaId, msa2.msaIdx, false, neighbours);
					if (llstr != null) System.out.println(llstr);
				}
			}
		}

		double avg = count > 0 ? sum / count : Double.POSITIVE_INFINITY;
		System.err.println("Average distance: " + avg + "\t number of pairs: " + count);
	}

	String getMolecule(String pdbId, String chain) {
		return pdbIdChainToMolecule.get(pdbId + ":" + chain);
	}

	/**
	 * Load file show a list of PDB entries having interactions (COMPOUND molecules)
	 * @return A list of pdbIds that have interacting molecules
	 */
	List<String> loadPdbInteracionList() {
		if (verbose) Timer.showStdErr("Loading PDB interactions from '" + pdbInteractionsFile + "'");

		pdbIdChainToMolecule = new HashMap<String, String>();
		pdbIdToChains = new HashMap<String, List<String>>();
		List<String> pdbIdsOk = new ArrayList<String>();

		int countOk = 0;
		String lines[] = Gpr.readFile(pdbInteractionsFile).split("\n");
		for (String line : lines) {
			String fields[] = line.split("\t");

			String pdbId = fields[0].toUpperCase();
			if (verbose) System.out.println(pdbId);

			boolean hasMapping = true;
			Set<String> molecules = new HashSet<>();

			// Make sure all molecules can be mapped to MSAs
			List<String> chains = new LinkedList<String>();
			for (int i = 2; i < fields.length - 1;) {
				String moleculeId = fields[i++];
				String moleculeName = fields[i++];
				i++;
				String chain = fields[i++];
				String organism = fields[i++];

				molecules.add(moleculeId);

				// Split chain IDs
				for (String ch : chain.split(",")) {
					ch = ch.trim(); // Chain IDs
					if (verbose) System.out.println("\t" + moleculeId + ":" + ch + "\t" + moleculeName + "\t" + organism);
					chains.add(ch);

					addChainToMolecule(pdbId, ch, moleculeName);

					// Map the pdbId:chain to an MSA
					List<IdMapperEntry> imes = idMapper.getByPdbId(pdbId, ch);
					if (imes == null || imes.isEmpty()) {
						hasMapping = false;
					} else {
						// Add molecule

						// Show all mappings
						if (verbose) {
							for (IdMapperEntry ime : imes)
								System.out.println("\t\t=>\t" + ime);
						}
					}
				}
			}

			// Is this pdbEntry suitable for interaction analysis?
			boolean ok = (molecules.size() > 1 && hasMapping);
			if (ok) {
				countOk++;
				pdbIdToChains.put(pdbId, chains);
				pdbIdsOk.add(pdbId);
			}

			if (verbose) System.out.println("\t" + (ok ? "OK" : "NO"));
		}

		if (verbose) Timer.showStdErr("Count OK: " + countOk + " / " + lines.length);
		Collections.sort(pdbIdsOk);
		return pdbIdsOk;
	}

	/**
	 * Parse PDB and find sites in contact between different molecules
	 * Note: These are PDB files of COMPOUND molecules
	 */
	void parsePdbFile(String pdbFile, String pdbId) {
		// Make sure Pdb entries map to genome and sequences match
		Set<String> confirmedMappings = new HashSet<String>();
		List<IdMapperEntry> idMaps = pdbGenomeMsas.checkSequencePdbGenome(pdbFile, verbose);

		for (IdMapperEntry ime : idMaps) {
			if (verbose) Gpr.debug("\tMapping:\t" + ime);
			confirmedMappings.add(ime.pdbId.toUpperCase() + ":" + ime.pdbChainId);
		}

		if (confirmedMappings.isEmpty()) {
			if (verbose) Gpr.debug("No mappings found for '" + pdbId + "'");
			return;
		} else if (verbose) Gpr.debug("Confirmed mappings: " + idMaps.size());

		Timer.showStdErr("Parsing PDB file '" + pdbFile + "'");
		Structure pdbStruct = pdbGenomeMsas.readPdbFile(pdbFile);
		if (verbose) Gpr.debug("PDB Structure:\n" + pdbStruct);

		// Get uniprot references
		Map<String, String> chain2uniproId = chainUniprotIds(pdbStruct);

		// Find Pdb chains
		List<String> chains = pdbIdToChains.get(pdbId);
		if (chains == null) throw new RuntimeException("Cannot find chains for pdbId '" + pdbId + "'");

		// Analyze distance between amino acids in different chains
		for (String chain1 : chains) {
			if (!confirmedMappings.contains(pdbId + ":" + chain1)) {
				if (verbose) Gpr.debug("No confirmed mapping available for '" + pdbId + ":" + chain1 + "'");
				continue;
			}

			for (String chain2 : chains) {
				if (chain1.compareTo(chain2) >= 0) continue; // Only calculate once

				// Only different molecules
				String molecule1 = getMolecule(pdbId, chain1);
				String molecule2 = getMolecule(pdbId, chain2);
				if (molecule1.equals(molecule2)) continue;

				// Compare uniprot IDs
				String uniprot1 = chain2uniproId.get(chain1);
				String uniprot2 = chain2uniproId.get(chain2);
				if (uniprot1 != null && uniprot2 != null && uniprot1.equals(uniprot2)) {
					Gpr.debug("Different molecule names but same UNIPROT IDs: '" + uniprot1 + "'\n\t" + molecule1 + "\n\t" + molecule2);
					continue;
				}

				// Get confirmed maps
				if (!confirmedMappings.contains(pdbId + ":" + chain2)) {
					if (verbose) Gpr.debug("No confirmed mapping available for '" + pdbId + ":" + chain2 + "'");
					continue;
				}

				// Get mappings (note that we only get the first available mapping).
				IdMapperEntry ime1 = idMapper.getByPdbId(pdbId, chain1).get(0);
				IdMapperEntry ime2 = idMapper.getByPdbId(pdbId, chain2).get(0);

				// Are the transcripts available in the genome?
				String trId1 = ime1.refSeqId;
				if (pdbGenomeMsas.getTranscript(trId1) == null) {
					if (verbose) Gpr.debug("Cannot find transcript '" + trId1 + "' for " + ime1);
					continue;
				}

				String trId2 = ime2.refSeqId;
				if (pdbGenomeMsas.getTranscript(ime2.refSeqId) == null) {
					if (verbose) Gpr.debug("Cannot find transcript '" + trId2 + "' for " + ime2);
					continue;
				}

				// Do we have MSAs for these transripts?
				if (msas.getMsasByTrId(trId1) == null) {
					if (verbose) Gpr.debug("Cannot find MSAs for transcript '" + trId1 + "'");
					continue;
				}

				if (msas.getMsasByTrId(trId2) == null) {
					if (verbose) Gpr.debug("Cannot find MSAs for transcript '" + trId2 + "'");
					continue;
				}

				// Compute inter-chain distances and likelihoods
				if (verbose) Gpr.debug(pdbId + "\t" + chain1 + ": '" + molecule1 + "' (" + trId1 + ")" + "\t" + chain2 + ": '" + molecule2 + "' (" + trId2 + ")");
				findInteracting(pdbId, pdbStruct, chain1, chain2, trId1, trId2);
			}
		}
	}

	public void run() {
		interactionLikelihood.precalcExps();

		// Find which PDB entries have more than one molecule (and molecules can be mapped to MSA)
		List<String> pdbIds = loadPdbInteracionList();

		// Read PDB entries and find mappings to MSA
		for (String pdbId : pdbIds) {
			// Read pdb file
			String pdbDir = pdbGenomeMsas.getPdbDir();
			String pdbFile = pdbDir + "/" + pdbId.toLowerCase() + ".pdb";
			if (Gpr.exists(pdbFile)) {
				parsePdbFile(pdbFile, pdbId);
			}
		}
	}
}
