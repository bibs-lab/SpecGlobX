package utility;

import java.util.HashMap;

/**
 * This class is used to manipulates masses and get the mass of peptide or amino
 * acids
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 */
public class AminoAcids {

	/**
	 * This HashMap contain the letter corresponding to amino acids in Key with
	 * their associated mass in Value
	 */
	private static HashMap<String, Double> massTable = new HashMap<>();

	static {

		massTable.put("G", 57.021463721083); // Glycine
		massTable.put("A", 71.037113785565); // Alanine
		massTable.put("S", 87.032028405125); // Serine
		massTable.put("P", 97.052763850047); // Proline
		massTable.put("V", 99.068413914529); // Valine
		massTable.put("T", 101.047678469607); // Threonine
		massTable.put("C", 103.009184785565); // Cysteine
		massTable.put("I", 113.084063979011); // Isoleucine
		massTable.put("L", 113.084063979011); // Leucine
		massTable.put("N", 114.042927442166); // Asparagine
		massTable.put("D", 115.026943024685); // Aspartic Acid
		massTable.put("Q", 128.058577506648); // Glutamine
		massTable.put("K", 128.094963016052); // Lysine
		massTable.put("E", 129.042593089167); // Glutamic Acid
		massTable.put("M", 131.040484914529); // Methionine
		massTable.put("H", 137.058911859647); // Histidine
		massTable.put("F", 147.068413914529); // Phenylalanine
		massTable.put("R", 156.101111025652); // Arginine
		massTable.put("Y", 163.063328534089); // Tyrosine
		massTable.put("W", 186.07931295157); // Tryptophan
		massTable.put("U", 168.964198469607); // Selenocysteine
		massTable.put("O", 255.158291550141); // Pyrrolysine
		massTable.put("X", 0.0);
		massTable.put("*", 0.0);

		massTable.put("NT", 1.007825032241); // H
		massTable.put("H+", 1.007276466879); // H+
		massTable.put("CT", 15.99491461956 + 1.007825032241); // OH
		massTable.put("", 0.0);

	}

	private AminoAcids() {

	}

	/**
	 * Function to get the mass of an amino acid for String
	 * 
	 * @param aminoAcidLetter
	 * @return mass of AA
	 */
	public static double getUnitMass(String aminoAcidLetter) {
		if (massTable.containsKey(aminoAcidLetter))
				return massTable.get(aminoAcidLetter);
		else return -1;
	}

	/**
	 * Function to get the mass of an amino acid for char
	 * 
	 * @param aminoAcidLetter
	 * @return mass of AA
	 */
	public static double getUnitMass(char aminoAcidLetter) {
		return massTable.get(String.valueOf(aminoAcidLetter));
	}

	/**
	 * Function to get the mass of a sequence of amino acids
	 * 
	 * @param pepSequence : The peptide sequence in letters
	 * @return mass of the Peptide
	 */
	public static double getSequenceMass(String pepSequence) {
		double massPep = 0;
		int lengthPep = pepSequence.length();
		int i = 0;

		while (i < lengthPep) {
			massPep += getUnitMass(String.valueOf(pepSequence.charAt(i)));
			i += 1;
		}

		return massPep;
	}

	/**
	 * Generalized getMass function to get mass without knowing length of the
	 * sequence
	 * 
	 * @param sequence : an amino acid sequence
	 * @return the mass of the sequence
	 */
	public static double getMass(String sequence) {
		int seqLength = sequence.length();
		if (seqLength == 1) {
			return getUnitMass(sequence);
		} else {
			return getSequenceMass(sequence);
		}

	}

	public static double getYBaseMass() {
		return getUnitMass("CT") + getUnitMass("H+") + getUnitMass("NT");
	}

	public static void updateAAMapMass() {
		massTable.put("G", 57.021463721083 + SGXProperties.AA_MODIFS.get("G")); // Glycine
		massTable.put("A", 71.037113785565 + SGXProperties.AA_MODIFS.get("A")); // Alanine
		massTable.put("S", 87.032028405125 + SGXProperties.AA_MODIFS.get("S")); // Serine
		massTable.put("P", 97.052763850047 + SGXProperties.AA_MODIFS.get("P")); // Proline
		massTable.put("V", 99.068413914529 + SGXProperties.AA_MODIFS.get("V")); // Valine
		massTable.put("T", 101.047678469607 + SGXProperties.AA_MODIFS.get("T")); // Threonine
		massTable.put("C", 103.009184785565 + SGXProperties.AA_MODIFS.get("C")); // Cysteine
		massTable.put("I", 113.084063979011 + SGXProperties.AA_MODIFS.get("I")); // Isoleucine
		massTable.put("L", 113.084063979011 + SGXProperties.AA_MODIFS.get("L")); // Leucine
		massTable.put("N", 114.042927442166 + SGXProperties.AA_MODIFS.get("N")); // Asparagine
		massTable.put("D", 115.026943024685 + SGXProperties.AA_MODIFS.get("D")); // Aspartic Acid
		massTable.put("Q", 128.058577506648 + SGXProperties.AA_MODIFS.get("Q")); // Glutamine
		massTable.put("K", 128.094963016052 + SGXProperties.AA_MODIFS.get("K")); // Lysine
		massTable.put("E", 129.042593089167 + SGXProperties.AA_MODIFS.get("E")); // Glutamic Acid
		massTable.put("M", 131.040484914529 + SGXProperties.AA_MODIFS.get("M")); // Methionine
		massTable.put("H", 137.058911859647 + SGXProperties.AA_MODIFS.get("H")); // Histidine
		massTable.put("F", 147.068413914529 + SGXProperties.AA_MODIFS.get("F")); // Phenylalanine
		massTable.put("R", 156.101111025652 + SGXProperties.AA_MODIFS.get("R")); // Arginine
		massTable.put("Y", 163.063328534089 + SGXProperties.AA_MODIFS.get("Y")); // Tyrosine
		massTable.put("W", 186.07931295157 + SGXProperties.AA_MODIFS.get("W")); // Tryptophan
		massTable.put("U", 168.964198469607 + SGXProperties.AA_MODIFS.get("U")); // Selenocysteine
		massTable.put("O", 255.158291550141 + SGXProperties.AA_MODIFS.get("O")); // Pyrrolysine

		massTable.put("NT", 1.007825032241 + SGXProperties.AA_MODIFS.get("NT")); // NTER modification
		massTable.put("H+", 1.007276466879 ); // H+
		massTable.put("Hy", 1.00782503224);
		massTable.put("CT", 15.99491461956 + 1.007825032241 + SGXProperties.AA_MODIFS.get("CT")); // OH
		massTable.put("", 0.0);

	}

}
