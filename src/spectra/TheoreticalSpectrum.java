package spectra;

import java.util.Map;
import java.util.TreeMap;

import utility.AminoAcids;

/**
 * Son class for theoretical spectra that implement specific things
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class TheoreticalSpectrum extends SGSpectrum {
	// Attributes
	/**
	 * The amino acid sequence of the peptide
	 */
	private String _peptideSequence;

	// Constructor

	/**
	 * Empty constructor for Theoretical Spectrum object
	 */
	public TheoreticalSpectrum() {

	}

	/**
	 * The constructor for Theoretical spectrum that contain mass list generated
	 * with the peptide sequence
	 * 
	 * @param peptide                  : The peptide sequence corresponding to PSM
	 * @param experimentalSpectrumData : Data from the experimental Spectrum
	 */
	public TheoreticalSpectrum(String peptide) {

		super(generatePeakList(peptide));
		setPeptideSequence(peptide);
		setMainMass(pepMassCalculator(peptide));

	}

	// Operators
	/**
	 * This is for calculated all theoretical perfect peak that this molecule
	 * produce in mass spectrometry
	 * 
	 * @param sequence : The peptide sequence
	 * @return : peak list in Map<m/z, intensity> format
	 */
	public static Map<Double, Double> generatePeakList(String sequence) {
		Map<Double, Double> tempMap = new TreeMap<>();

		for (int i = 0; i <= sequence.length(); i++) {
			// b peaks
			tempMap.put((AminoAcids.getSequenceMass(sequence.substring(0, i)) + AminoAcids.getUnitMass("H+")), 10000.0);
			// y peaks
//			tempMap.put((AminoAcids.getPeptideMass(sequence.substring(i)) + AminoAcids.getAAMass("H+")
//					+ AminoAcids.getAAMass("NT") + AminoAcids.getAAMass("CT")), 10000.0);
		}
		return tempMap;
	}

	public Double pepMassCalculator(String sequence) {
		// add values of mass for CT, NT and the 2 H+
		Double pepMass = AminoAcids.getUnitMass("NT") + AminoAcids.getUnitMass("CT");
		for (char aminoacid : sequence.toCharArray()) {
			pepMass += AminoAcids.getUnitMass(String.valueOf(aminoacid));
		}

		return pepMass;
	}

	@Override
	public String toString() {

		StringBuilder querry = new StringBuilder(getPeptideSequence() + "\n");

		for (Double mass : getMassList()) {
			querry.append(mass + "\t" + getPeakList().get(mass) + "\n");
		}
		return querry.toString();
	}

	// Getters and Setters
	public String getPeptideSequence() {
		return _peptideSequence;
	}

	public void setPeptideSequence(String peptideSequence) {
		_peptideSequence = peptideSequence;
	}

}
