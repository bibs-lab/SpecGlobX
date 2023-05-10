package spectra;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import utility.AminoAcids;
import utility.SGXProperties;

/**
 * Class that implement specificities to Experimental Spectra
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class ExperimentalSpectrum extends SGSpectrum {

	// Attributes
	/**
	 * Spectrum into output format from JMzReader that contain all needed
	 * information
	 */
	private Spectrum _jmzSpectrumData;

	/**
	 * A peak list that contain information about if it is initial peak,
	 * symmetrized, or both, to adapt the score system. Pair<Intensity, type>
	 */
	private Map<Double, Byte> _peakListSymmetrized;

	// Constructor
	/**
	 * The empty constructor
	 */
	public ExperimentalSpectrum() {

	}

	/**
	 * The constructor that use informations in the JMzReader Spectrum class
	 * 
	 * @param data : data load from the scan file via JMzReader
	 */
	public ExperimentalSpectrum(Spectrum data) {
		super(data.getPeakList());

		setJmzSpectrumData(data);

		// Use the chosen filter on spectrum to select peaks
		if (SGXProperties.FILTER_TYPE == 1) // Select most intense peaks
			filterMostIntense(SGXProperties.N_MOST_INTENSE);
		else
			filterIntensityRate(SGXProperties.INTENSITY_RATE);

		createSymetricPeakList();

		setMainMass(calculateMass());

	}

	// Operator

	/**
	 * Create a SymetricPeakList that contain symmetric peaks and the information if
	 * a peak is the initial peak, the symmetric peak, or both if already exist.
	 */
	public void createSymetricPeakList() {
		Map<Double, Byte> symPeakList = new TreeMap<>();
		Map<Double, Double> tempMap = new TreeMap<>();
		List<Double> tempSymList = new ArrayList<>();
		double tempMass;

		double precursorMass = getJmzSpectrumData().getPrecursorMZ();
		int charge = getJmzSpectrumData().getPrecursorCharge();
		double massProton = AminoAcids.getUnitMass("H+");

		// we do the symmetrization of peaks and put them if there are not already
		// present, and set the type
		for (Map.Entry<Double, Double> entry : intensityReverseClassified(getPeakList()).entrySet()) {
			symPeakList.put(entry.getKey(), (byte) 0);
			tempMap.put(entry.getKey(), entry.getValue());
			tempMass = ((precursorMass * charge) + (2 * massProton) - (charge * massProton)) - entry.getKey();
			tempSymList.add(tempMass);

			if (!SGSpectrum.contains(tempMass, getMassList())) {
				symPeakList.put(tempMass, (byte) 1);
				tempMap.put(tempMass, entry.getValue());
			}
		}

		// we add a peak with NT mass (1.0078) if it is not detected, to give better
		// chance to align the first amino acid in b if it is present
		double ntMass = AminoAcids.getUnitMass("NT");

		if (!SGSpectrum.contains(ntMass, tempSymList)) {
			symPeakList.put(ntMass, (byte) 0);
			tempMap.put(ntMass, 1000.0);
		}

		// We add the B peak corresponding to the precursor (complete peptide)
		double precusorBion = ((precursorMass * charge) + (2 * massProton) - (charge * massProton))
				- AminoAcids.getYBaseMass();

		if (!SGSpectrum.contains(precusorBion, tempSymList)) {
			symPeakList.put(precusorBion, (byte) 0);
			tempMap.put(precusorBion, 1000.0);
		}

		// we check if peak have already his symmetric, that correspond to both b and y.
		for (Double mass : getPeakList().keySet()) {
			if (SGSpectrum.contains(mass, tempSymList)) {
				symPeakList.put(mass, (byte) 2);

			}

		}

		setPeakListSymmetrized(symPeakList);
		setMassList(tempMap);
		// setPeakList(tempMap);

		// System.out.println(" Sym size = " + symPeakList.size() + " - peakList size =
		// " + tempMap.size());

	}

	@Override
	public String toString() {

		StringBuilder querry = new StringBuilder(getJmzSpectrumData().getId() + "\n");

		for (Double mass : getMassList()) {
			querry.append(mass + "\t" + getPeakList().get(mass) + "\n");
		}
		return querry.toString();
	}

	/**
	 * Function for calculate the Experimental mass from precursor and it charge
	 * 
	 * @return
	 */
	public Double calculateMass() {
		int charge = getJmzSpectrumData().getPrecursorCharge();

		return ((getJmzSpectrumData().getPrecursorMZ() * charge) - charge * AminoAcids.getUnitMass("H+"));
	}

	// Getters and Setters
	public Spectrum getJmzSpectrumData() {
		return _jmzSpectrumData;
	}

	public void setJmzSpectrumData(Spectrum jmzSpectrumData) {
		_jmzSpectrumData = jmzSpectrumData;
	}

	public Map<Double, Byte> getPeakListSymmetrized() {
		return _peakListSymmetrized;
	}

	public void setPeakListSymmetrized(Map<Double, Byte> peakListSymmetrized) {
		_peakListSymmetrized = peakListSymmetrized;
	}

}
