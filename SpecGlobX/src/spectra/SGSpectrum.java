package spectra;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import utility.AminoAcids;
import utility.SGXProperties;

/**
 * Mother class of Spectra that implement basic attributes
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SGSpectrum {

	// Attributes

	/**
	 * A Map that contain couple of double corresponding to m/z , Intensity
	 */
	private Map<Double, Double> _peakList;

	/**
	 * A Map that contain couple of double corresponding to m/z , Intensity
	 */
	private Map<Double, Double> _totalPeakList;

	/**
	 * The list of the mass to simplify filling matrix
	 */
	private List<Double> _massList;

	/**
	 * The mass of the peptide. It can be from precursor for experimental or from
	 * the peptide sequence for theoretical
	 */
	private Double _mainMass;

	// Constructor

	/**
	 * An empty constructor for the class GSSpectrum
	 */
	public SGSpectrum() {

	}

	/**
	 * The GSSpectrum constructor which take in arguments a Map of peak
	 * 
	 * @param peakList : Map<Double, Double> that contain all peaks
	 */
	public SGSpectrum(Map<Double, Double> peakList) {
		setPeakList(new TreeMap<>(peakList));
		setMassList(peakList);

	}

	// Operators

	/**
	 * This is used to classify peaks in function of their Intensity value
	 * 
	 * @param peakList List of peak in Map format <m/z , Intensity>
	 * @return Intensity value tried peak list map
	 */
	public Map<Double, Double> intensityReverseClassified(Map<Double, Double> peakList) {

		LinkedHashMap<Double, Double> reversedValueMap = new LinkedHashMap<>();
		peakList.entrySet().stream().sorted(Map.Entry.comparingByValue(Comparator.reverseOrder()))
				.forEachOrdered(x -> reversedValueMap.put(x.getKey(), x.getValue()));

		return reversedValueMap;
	}

	/**
	 * Filter and select peak that have the best value of intensity
	 * 
	 * @param numberOfPeak The number of accepted peak that have the best value of
	 *                     intensity
	 * @return The Map peak list with best peaks
	 */
	public void filterMostIntense(int numberOfPeak) {

		int i = 1;
		TreeMap<Double, Double> tempMap = new TreeMap<>();
		for (Map.Entry<Double, Double> entry : intensityReverseClassified(getPeakList()).entrySet()) {
			if (i > numberOfPeak) {
				setPeakList(tempMap);
				setMassList(tempMap);
				break;
			}
			tempMap.put(entry.getKey(), entry.getValue());
			// System.out.println("Key = " + entry.getKey() + ", Value = " +
			// entry.getValue());
			i++;
		}
	}

	/**
	 * Filter and select peak that have an Intensity superior than the given
	 * Threshold value
	 * 
	 * @param intensityThreshold The Threshold intensity value of accepted peak
	 * @return The Map peak list without peak with an intensity bellow the Threshold
	 */
	public void filterMoreThanIntensity(double intensityThreshold) {
		TreeMap<Double, Double> tempMap = new TreeMap<>();
		for (Map.Entry<Double, Double> entry : intensityReverseClassified(getPeakList()).entrySet()) {
			if (entry.getValue() < intensityThreshold) {
				setPeakList(tempMap);
				setMassList(tempMap);
				break;
			}
			tempMap.put(entry.getKey(), entry.getValue());
			// System.out.println("Key = " + entry.getKey() + ", Value = " +
			// entry.getValue());
		}

	}

	/**
	 * Filter and select peak that have an Intensity superior than the given rate of
	 * the maximal intensity value
	 * 
	 * @param intensityRate The rate of the maximal intensity (in %)
	 * @return The Map peak list without peak with an intensity bellow the value of
	 *         max rate
	 */
	public void filterIntensityRate(int intensityRate) {
		int i = 0;
		Double bestIntensity = 1000.0;
		TreeMap<Double, Double> tempMap = new TreeMap<>();
		for (Map.Entry<Double, Double> entry : intensityReverseClassified(getPeakList()).entrySet()) {
			if (i == 0) {
				bestIntensity = entry.getValue();
				// System.out.println(bestIntensity * intensityRate / 100d);
				i++;
			}
			if (entry.getValue() < bestIntensity * intensityRate / 100d) {
				setPeakList(tempMap);
				setMassList(tempMap);
				break;
			}
			tempMap.put(entry.getKey(), entry.getValue());

			// System.out.println("Key = " + entry.getKey() + ", Value = " +
			// entry.getValue());
		}

	}

	// Utility method

	/**
	 * A static method useful to check if the listMass contain the proposed mass
	 * value and to remove some peaks
	 * 
	 * @param value     : The mass value to check if it is already inside the
	 *                  massList
	 * @param list      : The massList to check in
	 * @param precision : The precision of values in the scan
	 * @return
	 */
	public static boolean contains(Double value, List<Double> list) {
		list.add(AminoAcids.getYBaseMass());
		for (Double e : list) {
			if (value.equals(e) || Math.abs(value - e) < SGXProperties.PRECISION) {
				// for debug:
				// System.out.println("We are same at a precision of " + getPrecision() + " : "
				// + value + " - " + e);
				return true;
			}
		}
		return false;
	}

	// Getters and Setters
	public Map<Double, Double> getPeakList() {
		return _peakList;
	}

	public void setPeakList(Map<Double, Double> map) {
		_peakList = map;
	}

	public Map<Double, Double> getTotalPeakList() {
		return _totalPeakList;
	}

	public void setTotalPeakList(Map<Double, Double> totalPeakList) {
		_totalPeakList = totalPeakList;
	}

	public List<Double> getMassList() {
		return _massList;
	}

	public void setMassList(List<Double> massList) {
		_massList = massList;
	}

	/**
	 * Create the Mass list from the peakList
	 * 
	 * @param peakList
	 */
	public void setMassList(Map<Double, Double> peakList) {
		List<Double> massList = new ArrayList<>(peakList.keySet());
		_massList = massList;
	}

	public Double getMainMass() {
		return _mainMass;
	}

	public void setMainMass(Double mainMass) {
		this._mainMass = mainMass;
	}

}
