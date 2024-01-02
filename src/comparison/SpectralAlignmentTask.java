package comparison;

import java.util.HashMap;
import java.util.concurrent.CountDownLatch;

import spectra.ExperimentalSpectrum;
import spectra.TheoreticalSpectrum;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import utility.SGXProperties;

/**
 * Creates a thread to manage the alignment of a list of PSMs
 * Uses an instance of SpectraAlignment to align the spectra list provided as a parameter
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SpectralAlignmentTask extends Thread {

	// ATTRIBUTES

	/**
	 * The Object Spectral Alignment that will proceed alignments
	 */
	private SpectralAlignment _specAlign;

	/**
	 * Spectra data
	 */
	private ExperimentalSpectrum[] _experimentalSpectraData;

	/**
	 * A map that associate the spectrum Title and the given ID to call it
	 */
	private HashMap<String, Integer> _idScansMap;

	/**
	 * Ordered List of spectra titles to analyze
	 */
	private String[] _expeSpectraTitles;

	/**
	 * Ordered List of PSM to align with theoretical spectra
	 */
	private String[] _theoSpectraPsm;

	/**
	 * String that will get the final result of all alignments
	 */
	private String _result;

	/**
	 * The countdown for task
	 */
	private CountDownLatch _latch;

	/**
	 * Size of the list
	 */
	private int _size;

	// CONSTRUCTOR

	/**
	 * The constructor of the Object to initialize all needed to do alignments in a
	 * thread
	 * 
	 * @param titleList               : List of spectra titles that must be aligned by this
	 *                                thread
	 * @param experimentalSpectraData : Experimental spectrum list that must be aligned by this
	 *                                thread
	 * @param idScansMap              : the association of Title and Ids for spectra
	 *                                pick in the JMZ Object
	 * @param psmList                 : List of PSM associated to the Titles to
	 *                                align with experimental spectra
	 * @param specAlign               : One instance of the SpectralAlignment class to use
	 *                                for managing alignments in this thread
	 * @param latch                   : Object that permits to indicate when
	 *                                the current Thread has finished
	 * @param size					  : number of spectra that must be aligned by this
	 *                                thread
	 */
	public SpectralAlignmentTask(String[] titleList, ExperimentalSpectrum[] experimentalSpectraData,
			HashMap<String, Integer> idScansMap, String[] seqList, SpectralAlignment specAlign, CountDownLatch latch,
			int size) {
		setExpeSpectraTitles(titleList);
		setExperimentalSpectraData(experimentalSpectraData);
		setIDScans(idScansMap);
		setTheoSpectraPsm(seqList);
		setSpecAlign(specAlign);
		setResult("");
		setLatch(latch);
		setSize(size);

	}

	// OPERATORS

	@Override
	public void run() {
		StringBuilder output = new StringBuilder();
		String prevTitleScan = "";

		float progressStep = (float) (500.0 / getSize());

		for (int i = 0; i < getSize(); i++) {
            if (!SpecGlobXGUI.commandMode)
            	SpecGlobXGUI.progressBar.setValue(Math.round(400 + i * progressStep));
            
			String titleScan = getExpeSpectraTitles(i);
			String psm = getTheoSpectraPsm(i);
			if (getIDScans().containsKey(titleScan)) {

				if (!titleScan.equals(prevTitleScan)) {
					prevTitleScan = titleScan;
				}

				getSpecAlign().setExpeSpec(getExperimentalSpectraData()[i]);
				getSpecAlign().setTheoSpec(new TheoreticalSpectrum(psm));

				getSpecAlign().completeAlignment();
				// if the score is less than filter, result is not write on the output
				if (getSpecAlign().getMaxScore() >= SGXProperties.SCORE_MIN_DISPLAY) {
					output.append(titleScan + ";" + psm + ";" + getSpecAlign().getFinalResult() + "\n");
				}
			} else {
				output.append(titleScan + ";" + psm + ";Not Good Title\n");
				 if (!SpecGlobXGUI.commandMode)
				      SpecGlobXGUI.LOG.append("Title doesn't correspond ... Please check in both files\n");
				 else
					 System.out.println("Title doesn't correspond ... Please check in both files");
			}
			// System.out.println("Finish align number " + i);

		}
		setResult(output.toString());
		getLatch().countDown();
	}

	// GETTERS AND SETTERS
	public SpectralAlignment getSpecAlign() {
		return _specAlign;
	}

	public void setSpecAlign(SpectralAlignment specAlign) {
		_specAlign = specAlign;
	}

	public ExperimentalSpectrum[] getExperimentalSpectraData() {
		return _experimentalSpectraData;
	}

	public void setExperimentalSpectraData(ExperimentalSpectrum[] experimentalSpectraData) {
		_experimentalSpectraData = experimentalSpectraData;
	}

	public HashMap<String, Integer> getIDScans() {
		return _idScansMap;
	}

	public void setIDScans(HashMap<String, Integer> iDScans) {
		_idScansMap = iDScans;
	}

	public String[] getExpeSpectraTitles() {
		return _expeSpectraTitles;
	}

	public String getExpeSpectraTitles(int i) {
		return _expeSpectraTitles[i];
	}

	public void setExpeSpectraTitles(String[] expeSpectraTitles) {
		_expeSpectraTitles = expeSpectraTitles;
	}

	public String[] getTheoSpectraPsm() {
		return _theoSpectraPsm;
	}

	public String getTheoSpectraPsm(int i) {
		return _theoSpectraPsm[i];
	}

	public void setTheoSpectraPsm(String[] theoSpectraPsm) {
		_theoSpectraPsm = theoSpectraPsm;
	}

	public String getResult() {
		return _result;
	}

	public void setResult(String result) {
		_result = result;
	}

	public CountDownLatch getLatch() {
		return _latch;
	}

	public void setLatch(CountDownLatch latch) {
		_latch = latch;
	}

	public int getSize() {
		return _size;
	}

	public void setSize(int size) {
		_size = size;
	}
}
