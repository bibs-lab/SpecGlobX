package comparison;

import java.util.HashMap;
import java.util.concurrent.CountDownLatch;

import spectra.ExperimentalSpectrum;
import spectra.TheoreticalSpectrum;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import utility.SGXProperties;

/**
 * This class is to create a process that will use an instance of Spectral
 * Alignment to align a list of spectrum in a Thread to allow a multiprocessing
 * by multithreading
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
	private JMzReader _experimentalSpectraData;

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
	 * @param titleList               : List of titles of scans to analyze in this
	 *                                thread
	 * @param experimentalSpectraData : The Spectra datas from the File in JMZreader
	 *                                format
	 * @param idScansMap              : the association of Title and Ids for spectra
	 *                                pick in the JMZ Object
	 * @param psmList                 : List of PSM associated to the Titles to
	 *                                align with experimental spectra
	 * @param specAlign               : The instance of Spectral Alignment to use
	 *                                for doing alignments
	 * @param latch                   : The Object that permit to check status of
	 *                                Thread and to count finished thread to wait
	 *                                all finished before writing results in file
	 */
	public SpectralAlignmentTask(String[] titleList, JMzReader experimentalSpectraData,
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
					try {

						getSpecAlign().setExpeSpec(new ExperimentalSpectrum(
								getExperimentalSpectraData().getSpectrumByIndex(getIDScans().get(titleScan))));
					} catch (JMzReaderException e) {
						e.printStackTrace();
					}
				}

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

	public JMzReader getExperimentalSpectraData() {
		return _experimentalSpectraData;
	}

	public void setExperimentalSpectraData(JMzReader experimentalSpectraData) {
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
