package comparison;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.CountDownLatch;

import spectra.ExperimentalSpectrum;
import spectra.TheoreticalSpectrum;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzdata_wrapper.MzMlWrapper;
import utility.AminoAcids;
import utility.InputCSVLoader;
import utility.SGXProperties;

/**
 * Class that load files (mgf, csv and parameters) and where all task are
 * launched
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SpecGlobX {

	// Global variable

	// Amino Acids Modifications:
	/**
	 * This is the HashMap of all modification which are present and set from config
	 * file
	 */
	public static HashMap<String, Double> AA_MODIFS = new HashMap<>();
	static {
		AA_MODIFS.put("G", 0.0);
		AA_MODIFS.put("A", 0.0);
		AA_MODIFS.put("S", 0.0);
		AA_MODIFS.put("P", 0.0);
		AA_MODIFS.put("V", 0.0);
		AA_MODIFS.put("T", 0.0);
		AA_MODIFS.put("C", 0.0);
		AA_MODIFS.put("I", 0.0);
		AA_MODIFS.put("L", 0.0);
		AA_MODIFS.put("N", 0.0);
		AA_MODIFS.put("D", 0.0);
		AA_MODIFS.put("Q", 0.0);
		AA_MODIFS.put("K", 0.0);
		AA_MODIFS.put("E", 0.0);
		AA_MODIFS.put("M", 0.0);
		AA_MODIFS.put("H", 0.0);
		AA_MODIFS.put("F", 0.0);
		AA_MODIFS.put("R", 0.0);
		AA_MODIFS.put("Y", 0.0);
		AA_MODIFS.put("W", 0.0);
		AA_MODIFS.put("U", 0.0);
		AA_MODIFS.put("O", 0.0);
	}

	// attributes

	/**
	 * The file Object where are data of scans . It can be .mgf or .mzml, and will
	 * be differentiate after.
	 */
	private File _scanFile;

	/**
	 * Indication of the spectra file format
	 */
	private String _dataFormat;

	/**
	 * The file Object that correspond to where will be stored information of
	 * alignment
	 */
	private String _outputFilePath="output.csv";

	/**
	 * The object that contain informations about the CSV file, path and column id
	 * that contain needed informations
	 */
	private InputCSVLoader _infoFileCSV;

	/**
	 * Object that contain all data and information of all experimental spectra in
	 * the mgf or mzml file.
	 */
	private JMzReader _experimentalSpectraData;
	private MgfFile _expeSpectraMgf;

	/**
	 * The maximum number of peak in all loaded spectra to define the size of matrix
	 */
	private int _maxLengthSpectrum = 2000;

	/**
	 * A map that associate the spectrum Title and the given ID to call it
	 */
	private HashMap<String, Integer> _idScansMap;

	/**
	 * program launched in command mode or in GUI mode
	 */
	private static boolean _commandMode = false;

	/**
	 * The command mode ms input file path (MGF of MZML)
	 */
	private static String _msFilePath = "";

	/**
	 * The command mode csv info file path
	 */
	private static String _csvFilePath = "";

	/**
	 * The command set columns with the properties file
	 */
	private int _scanIDColumn;
	private int _peptideSeqColumn;

	// constructor

	/**
	 * The constructor of the SpecGlob Object that load the spectra scan data and
	 * launch alignment. It is used with the command mode
	 * 
	 * @param args : arguments of the command line
	 * @throws JMzReaderException
	 */
	public SpecGlobX(String[] args) throws JMzReaderException {

		SGXProperties.setConfigFromFile();
		if (SGXProperties.FILTER_TYPE == 1) {
			setMaxLengthSpectrum(SGXProperties.N_MOST_INTENSE * 2);
		}

		String titleIDcol = "a";
		String pepIDcol = "b";

		for (int i = 0; i < args.length; i++) {
			switch (args[i]) {
			case "-msfile":
				i++;
				setMsFilePath(args[i]);
				setScanFile(new File(args[i]));
				break;
			case "-csvfile":
				i++;
				setCsvFilePath(args[i]);
				break;
			case "-outfile":
				i++;
				setOutputFile(args[i]);
				break;
			case "-pepcol":
				i++;
				pepIDcol = args[i];
				break;
			case "-titlecol":
				i++;
				titleIDcol = args[i];
				break;
			}
		}

		setExperimentalSpectraData(commandLoadScanData());
		setInfoFileCSV(new InputCSVLoader(getCsvFilePath(), titleIDcol, pepIDcol));

		associateIDTitleScan();

	}

	/**
	 * SpecGlob Constructor that be in association with the application interface.
	 * It will take arguments to create the instance of alignment.
	 * 
	 * @param dataFile     : the File object of Spectra data (mgf or mzml)
	 * @param csvFile      : the File object of the list of pair Titlescan / peptide
	 *                     PSM
	 * @param outputFile   : The String of the path of the output File
	 * @param dataFileType : the type of spectra file (mgf or mzml)
	 * @param parallel     : boolean that if multiprocess is true or false
	 * @param nbThread     : nb of thread used for multiprocess
	 * @param filterType   : The int value of the filter 0 (intensity rate) or 1 (n
	 *                     most intense peak)
	 * @param valueFilter  : the value associated to the filter % intensity or n
	 *                     quantity of peaks
	 * @throws JMzReaderException
	 */
	public SpecGlobX(File dataFile, String csvFilePath, String titleScanCol, String peptideCol, String outputFile,
			String dataFileType, boolean parallel, byte nbThread, byte filterType, int valueFilter, double precision)
			throws JMzReaderException {

		SGXProperties.setConfigFromFile();

		// setup the Scan file data
		setScanFile(dataFile);
		setDataFormat(dataFileType);
		switch (dataFileType) {
		case "MGF":
			setExperimentalSpectraData(new MgfFile(dataFile));
			setExpeSpectraMgf(new MgfFile(dataFile));
			if (SpecGlobXGUI.commandMode)
				System.out.println("MGF File loaded");
			else
			     SpecGlobXGUI.LOG.append("MGF File Loaded\n");
			break;
		case "MZML":
			setExperimentalSpectraData(new MzMlWrapper(dataFile));
			if (SpecGlobXGUI.commandMode)
				System.out.println("mzML file loaded");
				else
				SpecGlobXGUI.LOG.append("mzML File Loaded\n");
			break;
		}

		setInfoFileCSV(new InputCSVLoader(csvFilePath, titleScanCol, peptideCol));

		setOutputFile(outputFile);

		SGXProperties.FILTER_TYPE = filterType;
		switch (filterType) {
		case 0:
			SGXProperties.INTENSITY_RATE = (byte) valueFilter;
			break;
		case 1:
			SGXProperties.N_MOST_INTENSE = valueFilter;
			setMaxLengthSpectrum(valueFilter * 2);
			break;
		}

		SGXProperties.IS_PARALLELIZED = parallel;
		SGXProperties.NB_THREADS = nbThread;

		SGXProperties.PRECISION = precision;

		associateIDTitleScan();
	}

	/**
	 * This is the load file function to use when it is launch in command mode
	 * 
	 * @return The JMzReader Object
	 * @throws JMzReaderException
	 */
	public JMzReader commandLoadScanData() throws JMzReaderException {
		File selectedFile = new File(getMsFilePath());

		// check if it is MGF or MZML
		if (getMsFilePath().endsWith(".mgf")) {
			System.out.println("I'm a MGF file !!!");
			setDataFormat("MGF");
			MgfFile mgfData = new MgfFile(selectedFile);
			setExpeSpectraMgf(mgfData);
			return mgfData;

		} else if (getMsFilePath().endsWith(".mzML")) {
			setDataFormat("MZML");
			System.out.println("I'm a mzML file !!!");
			return new MzMlWrapper(selectedFile);

		} else
			System.out.println("NOT VALID DATA FORMAT");
		System.exit(0);
		return (null);
	}

	/**
	 * Create an instance of SpectralAlignment and pass all spectra couple inside
	 * 
	 * @throws JMzReaderException
	 */
	public void launchAlignments() throws JMzReaderException, FileNotFoundException {

		SpectralAlignment specAlign = new SpectralAlignment(null, null, getMaxLengthSpectrum() * 2);

		Path pathToFile = Paths.get(getInfoFileCSV().getFilePath().getAbsolutePath());

		long nbLine = countLineNumberReader() - 1;
		if (SpecGlobXGUI.commandMode)
			System.out.println("There are " + nbLine + " alignments to do !\n");
		else
			SpecGlobXGUI.LOG.append("There are " + nbLine + " alignments to do !\n");

		float progressStep = (float) (500.0 / nbLine);
		System.out.println("progressInLaunchAlignments "+progressStep);
		int actualAlign = 0;

		// Here is the format of reading
		try (BufferedReader br = Files.newBufferedReader(pathToFile, StandardCharsets.UTF_8)) {
			// this is for write in the CSV as things progress
			try (PrintWriter writerCSV = new PrintWriter(getOutputFile())) {

				writerCSV.write(
						"Title;Peptide;MassDelta;SharedPeaksBeforeAlign;SharedPeaksAfterAlign;PreAlignedPeptide;AlignedPeptide;NbShift;NotAlignedMass;ScoreAlign;IntensityExplained\n");

				String prevTitleScan = "";

				// read the first line from the text file
				String line = br.readLine();
				// we pass the first line to avoid interpret first line with titles as
				// information for scan selection
				line = br.readLine();

				// loop until all lines are read
				while (line != null) {

					actualAlign += 1;
					if (!SpecGlobXGUI.commandMode)
						SpecGlobXGUI.progressBar.setValue(Math.round(400 + actualAlign * progressStep));
					else 
						System.out.println("Progress....." + Math.round(400 + actualAlign * progressStep));

					// use string.split to load a string array with the values from
					// each line of the file, using a ";" as the delimiter
					String[] attributes = line.split(";");
					String titleScan = attributes[getInfoFileCSV().getScanIdColumn()];
					String seqPeptide = attributes[getInfoFileCSV().getProteinSeqColumn()];

					// We check if the title scan is in the list to avoid error in searching key in
					// map
					if (getIDScans().containsKey(titleScan)) {

						if (!titleScan.equals(prevTitleScan)) {
							ExperimentalSpectrum expeSpec = new ExperimentalSpectrum(
									getExperimentalSpectraData().getSpectrumByIndex(getIDScans().get(titleScan)));
							specAlign.setExpeSpec(expeSpec);
						}

						TheoreticalSpectrum theoSpec = new TheoreticalSpectrum(seqPeptide);

						specAlign.setTheoSpec(theoSpec);

						prevTitleScan = titleScan;

						if (SGXProperties.DEBUG_MODE) {
							System.out.println("________________");
							System.out.println("|SPECTRUM TITLE|  :  " + titleScan);
							System.out.println("________________");
							System.out.println(
									"==========START ALIGN OF " + theoSpec.getPeptideSequence() + "==========");
						}

						specAlign.completeAlignment();
					   if (specAlign.getMaxScore() >= SGXProperties.SCORE_MIN_DISPLAY) 
					   {
						// write the result of the alignment in the CSV file only if above minScore
						writerCSV.write(titleScan + ";" + theoSpec.getPeptideSequence() + ";"
								+ specAlign.getFinalResult() + "\n");}
						
						// if the title is not in the map, we indicate it
					} else {
						writerCSV.write(titleScan + ";" + seqPeptide + ";Not Good Title\n");
						if (SpecGlobXGUI.commandMode)
							System.out.println("Title doesn't correspond ... Please check in both files");
						else
							SpecGlobXGUI.LOG.append("Title doesn't correspond ... Please check in both files\n");
					}
					// read next line before looping
					line = br.readLine();
					// if end of file reached, line would be null and stop while
				}
			}
		} catch (

		IOException ioe) {
			if (SpecGlobXGUI.commandMode)
				System.out.println("Can't open result file "+ pathToFile + "\n");
			else
				SpecGlobXGUI.LOG.append("Can't open result file "+ pathToFile + "\n");
			ioe.printStackTrace();
		}

	}

	/**
	 * Create multiple instances of Spectral Alignment in multiple Thread to make
	 * alignments in parallel to decrease execution time
	 * 
	 * @throws InterruptedException That permit to avoid issues for Thread
	 *                              interruption
	 */
	public void parallelAlignmentLaunch() throws InterruptedException {
		// Read input file and cut into different input list for all Threads
		Path pathToFile = Paths.get(getInfoFileCSV().getFilePath().getAbsolutePath());

		byte nbThread = SGXProperties.NB_THREADS;

		long nbLine = countLineNumberReader() - 1;

		// avoid error linked to more line than thread so no input for some threads
		if (nbLine < nbThread)
			nbThread = (byte) nbLine;

		if (SpecGlobXGUI.commandMode)
			System.out.println("There are " + nbLine + " alignments to do !");
		else
			SpecGlobXGUI.LOG.append("There are " + nbLine + " alignments to do !\n");

		if (SpecGlobXGUI.commandMode)
			System.out.println( nbThread + " threads used to execute alignments !");
		else
			SpecGlobXGUI.LOG.append( nbThread + " threads used to execute alignments !\n");

		int step = (int) (nbLine / nbThread);
		long modulo = nbLine % nbThread;
		List<SpectralAlignmentTask> tasks = new ArrayList<>();
		CountDownLatch latch = new CountDownLatch(nbThread);

		// Here is the format of reading
		try (BufferedReader br = Files.newBufferedReader(pathToFile, StandardCharsets.UTF_8)) {

			// read the first line from the text file
			String line = br.readLine();
			// we pass the first line to avoid interpret first line with titles as
			// information for scan selection
			line = br.readLine();

			byte actualThread = 0;

			// loop until all lines are read
			while (line != null) {

				// increment number of threads to get count
				actualThread += 1;

				// to add missing lines if number of line is not a multiple of number of threads
				if (actualThread == nbThread) {
					step += modulo;
				}

				// Initialization of input arrays for title scan and psm
				String[] inputArrayTitle = new String[step];
				String[] inputArraySequence = new String[step];

				// Keep the step number of line in one list to create one thread
				for (int i = 0; i < step; i++) {
					// use string.split to load a string array with the values from
					// each line of the file, using a ";" as the delimiter
					String[] attributes = line.split(";");

					inputArrayTitle[i] = (attributes[getInfoFileCSV().getScanIdColumn()]);
					inputArraySequence[i] = (attributes[getInfoFileCSV().getProteinSeqColumn()]);
					line = br.readLine();

				}

				// initialization and launch of the process in a Thread
				SpectralAlignmentTask task = new SpectralAlignmentTask(inputArrayTitle, getExperimentalSpectraData(),
						getIDScans(), inputArraySequence, new SpectralAlignment(null, null, getMaxLengthSpectrum() * 2),
						latch, step);
				Thread t = new Thread(task);
				t.start();
				tasks.add(task);

			}

			// waiting for all Thread to finish to write results in the CSV file
			latch.await();

			// print the result in a csv file
			try (PrintWriter writerCSV = new PrintWriter(getOutputFile())) {
				writerCSV.write(
						"Title;Peptide;MassDelta;SharedPeaksBeforeAlign;SharedPeaksAfterAlign;PreAlignedPeptide;AlignedPeptide;NbShift;NotAlignedMass;ScoreAlign;IntensityExplained\n");
				for (int i = 0; i < nbThread; i++) {
					writerCSV.write(tasks.get(i).getResult());
				}
			}

		} catch (IOException ioe) {
			if (SpecGlobXGUI.commandMode)
				System.out.println("Can't open result file "+ pathToFile+"\n");
			else
				SpecGlobXGUI.LOG.append("Can't open result file "+ pathToFile + "\n");
			ioe.printStackTrace();
		}
	}

	/**
	 * Function that browses the peakList to found the spectrum with the maximal
	 * length and in parallel associate given ID with TITLE
	 * 
	 * @throws JMzReaderException
	 */
	public void associateIDTitleScan() throws JMzReaderException {

		System.out.println("START ASSOCIATION");
		if (SpecGlobXGUI.commandMode)
			System.out.println("Indexation of spectra");
		else
		     SpecGlobXGUI.LOG.append("Indexation of spectra\n");
		HashMap<String, Integer> tempIdScans = new HashMap<>();

		// we want create a map to associate ID and TITLE for all spectra

		// better way if we are in MGF file, we can get the title better
		if (getDataFormat().equals("MGF")) {
			MgfFile spectraData = getExpeSpectraMgf();

			for (int i = 0; i < spectraData.getMs2QueryCount(); i++) {
				// need i + 1, because indexes are not based same.
				// Spectrum in based to 1 , MS2Query is based to 0.
				tempIdScans.put(spectraData.getMs2Query(i).getTitle(), i + 1);

			}

		} else {
			JMzReader spectraData = getExperimentalSpectraData();

			for (int i = 1; i <= spectraData.getSpectraCount(); i++) {
				tempIdScans.put(String.valueOf(i), i);
			}

		}

		setIDScans(tempIdScans);

		System.out.println("END ASSOCIATION");
		if (SpecGlobXGUI.commandMode)
			System.out.println("Indexation DONE !");
		else
			SpecGlobXGUI.LOG.append("Indexation DONE !\n");
	}

	/**
	 * A function that count the number of line in the input csv file
	 * 
	 * @return the number of read lines in the file
	 */
	public long countLineNumberReader() {

		File file = getInfoFileCSV().getFilePath();
		long lines = 0;

		// LineNumberReader will automatically count number of line in the file
		try (LineNumberReader lnr = new LineNumberReader(new FileReader(file))) {

			while (lnr.readLine() != null)
				;
			lines = lnr.getLineNumber();

		} catch (IOException e) {
			e.printStackTrace();
		}

		return lines;

	}

	/**
	 * Set all modification to apply to amino acids from input in config.properties
	 * file and update new values
	 * 
	 * @param prop
	 */
	public void setModif(Properties prop) {
		AA_MODIFS.put("G", Double.valueOf(prop.getProperty("sg.modif.G")));
		AA_MODIFS.put("A", Double.valueOf(prop.getProperty("sg.modif.A")));
		AA_MODIFS.put("S", Double.valueOf(prop.getProperty("sg.modif.S")));
		AA_MODIFS.put("P", Double.valueOf(prop.getProperty("sg.modif.P")));
		AA_MODIFS.put("V", Double.valueOf(prop.getProperty("sg.modif.V")));
		AA_MODIFS.put("T", Double.valueOf(prop.getProperty("sg.modif.T")));
		AA_MODIFS.put("C", Double.valueOf(prop.getProperty("sg.modif.C")));
		AA_MODIFS.put("I", Double.valueOf(prop.getProperty("sg.modif.I")));
		AA_MODIFS.put("L", Double.valueOf(prop.getProperty("sg.modif.L")));
		AA_MODIFS.put("N", Double.valueOf(prop.getProperty("sg.modif.N")));
		AA_MODIFS.put("D", Double.valueOf(prop.getProperty("sg.modif.D")));
		AA_MODIFS.put("Q", Double.valueOf(prop.getProperty("sg.modif.Q")));
		AA_MODIFS.put("K", Double.valueOf(prop.getProperty("sg.modif.K")));
		AA_MODIFS.put("E", Double.valueOf(prop.getProperty("sg.modif.E")));
		AA_MODIFS.put("M", Double.valueOf(prop.getProperty("sg.modif.M")));
		AA_MODIFS.put("H", Double.valueOf(prop.getProperty("sg.modif.H")));
		AA_MODIFS.put("F", Double.valueOf(prop.getProperty("sg.modif.F")));
		AA_MODIFS.put("R", Double.valueOf(prop.getProperty("sg.modif.R")));
		AA_MODIFS.put("Y", Double.valueOf(prop.getProperty("sg.modif.Y")));
		AA_MODIFS.put("W", Double.valueOf(prop.getProperty("sg.modif.W")));
		AA_MODIFS.put("U", Double.valueOf(prop.getProperty("sg.modif.U")));
		AA_MODIFS.put("O", Double.valueOf(prop.getProperty("sg.modif.O")));
		AminoAcids.updateAAMapMass();
	}

	// Getter and Setters

	/**
	 * @return the file name and path to access
	 */
	public File getScanFile() {
		return _scanFile;
	}

	/**
	 * @param fileName
	 */
	public void setScanFile(File file) {
		_scanFile = file;
	}

	public String getDataFormat() {
		return _dataFormat;
	}

	public void setDataFormat(String dataFormat) {
		_dataFormat = dataFormat;
	}

	public String getOutputFile() {
		return _outputFilePath;
	}

	public void setOutputFile(String outputFile) {
		_outputFilePath = outputFile;
	}

	public InputCSVLoader getInfoFileCSV() {
		return _infoFileCSV;
	}

	public void setInfoFileCSV(InputCSVLoader infoFileCSV) {
		_infoFileCSV = infoFileCSV;
	}

	/**
	 * @return the data passed thought parser (MGF or MZXML) and it's a JMzReader
	 *         object
	 */
	public JMzReader getExperimentalSpectraData() {
		return _experimentalSpectraData;
	}

	/**
	 * @param spectraData
	 */
	public void setExperimentalSpectraData(JMzReader spectraData) {
		_experimentalSpectraData = spectraData;
	}

	public MgfFile getExpeSpectraMgf() {
		return _expeSpectraMgf;
	}

	public void setExpeSpectraMgf(MgfFile expeSpectraMgf) {
		_expeSpectraMgf = expeSpectraMgf;
	}

	public int getMaxLengthSpectrum() {
		return _maxLengthSpectrum;
	}

	public void setMaxLengthSpectrum(int maxlengthSpectrum) {
		_maxLengthSpectrum = maxlengthSpectrum;
	}

	public static boolean isCommandMode() {
		return _commandMode;
	}

	public static void setCommandMode(boolean commandMode) {
		_commandMode = commandMode;
	}

	public static String getMsFilePath() {
		return _msFilePath;
	}

	public static void setMsFilePath(String msFilePath) {
		_msFilePath = msFilePath;
	}

	public static String getCsvFilePath() {
		return _csvFilePath;
	}

	public static void setCsvFilePath(String csvFilePath) {
		_csvFilePath = csvFilePath;
	}

	public HashMap<String, Integer> getIDScans() {
		return _idScansMap;
	}

	public void setIDScans(HashMap<String, Integer> iDScans) {
		_idScansMap = iDScans;
	}

	public int getScanIDColumn() {
		return _scanIDColumn;
	}

	public void setScanIDColumn(int scanIDColumn) {
		_scanIDColumn = scanIDColumn;
	}

	public int getPeptideSeqColumn() {
		return _peptideSeqColumn;
	}

	public void setPeptideSeqColumn(int peptideSeqColumn) {
		_peptideSeqColumn = peptideSeqColumn;
	}

}
