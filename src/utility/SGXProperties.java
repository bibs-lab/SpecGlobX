package utility;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Map.Entry;

import comparison.SpecGlobXGUI;

/**
 * Class that stock properties loaded on the config.properties file
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SGXProperties {

	// Global variables for parameters
	// Variables can be modified in the config.properties file.

	/**
	 * Determines how decimal numbers are written (how many numbers after comma)
	 */
	public static DecimalFormat DECIMALFORM = new DecimalFormat("0.00");
	/**
	 * Fragment Precision with 0.02 by default but can be changed
	 */
	public static double PRECISION = 0.02;

	/**
	 * Filter type chosen : 
	 * 		0 : Peaks under a defined intensity rate are removed
	 *      1 : n most intense peaks are kept (n defined by the user)
	 */
	public static byte FILTER_TYPE = 0;

	/**
	 * Intensity rate if filter type = 0. Default = 2%.
	 */
	public static byte INTENSITY_RATE = 2;

	/**
	 * Number of most intense peaks selected if filter type = 1
	 */
	public static int N_MOST_INTENSE = 60;

	/**
	 * Minimum score above which SpecGlobX returns the alignment
	 */
	public static int SCORE_MIN_DISPLAY = 0;

	/**
	 * Defines if debug information is written (command mode)
	 * we can select more precisely what we want to show.
	 */
	public static boolean DEBUG_MODE = false;
	public static boolean SHOW_MATRICES = true;

	/**
	 * alignment with/without parallelization
	 */
	public static boolean IS_PARALLELIZED = false;

	/**
	 * number of threads useful in parallelization mode
	 */
	public static byte NB_THREADS = 8;

	// Amino Acids Fixe Modifications:
	/**
	 * HashMap of modifications provided in the config.properties
	 * 
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
	
	/**
	 * HashMap with the modification masses
	 */
	
	private static Map<String, String> fixedModifications = new HashMap<String, String>();

	// SCORE PARAMETERS
	/**
	 * Different scores used in the score matrix
	 */
	public static final String _nonAlign = "NonAlign";

	/**
	 * By default, we suggest to differentiate the score depending on the fragmentation pattern 
	 *      . "Native" : the peak was in the original spectrum
	 *      . "Sym" : the peak was added during the "symetrisation" process
	 *      . "Both" the fragmentation pattern contains both original and symetric peaks (b and y_ions in the original spectrum)
	 *      
	 */
	public static final String _reAlignNative = "ReAlignNative";
	public static final String _reAlignSym = "ReAlignSym";
	public static final String _reAlignBoth = "ReAlignBoth";

	public static final String _alignNative = "AlignNative";
	public static final String _alignSym = "AlignSym";
	public static final String _alignBoth = "AlignBoth";
	
	/**
	 * When one peak is absent, the algorithm can miss an alignment and consider it is a shift with a zero mass shift
	 * By default, SpecGlobX considers it as an alignment. But, we can add subtleties and considers that the score might be less that an alignment
	 * We tried this configuration and did not find that the results are really impacted. However, we leave this possibility in the source code
	 */

	private static final String _reAlignNativeNoOffset = "ReAlignNativeNO";
	private static final String _reAlignSymNoOffset = "ReAlignSymNO";
	public static final String _reAlignBothNoOffset = "ReAlignBothNO";

	/**
	 * A map that contain the score to apply according to the type of alignment. The
	 * key are the type of alignment (String) and Values are score. Scores are
	 * redefined by config.properties file
	 */
	public static Map<String, Integer> SCORETOAPPLY;
	static {
		SCORETOAPPLY = new HashMap<>();
		SCORETOAPPLY.put(_nonAlign, -4);
		SCORETOAPPLY.put(_reAlignNative, 2);
		SCORETOAPPLY.put(_reAlignSym, 1);
		SCORETOAPPLY.put(_reAlignBoth, 6);
		SCORETOAPPLY.put(_alignNative, 5);
		SCORETOAPPLY.put(_alignSym, 4);
		SCORETOAPPLY.put(_alignBoth, 10);
	}

	/**
	 * modifies the score on the last amino acid alignment (RA score
	 * rather than NA
	 * Not documented anyMore in the interface
	 */
	
	public static boolean BETTER_END_RA = false;

	// constructor
	private SGXProperties() {

	}

	public static void setConfigFromFile() {
		try (InputStream input = new FileInputStream("config.properties")) {
			Properties prop = new Properties();

			// Load properties
			prop.load(input);

			// get the property value and define them in program
			IS_PARALLELIZED = Boolean.valueOf(prop.getProperty("sg.parallelize"));
			NB_THREADS = Byte.valueOf(prop.getProperty("sg.nbthreads"));
			PRECISION = Double.valueOf(prop.getProperty("sg.precision"));

			INTENSITY_RATE = Byte.valueOf(prop.getProperty("sg.peakIntensityRate"));
			N_MOST_INTENSE = Integer.valueOf(prop.getProperty("sg.peakNumberKeeped"));
			FILTER_TYPE = Byte.valueOf(prop.getProperty("sg.filter"));

			// set Amino Acids mass modifications
			setModif(prop);

			// Set scores
			setScoreToApply(prop);
			BETTER_END_RA = Boolean.valueOf(prop.getProperty("sg.scoreBetterEndRA"));

			SCORE_MIN_DISPLAY = Integer.valueOf(prop.getProperty("sg.scoreMinDisplay"));

			// set the decimal format
			String decForm = "0.0";
			for (int i = 1; i < Integer.valueOf(prop.getProperty("sg.decimalFormat")); i++) {
				decForm += "0";
			}
			DECIMALFORM = new DecimalFormat(decForm);

			DEBUG_MODE = Boolean.valueOf(prop.getProperty("sg.debug"));

		} catch (IOException io) {
			io.printStackTrace();
			SpecGlobXGUI.LOG.append("config.properties file not found in the JAR folder");
		}
	}

	/**
	 * Set all modification to apply to amino acids from input in config.properties
	 * file and update new values
	 * 
	 * @param prop
	 */
	public static void setModif(Properties prop) {
		double massModif = 0.00;
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.G"));
		AA_MODIFS.put("G", massModif);
		if (massModif != 0)
			fixedModifications.put("G", prop.getProperty("sg.modif.G"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.A"));
		AA_MODIFS.put("A", massModif);
		if (massModif != 0)
			fixedModifications.put("A", prop.getProperty("sg.modif.A"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.S"));
		AA_MODIFS.put("S", massModif);
		if (massModif != 0)
			fixedModifications.put("S", prop.getProperty("sg.modif.S"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.P"));
		AA_MODIFS.put("P", massModif);
		if (massModif != 0)
			fixedModifications.put("P", prop.getProperty("sg.modif.P"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.V"));
		AA_MODIFS.put("V", massModif);
		if (massModif != 0)
			fixedModifications.put("V", prop.getProperty("sg.modif.V"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.T"));
		AA_MODIFS.put("T", massModif);
		if (massModif != 0)
			fixedModifications.put("T", prop.getProperty("sg.modif.T"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.C"));
		AA_MODIFS.put("C", massModif);
		if (massModif != 0)
			fixedModifications.put("C", prop.getProperty("sg.modif.C"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.I"));
		AA_MODIFS.put("I", massModif);
		if (massModif != 0)
			fixedModifications.put("I", prop.getProperty("sg.modif.I"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.L"));
		AA_MODIFS.put("L", massModif);
		if (massModif != 0)
			fixedModifications.put("L", prop.getProperty("sg.modif.L"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.N"));
		AA_MODIFS.put("N", massModif);
		if (massModif != 0)
			fixedModifications.put("N", prop.getProperty("sg.modif.N"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.D"));
		AA_MODIFS.put("D", massModif);
		if (massModif != 0)
			fixedModifications.put("D", prop.getProperty("sg.modif.D"));
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.Q"));
		AA_MODIFS.put("Q", massModif);
		if (massModif != 0)
			fixedModifications.put("Q", prop.getProperty("sg.modif.Q"));		
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.K"));
		AA_MODIFS.put("K", massModif);
		if (massModif != 0)
			fixedModifications.put("K", prop.getProperty("sg.modif.K"));		
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.E"));
		AA_MODIFS.put("E", massModif);
		if (massModif != 0)
			fixedModifications.put("E", prop.getProperty("sg.modif.E"));		
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.M"));
		AA_MODIFS.put("M", massModif);
		if (massModif != 0)
			fixedModifications.put("M", prop.getProperty("sg.modif.M"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.H"));
		AA_MODIFS.put("H", massModif);
		if (massModif != 0)
			fixedModifications.put("H", prop.getProperty("sg.modif.H"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.F"));
		AA_MODIFS.put("F", massModif);
		if (massModif != 0)
			fixedModifications.put("F", prop.getProperty("sg.modif.F"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.R"));
		AA_MODIFS.put("R", massModif);
		if (massModif != 0)
			fixedModifications.put("R", prop.getProperty("sg.modif.R"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.Y"));
		AA_MODIFS.put("Y", massModif);
		if (massModif != 0)
			fixedModifications.put("Y", prop.getProperty("sg.modif.Y"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.W"));
		AA_MODIFS.put("W", massModif);
		if (massModif != 0)
			fixedModifications.put("W", prop.getProperty("sg.modif.W"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.U"));
		AA_MODIFS.put("U", massModif);
		if (massModif != 0)
			fixedModifications.put("U", prop.getProperty("sg.modif.U"));	
		
		massModif = Double.valueOf(prop.getProperty("sg.modif.O"));
		AA_MODIFS.put("O", massModif);
		if (massModif != 0)
			fixedModifications.put("O", prop.getProperty("sg.modif.O"));	
		
		AminoAcids.updateAAMapMass();
	}
	
	/**
	 * 
	 * @return ArrayList<String> : list of amino acid defined with a fixed modification
	 */
	
	public static ArrayList<String> getAAFixedModifications() {
		ArrayList<String> aaList = new ArrayList<String>();
		Iterator<Entry<String, String>> iterator = fixedModifications.entrySet().iterator();
		while (iterator.hasNext()) {
			Map.Entry<String, String> entry = iterator.next();
			aaList.add(entry.getKey());
		}
		return aaList;
	}
	
	/*
	 * @param: amino acid
	 * @return: mass of the fixed modification applied on aa or "" if non
	 */
	
	public static String getAAFixedModification(String aa) {
		if (fixedModifications.containsKey(aa))
			return fixedModifications.get(aa);
		else
			return "";
	}

	/**
	 * set values of applied scores in according to the type of alignment
	 * 
	 * @param prop
	 */
	public static void setScoreToApply(Properties prop) {

		Map<String, Integer> tempMap = new HashMap<>();

		tempMap.put(_nonAlign, Integer.valueOf(prop.getProperty("sg.scoreNonAlign")));
		tempMap.put(_reAlignNative, Integer.valueOf(prop.getProperty("sg.scoreReAlignNative")));
		tempMap.put(_reAlignSym, Integer.valueOf(prop.getProperty("sg.scoreReAlignSym")));
		tempMap.put(_reAlignBoth, Integer.valueOf(prop.getProperty("sg.scoreReAlignBoth")));
		tempMap.put(_alignNative, Integer.valueOf(prop.getProperty("sg.scoreAlignNative")));
		tempMap.put(_alignSym, Integer.valueOf(prop.getProperty("sg.scoreAlignSym")));
		tempMap.put(_alignBoth, Integer.valueOf(prop.getProperty("sg.scoreAlignBoth")));
		tempMap.put(_reAlignNativeNoOffset, Integer.valueOf(prop.getProperty("sg.scoreReAlignNativeNO")));
		tempMap.put(_reAlignSymNoOffset, Integer.valueOf(prop.getProperty("sg.scoreReAlignSymNO")));
		tempMap.put(_reAlignBothNoOffset, Integer.valueOf(prop.getProperty("sg.scoreReAlignBothNO")));
		SCORETOAPPLY = tempMap;
	}

}
