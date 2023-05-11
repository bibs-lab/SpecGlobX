package utility;

import java.io.File;
import java.util.HashMap;

import javax.swing.Box;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * Class that stock informations about the input CSV file
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class InputCSVLoader {

	// Attributes
	/**
	 * Value that correspond to the Id of the table column that contain information
	 * about scan ID
	 */
	private int _scanIdColumn;

	/**
	 * Value that correspond to the Id of the table column that contain information
	 * about the associated peptide sequence
	 */
	private int _proteinSeqColumn;

	/**
	 * The information about the path of the File to open
	 */
	private File _filePath;

	// Constructor
	public InputCSVLoader() {

		setFilePath(selectCSVFile());

		setColumnsValues();
	}

	public InputCSVLoader(String csvPath, int scanCol, int pepcol) {
		setFilePath(new File(csvPath));
		setScanIdColumn(scanCol);
		setProteinSeqColumn(pepcol);
	}

	/**
	 * constructor if we already get the File Path and the 2 columns ID
	 * 
	 * @param csvPath : The String path of the CSV input file
	 * @param scanCol : The String from JText that contain information about ID of
	 *                the column that contain title Scan ID
	 * @param pepcol  : The String from JText that contain information about ID of
	 *                the column that contain of peptide
	 */
	public InputCSVLoader(String csvPath, String scanCol, String pepcol) {
		setFilePath(new File(csvPath));
		setValuesFromString(scanCol, pepcol);
	}

	/**
	 * Ask to the user to choose a CSV file where are stored information of scans
	 * 
	 * @return The CSV File to use
	 */
	public File selectCSVFile() {

		JFrame frame = new JFrame("Load CSV info file");
		JFileChooser fc = new JFileChooser(new File(System.getProperty("user.home") + "\\Documents"));
		fc.setDialogTitle("Load CSV info file");
		fc.addChoosableFileFilter(new FileNameExtensionFilter("CSV Files", "csv"));
		fc.setAcceptAllFileFilterUsed(false);

		while (true) {
			int result = fc.showOpenDialog(frame);
			if (result == JFileChooser.APPROVE_OPTION) {
				File selectedFile = fc.getSelectedFile();
				System.out.println("Selected file: " + selectedFile.getAbsolutePath());
				if (selectedFile.getAbsolutePath().endsWith(".csv")) {
					System.out.println("It's Okay !!!");
					return selectedFile;
				} else if ("".equals(selectedFile.getAbsolutePath())) {
					System.exit(0);
				} else {
					System.out.println("Not a good file format, try again !!!");
				}
			} else {
				System.exit(0);
			}
		}
	}

	/**
	 * This function will ask to the user the identification of column that contain
	 * needed information informations. Those informations are Spectra IDs
	 * (corresponding to TITLE) and PSM sequences.
	 */
	public void setColumnsValues() {

		JTextField scanIDColumn = new JTextField(5);
		JTextField proteinSeqColumn = new JTextField(5);

		JPanel myPanel = new JPanel();
		myPanel.add(new JLabel("IDspectra:"));
		myPanel.add(scanIDColumn);
		myPanel.add(Box.createHorizontalStrut(15));
		myPanel.add(new JLabel("Peptide:"));
		myPanel.add(proteinSeqColumn);

		int result = JOptionPane.showConfirmDialog(null, myPanel, "Please enter column for ScanID and PSM",
				JOptionPane.OK_CANCEL_OPTION);
		if (result == JOptionPane.OK_OPTION) {
			setValuesFromString(scanIDColumn.getText(), proteinSeqColumn.getText());
		}
		// System.out.println("ScanID = " + scanIDColumn.getText());
		// System.out.println("Protein = " + proteinSeqColumn.getText());
	}

	public void setValuesFromString(String scanIDCol, String peptideCol) {
		HashMap<Character, Integer> alphabet = new HashMap<>();
		int i = 0;
		for (char c = 'A'; c <= 'Z'; ++c) {
			alphabet.put(c, i);
			i++;
		}

		if (isParsable(scanIDCol)) {
			setScanIdColumn(Integer.parseInt(scanIDCol) - 1);
		} else {
			setScanIdColumn(alphabet.get(scanIDCol.toUpperCase().charAt(0)));
		}

		if (isParsable(peptideCol)) {
			setProteinSeqColumn(Integer.parseInt(peptideCol) - 1);
		} else {
			setProteinSeqColumn(alphabet.get(peptideCol.toUpperCase().charAt(0)));
		}

	}

	/**
	 * Check if the value of a String is Parsable into an Integer value
	 * 
	 * @param input a String to check if it can be transform into an integer
	 * @return true or false after check
	 */
	public boolean isParsable(String input) {
		try {
			Integer.parseInt(input);
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
	}

	// Getters and Setters
	public int getScanIdColumn() {
		return _scanIdColumn;
	}

	public void setScanIdColumn(int scanIdColumn) {
		_scanIdColumn = scanIdColumn;
	}

	public int getProteinSeqColumn() {
		return _proteinSeqColumn;
	}

	public void setProteinSeqColumn(int proteinSeqColumn) {
		_proteinSeqColumn = proteinSeqColumn;
	}

	public File getFilePath() {
		return _filePath;
	}

	public void setFilePath(File filePath) {
		_filePath = filePath;
	}
}
