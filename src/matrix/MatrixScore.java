package matrix;

/**
 * This is the class for the main Score Matrix It will be use to follow
 * alignment and establish a score to get the best possible alignment way
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class MatrixScore {

	// Attributes
	/**
	 * the number of column of the matrix (need the maximal size of experimental
	 * spectra)
	 */
	private int _nbColumn;
	/**
	 * the number of row of the matrix (need the maximal size of theoretical
	 * spectra)
	 */
	private int _nbRow;
	/**
	 * 2D array of integer that contain score at each step of alignment
	 */
	private int[][] _dataFrame;

	// constructor
	/**
	 * A constructor that is used to construct the score matrix
	 * 
	 * @param nbRow    number of peak for the bigger theoretical spectrum
	 * @param nbColumn number of peak for the bigger experimental spectrum
	 */
	public MatrixScore(int nbRow, int nbColumn) {
		setNbColumn(nbColumn);
		setNbRow(nbRow);
		_dataFrame = new int[nbRow][nbColumn];
	}

	/**
	 * A constructor that is used to construct the score matrix and auto fill the
	 * first column with a chain of NON ALIGN score addition
	 * 
	 * @param nbRow         number of peak for the bigger theoretical spectrum
	 * @param nbColumn      number of peak for the bigger experimental spectrum
	 * @param scoreNonAlign Score that is applied when amino acid is not found in
	 *                      experimental spectrum
	 */
	public MatrixScore(int nbRow, int nbColumn, int scoreNonAlign) {
		setNbColumn(nbColumn);
		setNbRow(nbRow);
		_dataFrame = new int[nbRow][nbColumn];
		// fillFirstCol(scoreNonAlign);
		// fillFirstLine(-1);
	}

	// Operators
	/**
	 * Show method to display the whole matrix
	 */
	public void show() {
		for (int i = 0; i < getNbRow(); i++) {
			for (int j = 0; j < getNbColumn(); j++) {
				System.out.printf("%5d", getData(i, j));
			}
			System.out.println();
		}
	}

	/**
	 * Show method with set of number of row and number of column to show
	 * 
	 * @param nbRow : number of row to show (theoretical size)
	 * @param nbCol : number of column to show (experimental size)
	 */
	public void show(int nbRow, int nbCol) {
		for (int i = 0; i < nbRow; i++) {
			for (int j = 0; j < nbCol; j++) {
				System.out.printf("%5d", getData(i, j));
			}
			System.out.println();
		}
	}

	/**
	 * Fill the first column with non align score to initialize the matrix
	 * 
	 * @param score : The NON ALIGN score
	 */
	public void fillFirstCol(int score) {
		for (int i = 1; i < getNbColumn(); i++) {
			setData(i, 0, score + getData(i - 1, 0));
		}
	}

	/**
	 * Fill the first line with 0 and value to small penalize align from offset
	 * 
	 * @param value : The value to put in fist line
	 */
	public void fillFirstLine(int value) {
		for (int j = 1; j < getNbRow(); j++) {
			setData(0, j, value);
		}
	}

	// Getters and Setters
	/**
	 * Get back the number of column in the matrix
	 * 
	 * @return number of column
	 */
	public int getNbColumn() {
		return _nbColumn;
	}

	/**
	 * Set the number of column in the matrix
	 * 
	 * @param nbCol number of column
	 */
	public void setNbColumn(int nbCol) {
		_nbColumn = nbCol;
	}

	/**
	 * Get back the number of row in the matrix
	 * 
	 * @return number of row
	 */
	public int getNbRow() {
		return _nbRow;
	}

	/**
	 * Set the number of row needed for the matrix
	 * 
	 * @param nbRow number of row
	 */
	public void setNbRow(int nbRow) {
		_nbRow = nbRow;
	}

	/**
	 * Get function to give back the value inside the case[row, column] of the data
	 * frame.
	 * 
	 * @param row    indices of the row corresponding to a peak rank of the
	 *               Theoretical Spectrum
	 * @param column indices of the column corresponding to a peak rank of the
	 *               Experimental Spectrum
	 * @return the score inside the specific case
	 */
	public int getData(int row, int column) {
		return _dataFrame[row][column];
	}

	/**
	 * Set function to put a value inside the case[row, column] of the data frame.
	 * 
	 * @param row    indices of the row corresponding to a peak rank of the
	 *               Theoretical Spectrum
	 * @param column indices of the column corresponding to a peak rank of the
	 *               Experimental Spectrum
	 * @param value  The int score value to put in the case
	 */
	public void setData(int row, int column, int value) {
		_dataFrame[row][column] = value;
	}

}
