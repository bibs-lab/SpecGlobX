package matrix;

/**
 * This is the class for the origin matrix. It will be use to stock informations
 * about where come from the score in the case and how it will be obtain
 * (Alignment (2 - AL) , Re Alignment (1 - RE) , NON Alignment (0 - NA)
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class MatrixOrigin {

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
	 * a 2D array of int[] that contain in int[0] the score origin column indices
	 * and in int[1] the type of alignment coded in 0 (Non align), 1(Re-align) and
	 * 2(Align)
	 */
	private int[][][] _dataFrame;

	// Constructor
	/**
	 * A constructor that is used to create the origin matrix
	 * 
	 * @param nbRow    number of peak for the bigger theoretical spectrum
	 * @param nbColumn number of peak for the bigger experimental spectrum
	 */
	public MatrixOrigin(int nbRow, int nbCol) {
		setNbColumn(nbCol);
		setNbRow(nbRow);
		_dataFrame = new int[nbRow][nbCol][2];

		for (int j = 1; j < nbCol; j++) {
			setDataAlign(0, j, 1);
			setDataOrigin(0, j, 0);
		}
		for (int i = 1; i < nbRow; i++) {
			setDataAlign(0, i, 0);
			setDataOrigin(0, i, 0);
		}
	}

	// Operators
	/**
	 * Show method to display the whole matrix
	 */
	public void show() {
		for (int i = 0; i < getNbRow(); i++) {
			for (int j = 0; j < getNbColumn(); j++) {
				System.out.printf("[%s,%2d]\t", getDataAlignType(i, j), getDataOrigin(i, j));
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
				System.out.printf("[%s,%d]\t", getDataAlignType(i, j), getDataOrigin(i, j));
			}
			System.out.println();
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
	 * Get an int value 0 , 1 or 2 that correspond to alignment type from origin.
	 * Alignment(2) Re-alignment(1), No-alignment (0)
	 * 
	 * @param row    indices of the row we need to investigate (indices of
	 *               theoretical peak)
	 * @param column indices of the column we need to investigate (indices of
	 *               experimental peak)
	 * @return String value of the type of alignment
	 */
	public int getDataAlign(int row, int column) {
		return _dataFrame[row][column][0];
	}

	/**
	 * Get a String value for the type of alignment for better understanding
	 * 
	 * @param row    indices of the row we need to investigate (indices of
	 *               theoretical peak)
	 * @param column indices of the column we need to investigate (indices of
	 *               experimental peak)
	 * @return String value of the type of alignment
	 */
	public String getDataAlignType(int row, int column) {
		switch (getDataAlign(row, column)) {
		case 0:
			return "NA"; // NON ALIGN
		case 1:
			return "RA"; // RE ALIGN
		case 2:
			return "AL"; // ALIGN
		}
		return "Err";
	}

	/**
	 * Set the Alignment type value in int type in the matrix case [row][column]
	 * Alignment(2) Re-alignment(1), No-alignment (0)
	 * 
	 * @param row    int indices of the row we need to investigate (indices of
	 *               theoretical peak)
	 * @param column int indices of the column we need to investigate (indices of
	 *               experimental peak)
	 * @param value  int the correspondent value for the alignment type :
	 *               Alignment(2) Re-alignment(1), No-alignment (0)
	 */
	public void setDataAlign(int row, int column, int value) {
		_dataFrame[row][column][0] = value;
	}

	/**
	 * Get back the indices of origin column where come from the calculated score
	 * (from row i-1)
	 * 
	 * @param row    int indices of the row we need to investigate (indices of
	 *               theoretical peak)
	 * @param column int indices of the column we need to investigate (indices of
	 *               experimental peak)
	 * @return the origin column of the calculated score at this case
	 */
	public int getDataOrigin(int row, int column) {
		return _dataFrame[row][column][1];
	}

	/**
	 * Set the indices of origin column where come from the calculated score (from
	 * row i-1)
	 * 
	 * @param row    int indices of the row we need to investigate (indices of
	 *               theoretical peak)
	 * @param column int indices of the column we need to investigate (indices of
	 *               experimental peak)
	 * @param value  int indices of the column where come from the calculated score
	 *               in this case
	 */
	public void setDataOrigin(int row, int column, int value) {
		_dataFrame[row][column][1] = value;
	}

}
