package matrix;

import java.util.List;

/**
 *  Delta mass matrix used to store mass differences between all peaks of the two spectrum
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class MatrixE {

	// Attributes
	/**
	 * Number of columns of the matrix (need the maximal size of experimental
	 * spectra)
	 */
	private int _nbColumn;
	/**
	 * Number of rows of the matrix (need the maximal size of theoretical
	 * spectra)
	 */
	private int _nbRow;
	/**
	 * a 2D array of double that contains mass difference between each peak of
	 * experimental and theoretical spectrum
	 */
	private double[][] _dataFrame;

	// Constructor
	public MatrixE(int nbRow, int nbColumn) {
		setNbColumn(nbColumn);
		setNbRow(nbRow);
		_dataFrame = new double[nbRow][nbColumn];
	}

	// Operators
	/**
	 * This function will be adapted in function of parser used. Function use to
	 * calculate all mass differences and auto fill the matrix with
	 * 
	 * @param theoreticalMassList  list of all masses from theoretical spectrum
	 * @param experimentalMassList list of all masses from experimental spectrum
	 */
	public void fillDMass(List<Double> theoreticalMassList, List<Double> experimentalMassList) {
		for (int i = 0; i < theoreticalMassList.size(); i++) {
			for (int j = 0; j < experimentalMassList.size(); j++) {
				setData(i, j, (experimentalMassList.get(j) - theoreticalMassList.get(i)));
			}
		}
	}

	/**
	 * Show method to display the whole matrix
	 */
	public void show() {
		for (int i = 0; i < getNbRow(); i++) {
			for (int j = 0; j < getNbColumn(); j++) {
				System.out.printf("%4.2f", getData(i, j) + "\t");
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
				System.out.printf("%4.2f", getData(i, j));
				System.out.printf("\t");
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
	 * Get and set functions to put and give back the value inside the case[row,
	 * column] of the data frame.
	 * 
	 * @param row    indices of the row corresponding to a peak rank of the
	 *               Theoretical Spectrum
	 * @param column indices of the column corresponding to a peak rank of the
	 *               Experimental Spectrum
	 * @return the mass offset between peak i and peak j
	 */
	public double getData(int row, int col) {
		return _dataFrame[row][col];
	}

	
	/**
	 * Set the value of mass difference between theoretical peak of indices [row]
	 * dans experimental peak of indices [col] inside the case[row, column] of the
	 * data frame.
	 * 
	 * @param row    indices of the row corresponding to a peak rank of the
	 *               Theoretical Spectrum
	 * @param column indices of the column corresponding to a peak rank of the
	 *               Experimental Spectrum
	 * @param value  mass difference between the two peaks in row and in columns
	 */
	public void setData(int row, int col, double value) {
		_dataFrame[row][col] = value;
	}

}
