package comparison;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import matrix.MatrixE;
import matrix.MatrixOrigin;
import matrix.MatrixScore;
import spectra.ExperimentalSpectrum;
import spectra.TheoreticalSpectrum;
import utility.AminoAcids;
import utility.SGXProperties;

/**
 * The Spectral Alignment class where are calculated scores and where the
 * Modified hit is estimated and generated
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SpectralAlignment {

	// Attributes
	/**
	 * Non Align String to access to the good score to apply
	 */
	private static final String _nonAlign = "NonAlign";

	/**
	 * Re Align String to access to the good score to apply
	 */
	private static final String _reAlignNative = "ReAlignNative";
	private static final String _reAlignSym = "ReAlignSym";
	private static final String _reAlignBoth = "ReAlignBoth";

	/**
	 * Align String to access to the good score to apply
	 */
	private static final String _alignNative = "AlignNative";
	private static final String _alignSym = "AlignSym";
	private static final String _alignBoth = "AlignBoth";

	/**
	 * Re Align String to access to the good score to apply
	 */
	private static final String _reAlignNativeNoOffset = "ReAlignNativeNO";
	private static final String _reAlignSymNoOffset = "ReAlignSymNO";
	private static final String _reAlignBothNoOffset = "ReAlignBothNO";

	/**
	 * The score matrix where is stored information of the score during alignment
	 */
	private MatrixScore _matScore;
	/**
	 * The mass offset matrix which contain all mass differences between
	 * experimental and theoretical spectra
	 */
	private MatrixE _matE;
	/**
	 * The matrix where are stored information of where does the alignment come from
	 * and the type of alignment
	 */
	private MatrixOrigin _matOrigin;
	/**
	 * This is the experimental spectrum in which are information of the
	 * spectrometric mass scan and data of all peaks
	 */
	private ExperimentalSpectrum _expeSpec;
	/**
	 * This is the theoretical spectrum created from a peptide sequence (PSM)
	 */
	private TheoreticalSpectrum _theoSpec;

	/**
	 * This is the Hit Modified sequence that contain mass modifications at good
	 * positions
	 */
	private String _hitModifiedSeq;
	/**
	 * This Hit modified is the version that occurs the best number of aligned peaks
	 */
	private String _modifiedAfterBestScore;
	/**
	 * A HitModified with a different sequence structure
	 */
	private String _otherModified;

	/**
	 * number of shared peak between the initial peptide and the experimental
	 * spectrum
	 */
	private int _nbPeakInitial;
	/**
	 * number of shared peak between hit_modified and experimental spectrum
	 */
	private int _nbPeakAfterAlign;
	/**
	 * number of shared peak for the better hit modified after reajust
	 */
	private int _nbReajustedPeak;

	/**
	 * The total mass differences between the hit modified peptide and the
	 * experimental peptide
	 */
	private double _totDeltaMass;

	/**
	 * The theoretical mass difference between the peptide hit sequence and the
	 * given peptide mass by experimental spectrum
	 */
	private double _trueDeltaMass;

	/**
	 * masses that are explained during post treatment and must be removed from the
	 * final residual mass
	 */
	private double _explainedMass;

	/**
	 * Maximal score for the best global alignment
	 */
	private int _maxScore;

	/**
	 * Count of the number of different modifications. Are counted as modification :
	 * Deletion bloc, insertion bloc, substitution
	 */
	private int _modificationNumber;

	private double _confidenceRate;

	/**
	 * The String that contain information of the line to write in the CSV
	 */
	private String _finalResult;

	// Constructor
	/**
	 * The constructor of SpectralAlignment Object with 2 maximal length for
	 * experimental and theoretical spectra
	 * 
	 * @param theoSpec
	 * @param expeSpec
	 * @param maxLengthTheoSpectrum
	 * @param maxLengthExpeSpectrum
	 */
	public SpectralAlignment(TheoreticalSpectrum theoSpec, ExperimentalSpectrum expeSpec, int maxLengthTheoSpectrum,
			int maxLengthExpeSpectrum) {
		setMatScore(new MatrixScore(maxLengthTheoSpectrum, maxLengthExpeSpectrum, getScoreToApply(_nonAlign)));
		setMatOrigin(new MatrixOrigin(maxLengthTheoSpectrum, maxLengthExpeSpectrum));
		setMatE(new MatrixE(maxLengthTheoSpectrum, maxLengthExpeSpectrum));
		setTheoSpec(theoSpec);
		setExpeSpec(expeSpec);

	}

	/**
	 * The constructor of SpectralAlignement Object with only one maximal length for
	 * all spectra, That make a squared matrix
	 * 
	 * @param theoSpec
	 * @param expeSpec
	 * @param maxLengthSpectrum
	 */
	public SpectralAlignment(TheoreticalSpectrum theoSpec, ExperimentalSpectrum expeSpec, int maxLengthSpectrum) {
		this(theoSpec, expeSpec, maxLengthSpectrum, maxLengthSpectrum);

	}

	// Operators
	/**
	 * This method is for doing the complete alignment for all peaks. Fill matrices
	 * at each coordinates and doing the backtrack
	 */
	public void completeAlignment() {
		double precision = SGXProperties.PRECISION;

		setExplainedMass(0);

		int maxScore = -1000;
		int maxPosI = 0;
		int maxPosJ = 0;

		getMatE().fillDMass(getTheoSpec().getMassList(), getExpeSpec().getMassList());

		// Make the alignment in filling the matrices

		for (int i = 1; i < getTheoSpec().getMassList().size(); i++) {
			for (int j = 1; j < getExpeSpec().getMassList().size(); j++) {
				fillMatrices(i, j, precision);
				// keep the Max Score during matrixes filling
				// the condition is to keep the max score in the last line (last theoretical
				// amino acid)
				// This is for try to keep all amino acids of the initial Theoretical sequence
				// (but can modify alignment choice)
				if (i == getTheoSpec().getMassList().size() - 1 && getMatScore().getData(i, j) > maxScore) {
					maxScore = getMatScore().getData(i, j);
					maxPosI = i;
					maxPosJ = j;
				}

			}

		}

		setMaxScore(maxScore);

		if (maxScore >= SGXProperties.SCORE_MIN_DISPLAY) {
			// for debug
			if (SGXProperties.DEBUG_MODE) {
				System.out.println("Best Score = " + getMatScore().getData(maxPosI, maxPosJ) + " - at i = " + maxPosI
						+ " j = " + maxPosJ);
				showMatrices();
			}

			// set the value of difference between Precursor mass and theoretical mass of
			// the proposed peptide
			setTrueDeltaMass(getExpeSpec().getMainMass() - getTheoSpec().getMainMass());

			backTrack(maxPosI, maxPosJ, precision);

			ArrayList<Double> expeMassList = new ArrayList<>(getExpeSpec().getPeakList().keySet());

			double offset = evaluate(getHitModifiedSeq(), expeMassList, precision, getExpeSpec().getMainMass(),
					getTheoSpec().getMainMass());

			setTotDeltaMass(offset);
			
			String sequenceWithFixedModif = getTheoSpec().getPeptideSequence();			
			setNbPeakInitial(numberSharedPeaks(sequenceWithFixedModif, expeMassList, precision));

			// remove brackets for amino acids having their correspondent peak realigned
			setModifiedAfterBestScore(getModifiedAfterBestScore());
			setOtherModified(removeBracketsForAlignedAA(getModifiedAfterBestScore(), expeMassList, precision));
			// set the number of final shared peak after all treatments
			setNbReajustedPeak(numberSharedPeaks(getModifiedAfterBestScore(), expeMassList, precision));
			// calculate a confidence rate that take the number of shared peaks and divide
			// by theoretical numbers of peaks (b and y)
			setConfidenceRate(
					explainedIntensityRate(getModifiedAfterBestScore(), getExpeSpec().getPeakList(), precision));
			// make the output of the actual alignment to store in the process thread of to
			// directly write on the output csv
			makeFinalResult();
		}

	}

	/**
	 * Backtrack is here to get the HitModified Sequence during the backtrack of
	 * score calculation. It start from best score align case in matrix and use
	 * origin matrix to go back to the start of the alignment
	 * 
	 * @param row       : the row indices of the best score
	 * @param column    : the column indices of the best score
	 * @param precision : the precision of measures
	 */
	public void backTrack(int row, int column, double precision) {
		String pepModified = "";
		String theoSequence = getTheoSpec().getPeptideSequence();
		int actualI = row;
		int prevI;
		int actualJ = column;
		int prevJ;
		double actualDelMass;
		double prevDelMass;
		int modifCount = 0;
		double totExplainMass = 0.0;

		// System.out.println(theoSequence);
		// System.out.println("actualI is : " + actualI + " and actualJ is : " +
		// actualJ);

		while (actualI > 0) {
			// System.out.println(actualI);
			// define the actual mass delta to the actual i and j
			String tempPepSeq = "";
			String tempAA = "";
			String aminoAcid = "";
			actualDelMass = getMatE().getData(actualI, actualJ);
			prevI = actualI - 1;
			prevJ = getMatOrigin().getDataOrigin(actualI, actualJ);
			prevDelMass = getMatE().getData(prevI, prevJ);

			// We checking first if last Amino acid are not aligned

			if (getMatOrigin().getDataAlignType(actualI, actualJ).equals("NA")) {
				// System.out.println("I'm NON ALIGN");
				tempPepSeq = "[" + theoSequence.charAt(actualI - 1) + "]";

				while (getMatOrigin().getDataAlignType(prevI, prevJ).equals("NA") && prevI > 0) {
					tempPepSeq = "[" + theoSequence.charAt(prevI - 1) + "]" + tempPepSeq;

					prevI--;
				}
				// modifCount++;
				actualI = prevI;
				actualJ = prevJ;

				pepModified = tempPepSeq + pepModified;

				// if not, we are in the case where there is an alignment or a realignment
				// if there is Align or a re-align, put the letter of Amino Acid and check what
				// we get before to see if we have a mass change
			} else {
				// we put the actual amino acid because he is found, and in function of
				// alignment type, we can have a mass offset
				tempPepSeq = String.valueOf(theoSequence.charAt(actualI - 1));
				aminoAcid = tempPepSeq;

				// we check if there is a deletion before the founded aminoAcid
				if (prevI > 0 && getMatOrigin().getDataAlignType(prevI, prevJ).equals("NA")) {

					// modifCount++;

					// we continue to check if this is a deletion bloc to keep the mass difference
					// due to the deletion of the block
					while (prevI > 0 && getMatOrigin().getDataAlignType(prevI, prevJ).equals("NA")) {

						tempPepSeq = "[" + theoSequence.charAt(prevI - 1) + "]" + tempPepSeq;
						tempAA = theoSequence.charAt(prevI - 1) + tempAA;

						prevI--;
						prevDelMass = getMatE().getData(prevI, prevJ);
						if (prevI == 0) {
							prevDelMass = 0.0;
						}

					}

					// check the mass delta to avoid showing a null mass delta
					if (Math.abs(actualDelMass - prevDelMass) > precision) {

						tempPepSeq = tempPepSeq.substring(0, tempPepSeq.length() - 1);

						tempPepSeq += "["
								+ SGXProperties.DECIMALFORM.format(Double.valueOf(actualDelMass - prevDelMass)) + "]"
								+ aminoAcid;
						totExplainMass += (actualDelMass - prevDelMass);
						modifCount++;

					}

					// if there this is just a re-align, we need to indicate the mass offset
				} else if (getMatOrigin().getDataAlignType(actualI, actualJ).equals("RA")) {

					tempPepSeq = "[" + SGXProperties.DECIMALFORM.format(Double.valueOf(actualDelMass - prevDelMass))
							+ "]" + tempPepSeq;
					modifCount++;
					totExplainMass += (actualDelMass - prevDelMass);

					// the fact when you align the first amino acid, but there is a leak of Amino
					// acid in OMS solution before
				} else if (actualI == 1 && getMatOrigin().getDataAlignType(actualI, actualJ).equals("AL")
						&& (Math.abs(actualDelMass) > precision)) {
					tempPepSeq = "[" + SGXProperties.DECIMALFORM.format(actualDelMass) + "]" + tempPepSeq;
					totExplainMass += actualDelMass;
					modifCount++;
				}

				pepModified = tempPepSeq + pepModified;
				// System.out.println(pepModified);

				actualI = prevI;
				actualJ = prevJ;

			}

		}

		setHitModifiedSeq(
				pepModified + "_[" + SGXProperties.DECIMALFORM.format(getTrueDeltaMass() - totExplainMass) + "]");

		setModificationNumber(modifCount);

	}

	/**
	 * This method do the alignment of the 2 Spectra and fill matrices at actual
	 * coordinates i(row - theoretical) and j(column - experimental)
	 * 
	 * @param theoIndiceI : The actual theoretical peak indices
	 * @param expeIndiceJ : The actual experimental peak indices
	 * @param precision   : The Precision or measures
	 *
	 */
	public void fillMatrices(int theoIndiceI, int expeIndiceJ, double precision) {
		// in first time, we set score corresponding to type of peak (initial/mirror or
		// both)
		// we set initial score to peak mirror or initial
		byte expePeakType = getExpeSpec().getPeakListSymmetrized().get(getExpeSpec().getMassList().get(expeIndiceJ));

		int alignScoreToAdd = getScoreToApply(_alignNative);
		int reAlignScoreToAdd = getScoreToApply(_reAlignNative);
		int reAlignScoreNOToAdd = getScoreToApply(_reAlignNativeNoOffset);

		switch (expePeakType) {

		case 0:
			// this is a native peak
			alignScoreToAdd = getScoreToApply(_alignNative);
			reAlignScoreToAdd = getScoreToApply(_reAlignNative);
			reAlignScoreNOToAdd = getScoreToApply(_reAlignNativeNoOffset);
			if (SGXProperties.BETTER_END_RA && theoIndiceI == getTheoSpec().getPeptideSequence().length()) {
				reAlignScoreToAdd = getScoreToApply(_nonAlign) + getScoreToApply(_alignNative);
			}
			break;
		case 1:
			// this is a symmetric peak
			alignScoreToAdd = getScoreToApply(_alignSym);
			reAlignScoreToAdd = getScoreToApply(_reAlignSym);
			reAlignScoreNOToAdd = getScoreToApply(_reAlignSymNoOffset);
			if (SGXProperties.BETTER_END_RA && theoIndiceI == getTheoSpec().getPeptideSequence().length()) {
				reAlignScoreToAdd = getScoreToApply(_nonAlign) + getScoreToApply(_alignSym);
			}
			break;
		case 2:
			// this is a native peak with a symmetric corresponding
			alignScoreToAdd = getScoreToApply(_alignBoth);
			reAlignScoreToAdd = getScoreToApply(_reAlignBoth);
			reAlignScoreNOToAdd = getScoreToApply(_reAlignBothNoOffset);
			if (SGXProperties.BETTER_END_RA && theoIndiceI == getTheoSpec().getPeptideSequence().length()) {
				reAlignScoreToAdd = getScoreToApply(_nonAlign) + getScoreToApply(_alignBoth);
			}
			break;
		default:
			break;

		}

		int k = getkValue(theoIndiceI, expeIndiceJ, precision);
		if (k == -1) {
			setMatricesData(theoIndiceI, expeIndiceJ,
					getMatScore().getData(theoIndiceI - 1, expeIndiceJ) + (getScoreToApply(_nonAlign)), expeIndiceJ, 0);
		} else {
			int scoreAlignK = getMatScore().getData(theoIndiceI - 1, k) + alignScoreToAdd;
			// if it come from non align, we must verify that there is no offset from
			// previous align
			if (getMatOrigin().getDataAlign(theoIndiceI - 1, k) == 0) {
				int l;
				for (l = theoIndiceI - 1; l > 0; l--) {
					if (getMatOrigin().getDataOrigin(l, k) != k)
						break;

				}
				if (Math.abs(getMatE().getData(l, k) - getMatE().getData(theoIndiceI, expeIndiceJ)) > precision)
					scoreAlignK = getMatScore().getData(theoIndiceI - 1, k) + reAlignScoreToAdd;

			}

			// int[0] = the j value m and int[1] = the score value
			int[] reAlignBestScore = getBestRealignScore(theoIndiceI, k, expeIndiceJ, reAlignScoreToAdd,
					reAlignScoreNOToAdd, precision);

			// For debug to see value for any match
			// System.out.println("score k = " + scoreAlignK + " - score m = " +
			// reAlignBestScore[1] + " - origin m = "
			// + reAlignBestScore[0]);

			if (scoreAlignK >= reAlignBestScore[1]) {
				setMatricesData(theoIndiceI, expeIndiceJ, scoreAlignK, k, 2);

			} else {
				setMatricesData(theoIndiceI, expeIndiceJ, reAlignBestScore[1], reAlignBestScore[0],
						reAlignBestScore[2]);
			}

		}
	}

	/**
	 * This function is use to check if there is a founded AA between experimental
	 * peaks at indices j and indices K Theoretically, the AA is between peaks of
	 * indices i and i-i
	 * 
	 * @param theoIndiceI   the i value (row) of the theoretical peak i
	 * @param expeIndiceJ   the j value (column) of the experimental peak j
	 * @param expePtentialK the potential k value (column) where mass[j] - mass[k] =
	 *                      mass[i] - mass[i-1]
	 * @param precision     : Precisions of measures
	 * 
	 * @return true if mass[j] - mass[k] = mass[i] - mass[i-1] so if AA is found
	 */
	public boolean aaFound(int theoIndicesI, int expeIndicesJ, int expeIndicesK, double precision) {
		double theoMass = (getTheoSpec().getMassList().get(theoIndicesI)
				- getTheoSpec().getMassList().get(theoIndicesI - 1));

		double expeMass = (getExpeSpec().getMassList().get(expeIndicesJ)
				- getExpeSpec().getMassList().get(expeIndicesK));

		return (Math.abs(theoMass - expeMass) < precision);
	}

	/**
	 * Method to found the value of k that give the mass of an amino acid between
	 * experimental masses j and k
	 * 
	 * @param theoIndicesI : the row where we actually are
	 * @param expeIndicesJ : the column where we actually are
	 * @param precision    : Precision of measures
	 * @return a int value of k
	 */
	public int getkValue(int theoIndicesI, int expeIndicesJ, double precision) {
		double theoMassI = getTheoSpec().getMassList().get(theoIndicesI);
		double theoMassPrev = getTheoSpec().getMassList().get(theoIndicesI - 1);
		double expeMassJ = getExpeSpec().getMassList().get(expeIndicesJ);
		// We check all column from j to 0
		for (int k = expeIndicesJ; k >= 0; k--) {
			// if we found that j-k give an amino acid, we keep k value
			if (aaFound(theoIndicesI, expeIndicesJ, k, precision))
				return k;
			// if we pass the mass of the theoretical amino acid, we stop to not over
			// calculate
			else if ((expeMassJ - getExpeSpec().getMassList().get(k)) > (theoMassI - theoMassPrev + precision))
				return -1;
		}
		return -1; // a value of -1 to show that there is no k value
	}

	/**
	 * Method to get a realigned j value of where come from the realignment and get
	 * the associated realigned score
	 * 
	 * @param theoIndicesI    : The row where we are actually
	 * @param expeIndicesK    : The founded k indices
	 * @param expeIndicesJ    : The column where we are actually
	 * @param reAlignScore    : Score to add if it is needed to add offset to
	 *                        realign
	 * @param alignScoreToAdd : Score to add if there is an alignment with peak k
	 * @param precision       : Precision of measures
	 * @return an int[] where int[0] is the origin column (m) where come from the
	 *         re-alignment and int[1] is the calculated re-aligned score and int[3]
	 *         is the alignment type code
	 * 
	 */
	public int[] getBestRealignScore(int theoIndicesI, int expeIndicesK, int expeIndicesJ, int reAlignScore,
			int alignScoreToAdd, double precision) {
		int bestScore = -10000;
		int origin = -1;
		int[] result = new int[3];

		// m is a j value between 0 and k where a realign can be do if we accept mass
		// offset
		for (int m = expeIndicesK; m > -1; m--) {
			// the >= here is for keep the highest S value if there is multiple S with the
			// same best score
			if (getMatScore().getData(theoIndicesI - 1, m) > bestScore) {
				bestScore = getMatScore().getData(theoIndicesI - 1, m);
				origin = m;
			}
		}

		result[0] = origin;
		result[1] = bestScore + reAlignScore;
		result[2] = 1;

		if (origin == -1)
			return result;

		// We check for the last alignment if we have chain of Non Alignment to compare
		// the last mass offset found
		int lastAlign = theoIndicesI - 1;
		for (int l = theoIndicesI - 1; l > 0; l--) {
			if (getMatOrigin().getDataAlign(l, origin) != 0) {
				lastAlign = l;
				break;
			}
		}

		// if the difference of mass offset between actual state and last align (or
		// realign) is null, we consider that to an alignment
		if ((lastAlign != (theoIndicesI - 1)) && (theoIndicesI > 1)
				&& (Math.abs(getMatE().getData(theoIndicesI - 1, expeIndicesK)
						- getMatE().getData(lastAlign, origin)) < precision)) {
			result[1] = bestScore + alignScoreToAdd;
			result[2] = 2;
		}

		// we return the origin (value of m) and the associate calculated score and the
		// type of Alignment
		return result;

	}

	/**
	 * ==================================================================
	 * RETREATMENT OF THE HIT_MODIFIED AFTER ALIGNMENT
	 * ==================================================================
	 */

	/**
	 * Evaluate the best alignment with the maximum number of shared peaks when a
	 * modification is not-aligned of after a series of unfound amino acids. Several
	 * locations are tested
	 * 
	 * @param modified1            : The hit modified found after alignment
	 * @param experimentalMassList : Mass list of the experimental spectrum
	 * @param precision            : Precision of the Mass spectrometer
	 * @param experimentalMass     : Mass of the experimental peptide giver by the
	 *                             precursor
	 * @param peptideMass          : Mass of the OMS given peptide
	 * @return
	 */
	public double evaluate(String modified1, ArrayList<Double> experimentalMassList, double precision,
			double experimentalMass, double peptideMass) {

		String modified2 = modified1;
		String unExplained = "";
		double massBetterModified;
		int nbSharedPeaksAfterAlign;

		// First, we try to explain not-aligned block of modification, and we evaluate
		// the actual quality of alignment
		int index = modified1.indexOf("_");
		if (index != -1) {
			modified1 = modified1.substring(0, index);
			unExplained = modified2.substring(index + 1);

			// calculate the number of peaks without the not-aligned mass left at the end
			// to evaluate if the unexplained mass is really not-aligned or on the last
			// amino acids
			massBetterModified = calculateModifiedPepMass(modified1, precision);
			nbSharedPeaksAfterAlign = numberSharedPeaks(modified1, experimentalMassList, precision);

			// calculate the number of shared peaks with the deltaM on the last aa

			String modified3 = modified1.concat(unExplained);
			int nbAlignAtEnd = numberSharedPeaks(modified3, experimentalMassList, precision);

			if (nbAlignAtEnd > nbSharedPeaksAfterAlign) {
				nbSharedPeaksAfterAlign = nbAlignAtEnd;
				massBetterModified = calculateModifiedPepMass(modified3, precision);
				modified1 = modified3;
			}

		} else {
			// initialize the number of shared masses if there is no not-aligned mass
			massBetterModified = calculateModifiedPepMass(modified1, precision);
			nbSharedPeaksAfterAlign = numberSharedPeaks(modified1, experimentalMassList, precision);
		}

		setNbPeakAfterAlign(nbSharedPeaksAfterAlign);

		// The not-aligned mass is now explain, or removed from modified1

		String bestModified = modified1;

		// Try to remove complementary mass offset 
		// Not done if the number of peaks is increased by more than two peaks
		// If the number of peaks > +2, it is considered as an information that must be kept in the alignment
		// Could be modified...the information is in the preAligned column
		
		modified1 = eliminateComplementaryDelta(modified1, precision);
		if (numberSharedPeaks(modified1, experimentalMassList, precision) < nbSharedPeaksAfterAlign + 2)
			bestModified = modified1;
		else
			modified1 = bestModified;

		// setModifiedAfterBestScore(modified1);

		// if bestModified contains negative offSet, try to remove amino acids to
		// explain it
		modified1 = eliminateNegativeOffset(modified1, precision);

		if (numberSharedPeaks(modified1, experimentalMassList, precision) >= nbSharedPeaksAfterAlign)
			bestModified = modified1;
		else
			modified1 = bestModified;

		// We'll try to evaluate for each deltaM if it can be a neutral loss or not
		// For this, we compare the number of shared peaks with and without the shift

		Pattern p = Pattern.compile("-?\\d+,\\d+");
		Matcher m = p.matcher(modified1);
		int beginIndex = 0;

		// If there is no staying offset, we just have to recalculate parameters like
		// shared peaks, mass of the last modified

		if (!m.find()) {
			nbSharedPeaksAfterAlign = numberSharedPeaks(modified1, experimentalMassList, precision);
			massBetterModified = calculateModifiedPepMass(modified1, precision);
			bestModified = modified1;
		} else {
			// For each deltaM
			while (beginIndex < bestModified.length()) {
				if (m.find(beginIndex)) {
					beginIndex = m.start();
					int endIndex = modified1.indexOf(']', beginIndex) + 1;

					// there is an offset at beginIndex position, but where is it best if there are
					// several unfound amino acids just before?

					int bestPosition = beginIndex - 7;

					while (bestPosition > 7) {

						modified2 = modified1;

						if ((modified2.charAt(bestPosition) == '[') && (modified2.charAt(bestPosition + 2) == ']')) {
							String modified3 = modified1.substring(0, bestPosition + 3);
							modified3 = modified3.concat(modified1.substring(beginIndex - 1, endIndex));
							modified3 = modified3.concat(modified1.substring(bestPosition + 3, beginIndex - 1));
							modified3 = modified3.concat(modified1.substring(endIndex));
							int nbPeaksWithModif = numberSharedPeaks(modified3, experimentalMassList, precision);
							if (nbPeaksWithModif > nbSharedPeaksAfterAlign) {
								nbSharedPeaksAfterAlign = nbPeaksWithModif;
								massBetterModified = calculateModifiedPepMass(modified3, precision);
								bestModified = modified3;

							}
							bestPosition -= 3;
						} else
							break;
					}
					// try in the case we remove the offset
					modified2 = modified1.substring(0, beginIndex - 1).concat(modified1.substring(endIndex));
					int nbPeaksWithModifNeutre = numberSharedPeaks(modified2, experimentalMassList, precision);

					if (nbPeaksWithModifNeutre >= nbSharedPeaksAfterAlign) {
						nbSharedPeaksAfterAlign = nbPeaksWithModifNeutre;
						massBetterModified = calculateModifiedPepMass(modified2, precision);
						bestModified = modified2;
						modified1 = modified2;
						m = p.matcher(modified1);
					} else
						beginIndex = endIndex;
				} else
					break;
			}
		}

		double offset = getExpeSpec().getMainMass() - massBetterModified;

		if (offset < -precision) {
			bestModified = tryToCumulateOffSets(bestModified, offset, nbSharedPeaksAfterAlign, experimentalMassList,
					precision);
			nbSharedPeaksAfterAlign = numberSharedPeaks(bestModified, experimentalMassList, precision);
			massBetterModified = calculateModifiedPepMass(bestModified, precision);
			offset = getExpeSpec().getMainMass() - massBetterModified;
		}

		if (Math.abs(offset) > precision) {
			bestModified = bestModified + "_[";
			bestModified = bestModified + SGXProperties.DECIMALFORM.format(offset);
			bestModified = bestModified + "]";
		}

		setNbReajustedPeak(nbSharedPeaksAfterAlign);
		setModifiedAfterBestScore(bestModified);

		return offset;
	}

	/**
	 * Calculate the number of peak that are shared (aligned) between experimental
	 * spectrum and a generated spectrum from the hitModified sequence
	 * 
	 * @param modified             : The sequence of the hitModified
	 * @param experimentalMassList : The list of mass peak from experimental spetrum
	 * @param precision            : Precision of the mass spectrometer
	 * @return the number of shared peak between hitModified generated spectrum and
	 *         experimental spectrum
	 */
	public int numberSharedPeaks(String modified, ArrayList<Double> experimentalMassList, double precision) {

		//TODO: should be optimized!!
		
		double variablePrecision = precision;

		ArrayList<Double> modifiedPeaks = generatePeaks(modified, precision);
		if (modifiedPeaks == null)
			return 0;

/*		System.out.println("Theo List : "+modified);
		for (double e : modifiedPeaks)
			System.out.println(e);
		System.out.println("Expe List :");
		for (double e : experimentalMassList)
			System.out.println(e);
*/
		int nbSharedPeaks = 0;
		double prevPeak = 0.0;

		Collections.sort(experimentalMassList);
		variablePrecision += 0.005; // value depending on the precision of the amino acids

		for (Double modPeak : modifiedPeaks) {
			
			if (Math.abs(modPeak - prevPeak) > variablePrecision) {

				for (Double expePeak : experimentalMassList) {
					if (Math.abs(modPeak - expePeak) < variablePrecision) {
						nbSharedPeaks += 1;
						prevPeak = modPeak;
						// System.out.println(modPeak + " ::::: " + expePeak);
						break;
					} else if (expePeak > modPeak)
						break;

				}
			}
		}

		return nbSharedPeaks;
	}

	/**
	 * Calculate the mass of the proposed modified
	 * 
	 * @param modified  : Sequence of hit modified
	 * @param precision : the precision of the instrument
	 * @return
	 */
	public double calculateModifiedPepMass(String modified, double precision) {

		// First, we try to remove the not-aligned mass if any
		
		double totalMass = 0.0;		
		
		int index = modified.indexOf("_");
		if (index != -1) 
			modified = modified.substring(0, index);
		

		// Second we add shifts
		Pattern p = Pattern.compile("-?\\d+,\\d+");
		Matcher m = p.matcher(modified);
		
		while (m.find()) {		
			int beginIndex = m.start();
			int endIndex = modified.indexOf(']', beginIndex) + 1;
			String mass = modified.substring(beginIndex,endIndex-1);
			// French numeric notation
			mass = mass.replace(",", ".");
			double massValue = Double.valueOf(mass);
			totalMass += massValue;	
		}
			
        // Remove brackets
			
		modified = modified.replaceAll("\\[", "");		
		modified = modified.replaceAll("\\]", "");
		int i = 0;
		double massPep=0;
		
		while (i < modified.length()) {
			double mass = AminoAcids.getUnitMass(String.valueOf(modified.charAt(i)));
			if (mass != -1)
				massPep+=mass;
			i += 1;
		}
		
		totalMass += massPep + AminoAcids.getUnitMass("NT") + AminoAcids.getUnitMass("CT");
		return totalMass;
	}

	/**
	 * Generate mass peak list to allow the comparison between hitModified spectrum
	 * and experimental spectrum
	 * 
	 * @param hitModified : sequence hitModified
	 * @param precision   : the precision of the instrument
	 * @return a list of peak mass
	 */
	private static ArrayList<Double> generatePeaks(String hitModified, double precision) {

		int length;

		int existNotExplained = hitModified.indexOf('_');
		if (existNotExplained == -1)
			length = hitModified.length();
		else
			length = existNotExplained;

		double bSum = AminoAcids.getUnitMass("H+");
		double ySum = AminoAcids.getUnitMass("CT") + AminoAcids.getUnitMass("NT") + AminoAcids.getUnitMass("H+");
		ArrayList<Double> peaks = new ArrayList<>();

		int debIndex = 0;
		double valueOffset = 0.0;

		for (int i = 0; i < length; i++) {

			if (hitModified.charAt(i) == '[') {
				debIndex = i;

				while (hitModified.charAt(i) != ']') {
					i++;
					if (i >= hitModified.length()) {
						System.out.println("Problem issue with hitModified :" + hitModified);
						return new ArrayList<>();
					}
				}

				// even if the amino acid is not present, we add the peak
				if (i == debIndex + 2) {

					bSum += AminoAcids.getUnitMass(hitModified.charAt(debIndex + 1));

					peaks.add(bSum);

				}

				else {

					String number = hitModified.substring(debIndex + 1, i);
					number = number.replace(',', '.');
					valueOffset = Double.valueOf(number);

					bSum += valueOffset;

					if (!peaks.isEmpty()) {
						peaks.add(bSum);
					} else {
						peaks.add(valueOffset);
					}

				}
			} else {

				int lastIndex = 0;
				// We add the 2 peaks that frame the amino acid if necessary
				if (!peaks.isEmpty()) {
					lastIndex = peaks.size() - 1;
					if (Math.abs(peaks.get(lastIndex) - bSum) > precision) {
						peaks.add(bSum);
					}
				} else {
					peaks.add(bSum);
				}

				bSum = bSum + AminoAcids.getUnitMass(hitModified.charAt(i));
				peaks.add(bSum);

			}
			// when we reach this point, modified.charAt(i) is an amino acid
		}

		int nbPeaks = peaks.size();
		double massTot = bSum + ySum;

		double y = massTot - AminoAcids.getUnitMass("H+");
		peaks.add(y);
		for (int i = 0; i < nbPeaks - 1; i++) {
			y = massTot - peaks.get(i);
			peaks.add(y);
		}

		Collections.sort(peaks);

		// elimination of redundant peaks and peaks that are out of range (negative or
		// more than total mass)
		double prevPeak = peaks.get(0);
		double actualPeak;
		for (int i = 1; i < peaks.size(); i++) {
			actualPeak = peaks.get(i);
			if (Math.abs(actualPeak - prevPeak) < precision) {
				peaks.remove(i);
				i--;
			}
			prevPeak = actualPeak;
		}

		return (peaks);
	}

	/**
	 * Use to delete offset when there is two opposite mass in the hit modified
	 * example : [150.2]SDS[-150.2]KR --> SDSKR
	 * 
	 * @param modified1 : The hit modified sequence
	 * @param precision : The mass spectrometer precision
	 * @return a new hit modified without complementary offset mass
	 */
	private String eliminateComplementaryDelta(String modified1, double precision) {
		// System.out.println("Elimination des complementaires");

		int index = modified1.indexOf("_");
		if (index == -1)
			index = modified1.length();

		Pattern p = Pattern.compile("-?\\d+,\\d+");
		Matcher m = p.matcher(modified1);
		ArrayList<Double> offSetInModified = new ArrayList<Double>();
		ArrayList<Double> offSetToRemove = new ArrayList<Double>();
		double offSet;

		while (m.find()) {
			int beginIndex = m.start();
			if (beginIndex < index) {

				int endIndex = modified1.indexOf(']', beginIndex);
				String value = modified1.substring(beginIndex, endIndex);
				value = value.replace(',', '.');
				offSet = Double.valueOf(value);
				offSetInModified.add(offSet);

				for (int i = 0; i < offSetInModified.size(); i++) {
					if (Math.abs(offSet + offSetInModified.get(i)) < precision) {
						offSetToRemove.add(offSetInModified.get(i));
						offSetToRemove.add(offSet);
						break;
					}
				}
			}
		}

		while (!offSetToRemove.isEmpty()) {
			offSet = offSetToRemove.get(0);
			String toRemove = "[".concat(String.valueOf(offSet));
			toRemove = toRemove.replace('.', ',');
			int begin = modified1.indexOf(toRemove);
			int end = modified1.indexOf("]", begin);
			if (begin > -1) {
				String modified = modified1;
				modified1 = modified1.substring(0, begin);
				modified1 = modified1.concat(modified.substring(end + 1));
			}
			offSetToRemove.remove(0);
		}

		if (modified1.charAt(modified1.length() - 1) == '_')
			modified1 = modified1.substring(0, modified1.length() - 1);

		return modified1;
	}

	/**
	 * Use to try to explain negative offset with deletion of amino acids
	 * 
	 * @param modified1 : The hit modified sequence
	 * @param precision : The mass spectrometer precision
	 * @return the new hit modified with some negative offset explain
	 */
	public String eliminateNegativeOffset(String modified1, double precision) {
		// System.out.println("Elimination des neg offsets");

		Pattern p = Pattern.compile("-\\d+,\\d+");
		Matcher m = p.matcher(modified1);
		double offSet = 0.0;
		int beginIndex = 0;
		int i = 0;

		while (m.find(beginIndex + 1)) { // do not process several times the same offSet
			beginIndex = m.start();
			int endIndex = modified1.indexOf(']', beginIndex);
			String value = modified1.substring(beginIndex, endIndex);
			value = value.replace(',', '.');
			offSet = -Double.valueOf(value);
			i = beginIndex;

			// While offSet negative and aa not_found, try to remove aminoacids
			// starting position i, find the first unfound aa

			while (i > 3) {
				if ((modified1.charAt(i - 2) == ']') && (modified1.charAt(i - 4) == '['))
					i = i - 3;
				else
					break;
			}

			int startUnSeen = i;
			while ((offSet > precision) && (i < beginIndex)) {
				// System.out.println("Pour le modified : " + modified1);
				// System.out.println("a i = " + i);
				// System.out.println(offSet);
				String aa = String.valueOf(modified1.charAt(i));
				double mass = AminoAcids.getUnitMass(aa);
				// System.out.println(aa + " of mass " + mass);
				if (offSet - mass > -precision) {
					offSet = offSet - mass;
					i = i + 3;
					setExplainedMass(getExplainedMass() + mass);
					// System.out.println("I cut a neg offset");
				} else {
					break;
				}
			}

			String modified = modified1;
			if (startUnSeen > 1)
				modified1 = modified1.substring(0, startUnSeen - 1);
			else
				modified1 = "";

			modified1 = modified1.concat(modified.substring(i - 1, beginIndex - 1));
			beginIndex = startUnSeen;

			if (Math.abs(offSet) > precision) {
				String addOffset = "[#";
				addOffset = addOffset.concat(SGXProperties.DECIMALFORM.format(offSet));
				addOffset = addOffset.concat("]");
				modified1 = modified1.concat(addOffset);
				beginIndex = beginIndex + addOffset.length();
			}

			if (endIndex + 1 < modified.length())
				modified1 = modified1.concat(modified.substring(endIndex + 1));

			m = p.matcher(modified1);

			if (beginIndex >= modified1.length())
				break;
		}

		return modified1.replace("#", "-");
	}

	/**
	 * This method is run only if there is not-aligned offSets. It may be explained
	 * by an inadequate alignment at the start because we give some advance to
	 * non-tryptic peptides. But, at the end of the alignment, we can consider it
	 * was not the best choice Trying to improve this alignment by adding
	 * not-aligned offSet with another one
	 * 
	 * @param bestModified            : Actual best hit modified
	 * @param offSet                  : The not-aligned staying mass
	 * @param nbSharedPeaksAfterAlign : Number of shared peak with spectra for the
	 *                                best hit modified
	 * @param experimentalMassList    : list of masses in experimental spectrum
	 * @param precision               : The precision of the Mass spectrometer
	 * @return a new hit modified
	 */
	public String tryToCumulateOffSets(String bestModified, double offSet, int nbSharedPeaksAfterAlign,
			ArrayList<Double> experimentalMassList, double precision) {

		Pattern p = Pattern.compile("\\[-?\\d+,\\d+");
		Matcher m = p.matcher(bestModified);
		int beginIndex = 0;
		int endIndex = 0;

		while (m.find(endIndex)) {
			beginIndex = m.start();
			endIndex = bestModified.indexOf(']', beginIndex);
			String value = bestModified.substring(beginIndex + 1, endIndex);
			value = value.replace(',', '.');
			double testOffSet = offSet + Double.valueOf(value);

			String modified1 = bestModified.substring(0, beginIndex + 1);
			modified1 = modified1.concat(String.format("%.3f", testOffSet));
			modified1 = modified1.concat(bestModified.substring(endIndex, bestModified.length()));

			int nbPeaksAndMass = numberSharedPeaks(modified1, experimentalMassList, precision);

			if (nbPeaksAndMass >= nbSharedPeaksAfterAlign) {
				nbSharedPeaksAfterAlign = nbPeaksAndMass;
				bestModified = modified1;
				m = p.matcher(bestModified);
			}
			if (endIndex >= bestModified.length())
				break;
		}
		return bestModified;
	}

	/**
	 * Use to remove brackets around not found Amino acids when the B or Y peak
	 * corresponding to the amino acid in the sequence is in the experimental peak
	 * list
	 *
	 * @param betterModified       : The better sequence of hit Modified
	 * @param experimentalMassList : The list of mass that are in the experimental
	 *                             spectrum
	 * @param precision            : The precision of the mass spectrometer
	 * @return the new sequence without brackets around found amino acids after all
	 *         post treatments
	 */
	public static String removeBracketsForAlignedAA(String betterModified, ArrayList<Double> experimentalMassList,
			double precision) {

		StringBuilder workModified = new StringBuilder();
		// we put the first char of the sequence to initiate
		workModified.append(betterModified.charAt(0));

		for (int i = 1; i < betterModified.length(); i++) {
			boolean founded = false;

			if (betterModified.charAt(i - 1) == '[' && betterModified.charAt(i + 1) == ']') {
				// cut the sequence to get easily b and y for the char at i
				String partB = betterModified.substring(0, i + 2);
				String partY = betterModified.substring(i - 1);
				// generate peakList to get the good peak from actual inside bracket amino acid
				// and check if the B or the Y is in the spectrum
				ArrayList<Double> peakListB = generatePeaks(partB, precision);
				double peakB = peakListB.get(peakListB.size() - 2);
				ArrayList<Double> peakListY = generatePeaks(partY, precision);
				double peakY = peakListY.get(peakListY.size() - 1);
				if (SGXProperties.DEBUG_MODE) {
					// System.out.println(betterModified.charAt(i) + " : " + partB + " - " + partY +
					// " == " + peakB
					// + " and " + peakY);
				}

				// we check if the peak is in the experimental peak list
				for (double mass : experimentalMassList) {
					// if the peak is in the experimental spectrum
					if (Math.abs(mass - peakB) < precision || Math.abs(mass - peakY) < precision) {
						workModified.deleteCharAt(workModified.length() - 1);
						workModified.append(betterModified.charAt(i));
						founded = true;
						i++;
						break;
					}
				}
				if (!founded) {
					workModified.append(betterModified.charAt(i));
				}
			} else {
				workModified.append(betterModified.charAt(i));
			}

		}

		return workModified.toString();
	}

	/**
	 * 
	 * @param betterModified       : The better sequence of hit Modified
	 * @param experimentalMassList : The list of mass that are in the experimental
	 *                             spectrum
	 * @param precision            : The precision of the mass spectrometer
	 * @return the rate of the alignment
	 */
	public double calculateConfidenceRate(String betterModified, ArrayList<Double> experimentalMassList,
			double precision) {

		int nbSharedPeaks = numberSharedPeaks(betterModified, experimentalMassList, precision);
		int nbTheoPeaks = generatePeaks(betterModified, precision).size();

		return (nbSharedPeaks * 1.0 / nbTheoPeaks * 1.0);

	}

	/**
	 * 
	 * @param betterModified
	 * @param experimentalMassList
	 * @param precision
	 * @return the rate of intensity explain on total intensity of experimental
	 *         peaks
	 */
	public double explainedIntensityRate(String betterModified, Map<Double, Double> experimentalPeakList,
			double precision) {

		double variablePrecision = precision;

		double totalIntensity = 0.0;
		double explainedIntensity = 0.0;

		ArrayList<Double> experimentalMassList = new ArrayList<>(experimentalPeakList.keySet());
		// add all intensities from the peak list in the variable
		for (Double mass : experimentalMassList) {
			totalIntensity += experimentalPeakList.get(mass);
		}

		ArrayList<Double> modifiedPeaks = generatePeaks(betterModified, precision);

		double prevPeak = 0.0;

		Collections.sort(experimentalMassList);

		boolean changePrecision = false;

		for (Double modPeak : modifiedPeaks) {
			// we add a variation to the precision due to a possible error increased by
			// amino acid mass addition
			if (!changePrecision && modPeak >= 1000.0) {
				changePrecision = true;
				variablePrecision += 0.05;
			}
			if (Math.abs(modPeak - prevPeak) > variablePrecision) {

				for (Double expePeak : experimentalMassList) {
					if (Math.abs(modPeak - expePeak) < variablePrecision) {
						explainedIntensity += experimentalPeakList.get(expePeak);
						prevPeak = modPeak;
						// System.out.println(modPeak + " ::::: " + expePeak);
						break;
					} else if (expePeak > modPeak)
						break;

				}
			}
		}

		return explainedIntensity / totalIntensity;
	}

	/**
	 * Use to put informations in Score and origin matrices
	 * 
	 * @param theoIndicesI  : the row coord to fill
	 * @param expeIndiciesJ : the column coord to fill
	 * @param score         : the score value to put in the score matrix
	 * @param origin        : the column origin to put in origin matrix
	 * @param alignType     : the type of alignment to put in origin matrix (-1 =
	 *                      NA, 0 = ReA, 1 = A)
	 */
	public void setMatricesData(int theoIndicesI, int expeIndiciesJ, int score, int origin, int alignType) {
		getMatScore().setData(theoIndicesI, expeIndiciesJ, score);
		getMatOrigin().setDataOrigin(theoIndicesI, expeIndiciesJ, origin);
		getMatOrigin().setDataAlign(theoIndicesI, expeIndiciesJ, alignType);
	}

	/**
	 * Show the 3 matrices : Score, Origin, E
	 */
	public void showMatrices() {
		System.out.println("\n" + getExpeSpec().getMassList());
		getMatScore().show(getTheoSpec().getMassList().size(), getExpeSpec().getMassList().size());
		System.out.println("");
		getMatOrigin().show(getTheoSpec().getMassList().size(), getExpeSpec().getMassList().size());
		System.out.println("");
		getMatE().show(getTheoSpec().getMassList().size(), getExpeSpec().getMassList().size());
	}

	/**
	 * Concat the result to put in the CSV result line WARNIG : Need to add : "Title
	 * ; PSM ;" and the final "\n"
	 */
	public void makeFinalResult() {
		setFinalResult(SGXProperties.DECIMALFORM.format(getTrueDeltaMass()) + ";" + getNbPeakInitial() + ";"
				+ getNbReajustedPeak() + ";" + getHitModifiedSeq() + ";" + getOtherModified() + ";"
				+ getModificationNumber() + ";" + SGXProperties.DECIMALFORM.format(getTotDeltaMass()) + ";"
				+ getMaxScore() + ";" + getConfidenceRate());
	}
	

	private String reWriteWithFixedModif(String pepModified) {

		/***
		 * We must rewrite amino acid with fixed modifications in the hitModified string
		 * to highlight these modifications just before returning the hitModified
		 ****/

		ArrayList<String> aaFixedModified = SGXProperties.getAAFixedModifications();
		for (String str : aaFixedModified) {
			String str1 = SGXProperties.getAAFixedModification(str);
			pepModified = pepModified.replaceAll(str, str1);
		}

		/*
		 * Due to fixed modifications, some modification masses must be cumulated
		 * 
		 */
		Pattern p = Pattern.compile("\\d+\\]\\[\\d+");
		Matcher matcher = p.matcher(pepModified);

		while (matcher.find()) {
			// concatenation of fixed and variable modifications et one place
			int indexDeb = matcher.start();
			int index = indexDeb;
			int indexEnd = matcher.end();
			while (pepModified.charAt(indexDeb) != '[') {
				indexDeb--;
				if (indexDeb < 0)
					return pepModified;
			}
			while (pepModified.charAt(index) != ']') {
				index++;
				if (index >= pepModified.length())
					return pepModified;
			}
			while (pepModified.charAt(indexEnd) != ']') {
				indexEnd++;
				if (indexEnd == pepModified.length())
					return pepModified;
			}

			String temp = pepModified.substring(indexDeb + 1, index).replaceAll(",", ".");
			double firstDelta = Double.parseDouble(temp);
			temp = pepModified.substring(index + 2, indexEnd).replaceAll(",", ".");
			double secondDelta = Double.parseDouble(temp);
			double delta = firstDelta + secondDelta;
			String pepModified1 = pepModified.substring(0, indexDeb);
			pepModified1 = pepModified1.concat("[").concat(String.valueOf(delta)).concat("]");
			pepModified1 = pepModified1.concat(pepModified.substring(indexEnd));
			pepModified = pepModified1;
			matcher = p.matcher(pepModified);
		}
		return pepModified;

	}

	// Getters and Setters

	public MatrixScore getMatScore() {
		return _matScore;
	}

	public void setMatScore(MatrixScore matScore) {
		_matScore = matScore;
	}

	public MatrixE getMatE() {
		return _matE;
	}

	public void setMatE(MatrixE matE) {
		_matE = matE;
	}

	public MatrixOrigin getMatOrigin() {
		return _matOrigin;
	}

	public void setMatOrigin(MatrixOrigin matOrigin) {
		_matOrigin = matOrigin;
	}

	public ExperimentalSpectrum getExpeSpec() {
		return _expeSpec;
	}

	public void setExpeSpec(ExperimentalSpectrum expeSpec) {
		_expeSpec = expeSpec;
	}

	public TheoreticalSpectrum getTheoSpec() {
		return _theoSpec;
	}

	public void setTheoSpec(TheoreticalSpectrum theoSpec) {
		_theoSpec = theoSpec;
	}

	public String getHitModifiedSeq() {
		return _hitModifiedSeq;
	}

	public void setHitModifiedSeq(String hitModifiedSeq) {
		_hitModifiedSeq = hitModifiedSeq;
	}

	public String getModifiedAfterBestScore() {
		return _modifiedAfterBestScore;
	}

	public void setModifiedAfterBestScore(String modifiedAfterBestScore) {
		_modifiedAfterBestScore = modifiedAfterBestScore;
	}

	public String getOtherModified() {
		return _otherModified;
	}

	public void setOtherModified(String otherModified) {
		_otherModified = otherModified;
	}

	public int getNbPeakInitial() {
		return _nbPeakInitial;
	}

	public void setNbPeakInitial(int nbPeakInitial) {
		_nbPeakInitial = nbPeakInitial;
	}

	public int getNbPeakAfterAlign() {
		return _nbPeakAfterAlign;
	}

	public void setNbPeakAfterAlign(int nbPeakAfterAlign) {
		_nbPeakAfterAlign = nbPeakAfterAlign;
	}

	public int getNbReajustedPeak() {
		return _nbReajustedPeak;
	}

	public void setNbReajustedPeak(int nbReajustedPeak) {
		_nbReajustedPeak = nbReajustedPeak;
	}

	public double getTotDeltaMass() {
		return _totDeltaMass;
	}

	public void setTotDeltaMass(double totDeltaMass) {
		_totDeltaMass = totDeltaMass;
	}

	public double getTrueDeltaMass() {
		return _trueDeltaMass;
	}

	public void setTrueDeltaMass(double trueDeltaMass) {
		_trueDeltaMass = trueDeltaMass;
	}

	public double getExplainedMass() {
		return _explainedMass;
	}

	public void setExplainedMass(double explainedMass) {
		_explainedMass = explainedMass;
	}

	public int getMaxScore() {
		return _maxScore;
	}

	public void setMaxScore(int maxScore) {
		_maxScore = maxScore;
	}

	public int getModificationNumber() {
		return _modificationNumber;
	}

	public void setModificationNumber(int modificationNumber) {
		_modificationNumber = modificationNumber;
	}

	public double getConfidenceRate() {
		return _confidenceRate;
	}

	public void setConfidenceRate(double confidenceRate) {
		_confidenceRate = confidenceRate;
	}

	public String getFinalResult() {
		return _finalResult;
	}

	public void setFinalResult(String finalResult) {
		_finalResult = finalResult;
	}

	public static int getScoreToApply(String alignType) {
		return SGXProperties.SCORETOAPPLY.get(alignType);
	}

}
