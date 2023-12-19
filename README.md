# SpecGlobX

## Presentation

Welcome to the SpecGlobX project !

SpecGlobX is a java program that can be used when interpreting mass spectra in proteomics. SpecGlobX aligns peptides and experimental mass spectra even though several shifts are required to perform good alignments. SpecGlobX generates a theoretical spectrum for each peptide associated with each PSM given as input and aligns this theoretical spectrum with its associated experimental spectrum, capturing dissimilarities arising from multiple modifications.<br> SpecGlobX supports spectra in mgf or mzML format and a comma-delimited spreadsheet must describe the list of PSMs. SpecGlobX returns the best alignment for each PSM, splitting the mass difference between spectra and peptides into one or several shifts.


## How to compile ?

Java VM is needed. The main of SpecGlobX is declared in SpecGlobXGUI class.
It will launch an instance of SpectralAlignement which contains the alignment algorithm called SpecGlobX.

 java8 JDK and more recent are required to compile.

We tested with success the compilation with java 17 JDK in changing some parameters in the *pom.xml*<br>
Change  
		``<source>1.8</source>``
        ``<target>1.8</target>``<br>
to
      ``<release>17</release>``

The program was tested with Eclipse java 17 version openJDK from adoptium and JRE8 oracle

Tested version :

* jre8
* jdk-17.0.1
* jdk-17.0.2

### Use SpecGlobX through its GUI Mode

A jar file can be downloaded directly from the target repository. Double-click on the jar file.<br><br>
A window appears, where you can select input and output files. The most important parameters can be chosen via the GUI. Other parameters are stored in the *config.properties* file that **must be in the same folder as the jar**. You can modify the parameters in this file before and while the jar is running.<br><br>
The spectra file can be in the MGF or mzML format.<br><br>
The CSV input file contains the list of PSMs (Peptide-to-Spectrum Matches) (obtained by any OMS or other).<br>
**BE AWARE: The CSV separator must be ";". The first row is not read (considered as header) **<br>

Choose the column where SpecGlobX can find "Spectrum Name" (identification of spectra) and "Peptide" (select letter or number of the excel column for example)<br>
**BE AWARE: It must be the title of the scan for MGF data file and the index (nb scan) for mzML**<br><br>
Select folder and give a name ending by ".csv" for the output file.<br>
Select your parameters :<br>

* Multiprocessing : activate multiprocessing to increase speed, and choose the number of threads to use.
* Peaks filter : choose the way to filter with "Most intense" (keep the n most intense peaks in each experimental spectrum) or "Max intensity percentage" (keep all peaks with more than n% of the maximum intensity).
* Mass spectrometer precision : choose the value of precision you need in Dalton

Click on the "Launch Alignments" button to launch the process.<br>
Wait until the sound and the message in log "=====PROCESS COMPLETED=====".
Logs are available at the end.

** If you encounter some issue, verify at first that titles are in the same format in the input csv file and in the spectra file.**


### Command Mode

The jar file can be launch in command mode with the option ```--c```.

Exemple :

``java -jar SpecGlobPub-1.0.0.jar --c -msfile fichierScan.mgf -csvfile data_information.csv -outfile outputResults.csv``

This command launches SpecGlobX in command mode (useful for pipeline implementations).
In this example, column for ScanID(Title) and peptide sequence(PSM) are set in the config file.<br>

Parameters :

```
* --c : enable command mode
* -msfile [] : Spectra file (absolute path preferred)
* -csvfile [] : CSV file that contain at least PSMs (spectrum, peptide)
* -outfile [] : the output file
* -titlecol [] : identification (a-Z) of column that contains titles of scans (or id if MzML) *(optional)*
* -pepcol [] : identification (a-Z) of column that contain the peptide sequences *(optional)*
```

## Configuration



The file **config.properties** contains the parameters configuration of SpecGlobX and it must be placed in same folder as the jar.

Parameters Description :

```
* parallelize : alignment launched in parallel mode or not (true or false) *(Is in GUI)*
* nbthread : Number of threads to launch for parallelization (be cautious about computer performances)(in GUI)
* precision : precision of fragmented ions provided by the mass spectrometer. Default value of 0.02 *(Is in GUI
* decimalFormat : Number of decimal written in results for masses. Default value = 4
* scoreMinDisplay : Minimum alignment score above which a result is returned. It is important to note that scores can be negative.
* filter : filter type applied on spectrum peaks (0 for intensity rate and 1 for number of maximal intensity peaks). Default = 1
  * peakIntensityRate : Minimal % of the best intensity used to filter peaks. Default = 1
  * peakNumberKeeped : Number of maximal intense peak to keep. Default = 60

* modif : Fixed modification on amino acids if required. Default 57.021464@C
* debug : Enable debug mode (matrices and information about alignments are written. Not recommended for non-experienced user)


* Dynamic programming scores : Modification of the scores  is not recommended for non-experienced user.

  * scoreNonAlign : Penalty for missing amino acid. Default = -4
  
  * scoreReAlignNative :  Penalty for a mass shift to align an original peak. Default = -8
  * scoreReAlignSym : Penalty for a mass shift to align a complementary peak. Default = -8
  * scoreReAlignBoth : Penalty for a mass shift to align a merged peak (original + complementary). Default = -6
  
  * scoreAlignNative : Bonus for an aligned original peak. Default = 7
  * scoreAlignSym : Bonus for an aligned complementary peak. Default = 7
  * scoreAlignBoth : Bonus for an aligned merged peak (original + complementary). Default = 10
  
  * scoreReAlignNativeNO : Bonus for a realignement with deltaM=0 (missing peak) with an original peak. Default = 2
  * scoreReAlignSymNO : Bonus for a realignement with deltaM=0 with a complementary peak. Default = 2
  * scoreReAlignBothN : Bonus for a realignement with deltaM=0 with a merged peak. Default = 5
  ```

## Results

Results are returned under the CSV format with one line per alignment. <br>


Detailed description of each column :

```
* Title = the spectrum title
* Peptide = Sequence given by OMS method (the giver peptide)
* MassDelta = mass difference observed between peptide and spectrum precursor
* SharedPeaksBeforeAlign = initial number of shared peaks between the theoretical spectrum of the peptide and the experimental spectrum
* SharedPeaksAfterAlign = number of shared peaks between experimental and theoretical spectrum after alignment by SpecGlobX
* PreAlignedPeptide = alignment with mass offset provided by SpecGlobX
* AlignedPeptide = alignment after the post-processing step
* NbShift = number of modifications returned by SpecGlobX
* NotAlignedMass = not-aligned mass remaining after the post-processing step
* ScoreAlign = alignment score computed by the dynamic programming step
* IntensityExplained = sum of intensity of aligned peaks / sum of intensity of the spectrum computed on the aligned peptide
```

SpecGlobX uses a specific syntax to express alignments as strings in the AlignedPeptide column. The aim is to summarize information about the alignment, providing a simplified fragmentation summary, highlighting stretches of detected (resp. unfound)  amino acids in the alignment.

The alignment is done with the filtered experimental spectrum (not all the peaks are considered).

* When both peaks of an amino acid are not used in the alignment, the amino acid is written between brackets
* When the peaks of an amino acid are used in the alignment, but this alignment requires a shift, then the amino acid is written, preceded by the value of the mass shift between brackets
* When the peaks of an amino acid are used in the alignment and no shift is required for this alignment, then the amino acid is reported as such.
At the end of the string, the AlignedPeptide column indicates the "not-aligned mass" after the underscore character.

The preAlignedPeptide column supports the same syntax as the alignedPeptide column, but this column gives information on how SpecGlobX managed the alignment (the dynamic programming part with the traceBack) before the post-processing step.

