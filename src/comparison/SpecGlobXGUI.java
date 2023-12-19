package comparison;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JLayeredPane;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingWorker;
import javax.swing.WindowConstants;
import javax.swing.border.EtchedBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.DefaultCaret;

import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import utility.SGXProperties;

/**
 * Main Class of program that implements the GUI
 * 
 * @author Gregoire Prunier, Albane Lysiak, Dominique Tessier
 *
 */
public class SpecGlobXGUI {

	private JFrame frmSpecglobtool;
	private JTextField textFieldInputMGF;
	private JTextField textFieldInputCSV;
	private JTextField textFieldOutput;
	private JTextField textFieldValueFilter;
	private JTextField txtA;
	private JTextField txtB;
	private JTextField textFieldPRECISION;

	private File spectraFile = new File("spectra.mgf");
	private String dataType = "MGF";

	private String spectraPath;

	private String titleScanCol;
	private String peptideCol;
	private String csvInputPath;
	private String csvOutputPath;
	private byte filterID = 0;
	private String filterChoosed;
	private boolean isMultiprocessed;
	private byte nbThread = 1;
	private int valueFilter;
	private double precision;

	private JButton buttonLauch;

	public static JTextArea LOG;

	public static JProgressBar progressBar;

	private Task task;

	private Color colorBIBS1 = new Color(102, 193, 191);
	private Color colorBIBS2 = new Color(0, 140, 142);
	private Color colorBIBS3 = new Color(39, 86, 98);
	public static boolean commandMode = false;

	class Task extends SwingWorker<Void, Void> {
		/*
		 * Main task. should be executed in background thread....
		 */
		@Override
		public Void doInBackground() {

			progressBar.setValue(100);
			setProgress(0);

			try {

				SpecGlobX executionSGT = new SpecGlobX(spectraFile, csvInputPath, titleScanCol, peptideCol,
						csvOutputPath, dataType, isMultiprocessed, nbThread, filterID, valueFilter, precision);
				progressBar.setValue(400);
				(new Thread() {
					public void run() {
						LOG.append("=====STARTING ALIGNMENTS=====\n");
					}
				}).start();

				//TODO the parallelised and monotread methods should be merged to avoid
				// maintenance and test issues
				
				if (SGXProperties.IS_PARALLELIZED) {
					executionSGT.parallelAlignmentLaunch();
				} else {
					executionSGT.launchAlignments();
				}

				progressBar.setValue(1000);

				LOG.append("=====PROCESS COMPLETED=====\n");

				buttonLauch.setEnabled(true);
				// System.exit(0);

			} catch (JMzReaderException | FileNotFoundException | InterruptedException e1) {

				if (SpecGlobXGUI.commandMode)
					System.out.println("ERROR: Can't open result file \n");
				else
					SpecGlobXGUI.LOG.append("ERROR Can't open result file \n");
				e1.printStackTrace();

			}
			setProgress(1000);

			return null;
		}
	

		/*
		 * Executed in event dispatching thread
		 */
		@Override
		public void done() {
			Toolkit.getDefaultToolkit().beep();
			buttonLauch.setEnabled(true);
			LOG.append("Done!\n");
			setProgress(0);
		}
	}

	/**
	 * The main function to launch the application and initialize the interface
	 * 
	 * @param args
	 */
	public static void main(String[] args) {


		// checking --c argument
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--c")) {
				commandMode = true;
				break;
			}
		}

		// if we are in command mode, we take args and use them to create an instance of
		// SpecGlob
		if (commandMode) {

			System.out.println("=====PROCESS STARTED IN COMMAND MODE=====");
			try {
				SpecGlobX executionSGT = new SpecGlobX(args);

				if (SGXProperties.IS_PARALLELIZED) {
					executionSGT.parallelAlignmentLaunch();
				} else {
					executionSGT.launchAlignments();
				
				}

			} catch (JMzReaderException | InterruptedException | FileNotFoundException e) {
				// TODO Auto-generated catch block
				System.out.println("exception JMzReader or FileNotFound or Interrupted " + e);
				e.printStackTrace();
			}

			System.out.println("=====PROCESS COMPLETED=====");

			System.exit(0);

			// if it's not in command mode, we simply launch the interface
		} else {
			EventQueue.invokeLater(new Runnable() {
				public void run() {
					try {
						SpecGlobXGUI window = new SpecGlobXGUI();
						window.frmSpecglobtool.setVisible(true);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			});
		}

	}

	/**
	 * Invoked when task's progress property changes.
	 */
	public void propertyChange(PropertyChangeEvent evt) {
		if ("progress" == evt.getPropertyName()) {
			int progress = (Integer) evt.getNewValue();
			
			progressBar.setValue(progress);
			LOG.append(String.format("Completed %d%% of task.\n", task.getProgress()));
		}
	}

	/**
	 * Create the application.
	 * 
	 * @wbp.parser.entryPoint
	 */
	public SpecGlobXGUI() {
		initialize();

	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {

		final String fontType = "Tahoma";
		final String baseFolderPath = System.getProperty("user.home") + "\\Documents";

		ImageIcon img = new ImageIcon(getClass().getResource("/ressources/IconSpecGlobTool.png"));
		frmSpecglobtool = new JFrame();
		frmSpecglobtool.setResizable(false);
		frmSpecglobtool.setTitle("SpecGlobX\r\n");
		frmSpecglobtool.setBounds(100, 100, 703, 580);
		frmSpecglobtool.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frmSpecglobtool.getContentPane().setLayout(null);
		frmSpecglobtool.setIconImage(img.getImage());

		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(10, 405, 669, 100);
		frmSpecglobtool.getContentPane().add(scrollPane);

		LOG = new JTextArea();
		LOG.setEditable(false);
		DefaultCaret caret = (DefaultCaret) LOG.getCaret();
		caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
		scrollPane.setViewportView(LOG);

		JLabel lblNewLabel = new JLabel("Welcome to SpecGlobX* !");
		lblNewLabel.setForeground(colorBIBS3);
		lblNewLabel.setFont(new Font("Tahoma", Font.BOLD, 16));
		lblNewLabel.setHorizontalAlignment(SwingConstants.CENTER);
		lblNewLabel.setBounds(150, 10, 384, 30);
		frmSpecglobtool.getContentPane().add(lblNewLabel);

		JLayeredPane layeredPane = new JLayeredPane();
		layeredPane.setForeground(Color.WHITE);
		layeredPane.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane.setBackground(colorBIBS2);
		layeredPane.setBounds(10, 47, 225, 239);
		frmSpecglobtool.getContentPane().add(layeredPane);

		textFieldInputMGF = new JTextField();
		textFieldInputMGF.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldInputMGF.setText("spectra.mgf");
		textFieldInputMGF.setBounds(10, 54, 161, 23);
		layeredPane.add(textFieldInputMGF);
		textFieldInputMGF.setColumns(10);

		// The button to select the input spectra file
		JButton buttonInputMGF = new JButton("Folder");
		buttonInputMGF.setBounds(181, 54, 23, 23);
		buttonInputMGF.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				System.out.println("Pressing button !!");
				JFrame frame = new JFrame("Load spectra scan data");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Load spectra scan data (MGF or MZML)");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("MGF Files", "mgf"));
				// fc.addChoosableFileFilter(new FileNameExtensionFilter("MZXML Files",
				// "mzxml"));
				fc.addChoosableFileFilter(new FileNameExtensionFilter("MZML Files", "mzML"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showOpenDialog(frame);

				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();

	
					spectraFile = selectedFile;
					textFieldInputMGF.setText(selectedFilePath);

					// if the input file is a MGF file :
					if (selectedFile.getAbsolutePath().endsWith(".mgf")) {

						dataType = "MGF";

						// if the input file is a MZML file :
					} else if (selectedFile.getAbsolutePath().endsWith(".mzML")) {

						dataType = "MZML";

					}
				}

			}
		});
		layeredPane.add(buttonInputMGF);

		textFieldInputCSV = new JTextField();
		textFieldInputCSV.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldInputCSV.setText("data.csv");
		textFieldInputCSV.setBounds(10, 114, 161, 23);
		layeredPane.add(textFieldInputCSV);
		textFieldInputCSV.setColumns(10);

		// The button to select the CSV input file
		JButton buttonInputCSV = new JButton("Folder");
		buttonInputCSV.setBounds(181, 114, 23, 23);
		buttonInputCSV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFrame frame = new JFrame("Load CSV info file");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Load CSV info file");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("CSV Files", "csv"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showOpenDialog(frame);
				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();
					textFieldInputCSV.setText(selectedFilePath);

				}
			}
		});
		layeredPane.add(buttonInputCSV);

		JLabel labelSptFile = new JLabel("Spectra file (.mgf) (.MzML)");
		labelSptFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelSptFile.setBounds(10, 32, 191, 23);
		layeredPane.add(labelSptFile);

		JLabel labelPSMFile = new JLabel("PSM list file (.csv)");
		labelPSMFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelPSMFile.setBounds(10, 92, 191, 23);
		layeredPane.add(labelPSMFile);

		JLabel labelINPUT = new JLabel("INPUT");
		labelINPUT.setBounds(10, 10, 205, 19);
		labelINPUT.setForeground(colorBIBS3);
		layeredPane.add(labelINPUT);
		labelINPUT.setHorizontalAlignment(SwingConstants.CENTER);
		labelINPUT.setFont(new Font("Tahoma", Font.BOLD, 14));

		JLabel labelColumnNb = new JLabel("Column number where to find :");
		labelColumnNb.setHorizontalAlignment(SwingConstants.CENTER);
		labelColumnNb.setFont(new Font(fontType, Font.PLAIN, 12));
		labelColumnNb.setBounds(10, 152, 194, 23);
		layeredPane.add(labelColumnNb);

		JLabel labelTITLE = new JLabel("Spectrum ID");
		labelTITLE.setHorizontalAlignment(SwingConstants.CENTER);
		labelTITLE.setFont(new Font(fontType, Font.PLAIN, 12));
		labelTITLE.setBounds(10, 174, 88, 23);
		layeredPane.add(labelTITLE);

		JLabel labelPEPTIDE = new JLabel("Peptide");
		labelPEPTIDE.setHorizontalAlignment(SwingConstants.CENTER);
		labelPEPTIDE.setFont(new Font(fontType, Font.PLAIN, 12));
		labelPEPTIDE.setBounds(130, 174, 57, 23);
		layeredPane.add(labelPEPTIDE);

		txtA = new JTextField();
		txtA.setFont(new Font(fontType, Font.PLAIN, 12));
		txtA.setHorizontalAlignment(SwingConstants.CENTER);
		txtA.setText("1");
		txtA.setColumns(10);
		txtA.setBounds(25, 194, 57, 23);
		layeredPane.add(txtA);

		txtB = new JTextField();
		txtB.setFont(new Font(fontType, Font.PLAIN, 12));
		txtB.setText("2");
		txtB.setHorizontalAlignment(SwingConstants.CENTER);
		txtB.setColumns(10);
		txtB.setBounds(130, 194, 57, 23);
		layeredPane.add(txtB);

		JLayeredPane layeredPane1 = new JLayeredPane();
		layeredPane1.setForeground(Color.WHITE);
		layeredPane1.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane1.setBackground(colorBIBS2);
		layeredPane1.setBounds(10, 296, 225, 99);
		frmSpecglobtool.getContentPane().add(layeredPane1);

		textFieldOutput = new JTextField();
		textFieldOutput.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldOutput.setText("output.csv");
		textFieldOutput.setColumns(10);
		textFieldOutput.setBounds(10, 54, 161, 23);
		layeredPane1.add(textFieldOutput);

		
		JButton buttonOutputCSV = new JButton("Folder");
		buttonOutputCSV.setBounds(181, 54, 23, 23);
		layeredPane1.add(buttonOutputCSV);
		buttonOutputCSV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				JFrame frame = new JFrame("Select CSV output file name");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Select CSV output file name");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("CSV Files", "csv"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showSaveDialog(frame);
				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();
//					LOG.append("Selected file to save: " + selectedFilePath + "\n");
//					csvOutputFile = selectedFile;
					textFieldOutput.setText(selectedFilePath);

				}

			}
		}); 
		

		JLabel labelAlignFile = new JLabel("Alignment file (.csv)\r\n");
		labelAlignFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelAlignFile.setBounds(10, 32, 191, 23);
		layeredPane1.add(labelAlignFile);

		JLabel labelOutput = new JLabel("OUTPUT");
		labelOutput.setHorizontalAlignment(SwingConstants.CENTER);
		labelOutput.setFont(new Font("Tahoma", Font.BOLD, 14));
		labelOutput.setForeground(colorBIBS3);
		labelOutput.setBounds(10, 10, 205, 19);
		layeredPane1.add(labelOutput);

		progressBar = new JProgressBar();
		progressBar.setMaximum(1000);
		progressBar.setValue(0);
		progressBar.setForeground(colorBIBS1);
		progressBar.setBounds(245, 374, 434, 21);
		frmSpecglobtool.getContentPane().add(progressBar);

		JLayeredPane layeredPane2 = new JLayeredPane();
		layeredPane2.setForeground(Color.WHITE);
		layeredPane2.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane2.setBackground(colorBIBS2);
		layeredPane2.setBounds(245, 47, 434, 239);
		frmSpecglobtool.getContentPane().add(layeredPane2);

		JLabel labelParam = new JLabel("PARAMETERS");
		labelParam.setHorizontalAlignment(SwingConstants.CENTER);
		labelParam.setFont(new Font("Tahoma", Font.BOLD, 14));
		labelParam.setForeground(colorBIBS3);
		labelParam.setBounds(10, 10, 409, 19);
		layeredPane2.add(labelParam);

		JComboBox<Integer> comboBoxMultiprocess = new JComboBox<>();
		comboBoxMultiprocess.setFont(new Font(fontType, Font.PLAIN, 13));
		comboBoxMultiprocess.setEnabled(false);
		comboBoxMultiprocess.setModel(new DefaultComboBoxModel<>(new Integer[] { 2, 4, 8, 10 }));
		comboBoxMultiprocess.setBounds(371, 39, 48, 23);
		layeredPane2.add(comboBoxMultiprocess);

		JCheckBox checkBoxMultiprocess = new JCheckBox("Multiprocessing");
		checkBoxMultiprocess.setFont(new Font(fontType, Font.PLAIN, 12));
		checkBoxMultiprocess.setBounds(20, 39, 139, 23);
		checkBoxMultiprocess.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				comboBoxMultiprocess.setEnabled(!comboBoxMultiprocess.isEnabled());

			}
		});
		layeredPane2.add(checkBoxMultiprocess);

		JLabel labelThreads = new JLabel("Number of threads :");
		labelThreads.setFont(new Font(fontType, Font.PLAIN, 12));
		labelThreads.setBounds(250, 39, 111, 23);
		layeredPane2.add(labelThreads);

		JLabel labelFilter = new JLabel("Peaks filter :");
		labelFilter.setFont(new Font(fontType, Font.PLAIN, 12));
		labelFilter.setBounds(20, 106, 123, 23);
		layeredPane2.add(labelFilter);

		JComboBox<String> filterChooseBox = new JComboBox<>();
		filterChooseBox.setModel(
				new DefaultComboBoxModel<>(new String[] { "Most intense (nb)", "Rel. intensity threshold (%)" }));
		filterChooseBox.setFont(new Font(fontType, Font.PLAIN, 13));
		filterChooseBox.setBounds(102, 106, 190, 23);
		layeredPane2.add(filterChooseBox);

		JLabel labelValue = new JLabel("Value :");
		labelValue.setFont(new Font(fontType, Font.PLAIN, 12));
		labelValue.setBounds(313, 106, 48, 23);
		layeredPane2.add(labelValue);

		textFieldValueFilter = new JTextField();
		textFieldValueFilter.setText("60");
		textFieldValueFilter.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldValueFilter.setBounds(371, 106, 48, 23);
		layeredPane2.add(textFieldValueFilter);
		textFieldValueFilter.setColumns(10);

		JLabel labelPrecision = new JLabel("Mass accuracy on MS/MS fragments :");
		labelPrecision.setFont(new Font(fontType, Font.PLAIN, 12));
		labelPrecision.setBounds(20, 166, 272, 23);
		layeredPane2.add(labelPrecision);

		textFieldPRECISION = new JTextField();
		textFieldPRECISION.setText("0.02");
		textFieldPRECISION.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldPRECISION.setColumns(10);
		textFieldPRECISION.setBounds(313, 166, 48, 23);
		layeredPane2.add(textFieldPRECISION);

		JLabel labelDa = new JLabel("Da");
		labelDa.setFont(new Font(fontType, Font.PLAIN, 12));
		labelDa.setBounds(371, 166, 21, 23);
		layeredPane2.add(labelDa);

		/**
		 * The main button to launch the alignment After set all parameters, launch
		 * alignment by clicking this button
		 */
		buttonLauch = new JButton("Launch Alignments !");
		buttonLauch.setBackground(colorBIBS3);
		buttonLauch.setFont(new Font("Tahoma", Font.BOLD, 12));
		buttonLauch.setForeground(Color.WHITE);
		buttonLauch.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		buttonLauch.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (buttonLauch.isEnabled()) {
					buttonLauch.setEnabled(false);
					
					LOG.append("Selected file: " + spectraFile.getAbsolutePath() + "\n");
					LOG.append("Selected csv input file: " + textFieldInputCSV.getText() + "\n");
					LOG.append("Selected file to save: " + textFieldOutput.getText() + "\n");
					LOG.setText("");

					// get informations about the spectra file
					spectraPath = spectraFile.getAbsolutePath();

					// get informations about the csv input file
					titleScanCol = txtA.getText();
					peptideCol = txtB.getText();
					csvInputPath = textFieldInputCSV.getText();
					csvOutputPath = textFieldOutput.getText();
					
					// get informations about selected filter
					filterChoosed = filterChooseBox.getSelectedItem().toString();

					switch (filterChoosed) {
					case "Most intense (nb)":
						filterID = 1;
						break;
					case "Rel. intensity threshold (%)":
						filterID = 0;
						break;
					}

					valueFilter = Integer.parseInt(textFieldValueFilter.getText());
					precision = Double.parseDouble(textFieldPRECISION.getText());

					// get informations about multiprocessing
					isMultiprocessed = comboBoxMultiprocess.isEnabled();
					if (isMultiprocessed) {
						nbThread = (byte) Integer.parseInt(comboBoxMultiprocess.getSelectedItem().toString());
					}

					(new Thread() {
						public void run() {
							LOG.append("Column for Title = " + titleScanCol + " ----- Column for peptides = "
									+ peptideCol + "\nThe spectra file is : " + spectraPath + " - Type : " + dataType
									+ "\nThe CSV input file is : " + csvInputPath + "\nThe output file is : "
									+ csvOutputPath + "\nPeaks filter choosen : " + filterChoosed + "(" + filterID + ")"
									+ " --> " + valueFilter + "\nMultiprocesing : " + isMultiprocessed + " --> "
									+ nbThread + " thread to use\n");

						}
					}).start();

					task = new Task();
					task.execute();

				}
			}
		});

		buttonLauch.setBounds(374, 304, 176, 60);
		frmSpecglobtool.getContentPane().add(buttonLauch);

		JLabel lblPleaseThanksTo = new JLabel("*How to cite us : https://doi.org/10.1101/2022.05.31.494131 ");
		lblPleaseThanksTo.setHorizontalAlignment(SwingConstants.LEFT);
		lblPleaseThanksTo.setBounds(10, 506, 429, 27);
		lblPleaseThanksTo.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		lblPleaseThanksTo.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {

				try {
					Desktop.getDesktop().browse(new URI("https://doi.org/10.1101/2022.05.31.494131"));
				} catch (IOException e1) {

					e1.printStackTrace();
				} catch (URISyntaxException e1) {

					e1.printStackTrace();
				}
			}

		});
		frmSpecglobtool.getContentPane().add(lblPleaseThanksTo);
		lblPleaseThanksTo.setForeground(new Color(0, 0, 0));
		lblPleaseThanksTo.setFont(new Font("Tahoma", Font.ITALIC, 12));

		JLabel logoINRAE = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoINRAE.png")));
		logoINRAE.setBounds(5, 5, 100, 40);
		frmSpecglobtool.getContentPane().add(logoINRAE);

		JLabel logoUNIV = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoUNIV.png")));
		logoUNIV.setBounds(102, 5, 119, 40);
		frmSpecglobtool.getContentPane().add(logoUNIV);

		JLabel logoBIBIS = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoBIBS.png")));
		logoBIBIS.setBounds(504, 5, 180, 40);
		logoBIBIS.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		logoBIBIS.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {

				try {
					Desktop.getDesktop().browse(new URI("https://www.bibs.inrae.fr/bibs_eng/"));
				} catch (IOException e1) {

					e1.printStackTrace();
				} catch (URISyntaxException e1) {

					e1.printStackTrace();
				}
			}

		});
		frmSpecglobtool.getContentPane().add(logoBIBIS);

	}
}
