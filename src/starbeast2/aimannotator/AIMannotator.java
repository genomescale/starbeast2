/*
 * Copyright (C) 2019 Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package starbeast2.aimannotator;

import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class AIMannotator extends TreeAnnotator {

    private enum SummaryStrategy { MEAN, MEDIAN }

    private static class AIMAnnotatorOptions {
        File inFile;
        File outFile = new File("summary.tree");
        File geneFlowFile = new File("");
//        File targetFile;
        double burninPercentage = 10.0;
        double minTreeSupport = 0.0;
        boolean useRank = true;
        SummaryStrategy summaryStrategy = SummaryStrategy.MEAN;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +                    
                    "Burn-in percentage: " + burninPercentage + "\n" +
                    "Minimal tree support for output: " + minTreeSupport + "\n" +
                    "Node height and conv. site summary: " + summaryStrategy + "\n"+
            		"Whether to consider ranked trees separated: " + useRank + "\n";

       }
    }

    public AIMannotator(AIMAnnotatorOptions options) throws IOException {
    	// define which attributes to get
    	Set<String> attributeNames = new HashSet<>();
    	attributeNames.add("length");
    	attributeNames.add("height");
    	attributeNames.add("Ne");
    	attributeNames.add("species");
    	attributeNames.add("to");
    	attributeNames.add("rates");
    	
    	
        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

        TreeSet treeSet = new FastTreeSet(options.inFile.toString(), (int) options.burninPercentage);

        if (options.useRank) {
	        // get the clades for each reassortment event in every network
	        RankedCladeSystem rankedCladeSystem = new RankedCladeSystem();
	        
	
	        // read in the tree clades
		        
	        // build the clades
	        treeSet.reset();
	        List<Integer> treeIndex = new ArrayList<>();
	        int totalTrees = 0;
	        while(treeSet.hasNext()) {
	        	treeIndex.add(rankedCladeSystem.add(treeSet.next(), true, attributeNames));
	        	totalTrees++;
	        }
	        
	        rankedCladeSystem.calculateCladeCredibilities(1);
	        int[] treeOrder = rankedCladeSystem.orderRankedTrees();
	       
	        
	        // print the trees to file
	        System.out.println("\nWriting output to " + options.outFile.getCanonicalPath()
		    	+ "...");
	        
	    	boolean nowTreePrinted = true;
	
		    try (PrintStream ps = new PrintStream(options.outFile)) {
	//	    	rankedCladeSystem.rankedTrees.get(0).tree.init(ps);
		    	ps.println("#NEXUS");
		    	ps.println("Begin trees;");
		    	for (int i = 0; i < treeOrder.length; i++) {
		    		if (rankedCladeSystem.rankedTrees.get(treeOrder[i]).credibility/totalTrees < options.minTreeSupport/100.0) {
		    			break;
		    		}else {
		    			nowTreePrinted = false;
			    		String base = "STATE_" + i + "_occurances_" + (int) rankedCladeSystem.rankedTrees.get(treeOrder[i]).credibility;
			    		Tree tree;
			    		if (options.summaryStrategy == SummaryStrategy.MEAN)
			    			tree = rankedCladeSystem.compactMetaData(treeOrder[i], attributeNames, true);
			    		else
			    			tree = rankedCladeSystem.compactMetaData(treeOrder[i], attributeNames, false);
			    		
				    	ps.println("tree " + base + 	" = " + tree.getRoot().toNewick() + ";");
				    	
				    	// try to print the meta data to a log file
				        // make filename for gene flow file
				    	System.out.println( options.outFile.getName().replace(".tree", ""));
				        File geneFlowFile = new File(options.outFile + "." + base + ".log");
				        System.out.println("\nWriting gene flow to " + geneFlowFile.getCanonicalPath() + "...");
	
				        try (PrintStream gps = new PrintStream(geneFlowFile)){
			    			rankedCladeSystem.printMetaData(gps, treeOrder[i], attributeNames);
			    		}
		    		}
		    	}
		    	
	    		ps.println("End;");
		    }	
		    System.out.println("\nDone!");
		    if (nowTreePrinted) {
		    	System.err.println("There was no tree that occured more often than " + options.minTreeSupport + " %, no tree was printed");
		    }
		    
	        treeSet.reset();
		    // print the ranked tree topologies to file
	        File rankedTreeFile = new File(options.outFile + ".ranked.log");
	
		    try (PrintStream ps = new PrintStream(rankedTreeFile)) {
		    	ps.print("Sample\tRankedTree\t\n");
		        for (int i = 0; i<treeIndex.size();i++) {
		        	
			    	ps.print(i + "\t" + arrayIndexOf(treeOrder, treeIndex.get(i)) + "\t" +"\n");
		        }
		    }
        }else {
	        // get the clades for each reassortment event in every network
	        UnrankedCladeSystem unrankedCladeSystem = new UnrankedCladeSystem();
	        
	
	        // read in the tree clades
		        
	        // build the clades
	        treeSet.reset();
	        List<Integer> treeIndex = new ArrayList<>();
	        int totalTrees = 0;
	        while(treeSet.hasNext()) {
	        	treeIndex.add(unrankedCladeSystem.add(treeSet.next(), true, attributeNames));
	        	totalTrees++;
	        }
	        
	        unrankedCladeSystem.calculateCladeCredibilities(1);
	        int[] treeOrder = unrankedCladeSystem.orderRankedTrees();
	       
	        
	        // print the trees to file
	        System.out.println("\nWriting output to " + options.outFile.getCanonicalPath()
		    	+ "...");
	        
	    	boolean nowTreePrinted = true;
	
		    try (PrintStream ps = new PrintStream(options.outFile)) {
		    	ps.println("#NEXUS");
		    	ps.println("Begin trees;");
		    	for (int i = 0; i < treeOrder.length; i++) {
		    		if (unrankedCladeSystem.unrankedTrees.get(treeOrder[i]).credibility/totalTrees < options.minTreeSupport/100.0) {
		    			break;
		    		}else {
		    			nowTreePrinted = false;
			    		String base = "STATE_" + i + "_occurances_" + (int) unrankedCladeSystem.unrankedTrees.get(treeOrder[i]).credibility;
			    		Tree tree;
			    		System.out.println("cred = " + unrankedCladeSystem.unrankedTrees.get(treeOrder[i]).credibility);
			    		if (options.summaryStrategy == SummaryStrategy.MEAN)
			    			tree = unrankedCladeSystem.compactMetaData(treeOrder[i], attributeNames, true);
			    		else
			    			tree = unrankedCladeSystem.compactMetaData(treeOrder[i], attributeNames, false);
			    		
				    	ps.println("tree " + base +	" = " + tree.getRoot().toNewick() + ";");
				    	
				    	// try to print the meta data to a log file
				        // make filename for gene flow file
				    	System.out.println( options.outFile.getName().replace(".tree", ""));
				        File geneFlowFile = new File(options.outFile + "." + base + ".log");
				        System.out.println("\nWriting gene flow to " + geneFlowFile.getCanonicalPath() + "...");
	
				        try (PrintStream gps = new PrintStream(geneFlowFile)){
				        	unrankedCladeSystem.printMetaData(gps, treeOrder[i], attributeNames);
			    		}
		    		}
		    	}		    	
	    		ps.println("End;");
		    }	
		    
		    System.out.println("\nDone!");
		    if (nowTreePrinted) {
		    	System.err.println("There was no tree that occured more often than " + options.minTreeSupport + " %, no tree was printed");
		    }
		    
	        treeSet.reset();
		    // print the ranked tree topologies to file
	        File rankedTreeFile = new File(options.outFile + ".unranked.log");
	
		    try (PrintStream ps = new PrintStream(rankedTreeFile)) {
		    	ps.print("Sample\tTopology\t\n");
		        for (int i = 0; i<treeIndex.size();i++) {
		        	
			    	ps.print(i + "\t" + arrayIndexOf(treeOrder, treeIndex.get(i)) + "\t" +"\n");
		        }
		    }

        }
	    
    }      
    
    private int arrayIndexOf(int[] array, int val) {
    	for (int i=0;i< array.length;i++)
    		if (array[i]==val)
    			return i;
    	return -1;
    }

 
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(AIMAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Isolation with Migration Annotator");

        JLabel logFileLabel = new JLabel("Isolation with migration species tree file:");
        JLabel outFileLabel = new JLabel("Output file:");
//        JLabel geneFlowLabel = new JLabel("Gene Flow file");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel minLabel = new JLabel("Minal tree support for output:");
        JLabel summaryMethodLabel = new JLabel("Node height summary method:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");
        
//        JTextField geneFlowFilename = new JTextField(20);
//        geneFlowFilename.setText(options.geneFlowFile.getName());
//        geneFlowFilename.setEditable(false);
//        JButton geneFileButton = new JButton("Choose File");

        
//        JTextField targetFilename = new JTextField(20);
//        targetFilename.setEditable(false);
//        JButton targetFileButton = new JButton("Choose File");

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);
        
        JSlider minTreeSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.minTreeSupport));
        minTreeSlider.setMajorTickSpacing(50);
        minTreeSlider.setMinorTickSpacing(10);
        minTreeSlider.setPaintTicks(true);
        minTreeSlider.setPaintLabels(true);
        minTreeSlider.setSnapToTicks(true);

        
//        JTextField geneFlow = new JTextField(20);
//        geneFlow.setText("none");
//        geneFlow.setEditable(true);


        JComboBox<SummaryStrategy> heightMethodCombo = new JComboBox<>(SummaryStrategy.values());

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
//                        .addComponent(targetFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(minLabel)
                        .addComponent(summaryMethodLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
//                        .addComponent(geneFlowFilename)
//                        .addComponent(targetFilename)
                        .addComponent(burninSlider)
                        .addComponent(minTreeSlider)
                        .addComponent(heightMethodCombo))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton)
//                        .addComponent(geneFileButton)
//                        .addComponent(targetFileButton)
                        ));

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(geneFlowLabel)
//                        .addComponent(geneFlowFilename,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE)
//                        .addComponent(geneFileButton))
//                .addGroup(layout.createParallelGroup()
//                        .addComponent(targetFileLabel)
//                        .addComponent(targetFilename,
//                                GroupLayout.PREFERRED_SIZE,
//                                GroupLayout.DEFAULT_SIZE,
//                                GroupLayout.PREFERRED_SIZE)
//                        .addComponent(targetFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(minLabel)
                        .addComponent(minTreeSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
               .addGroup(layout.createParallelGroup()
                        .addComponent(summaryMethodLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
        		);

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.summaryStrategy = (SummaryStrategy)heightMethodCombo.getSelectedItem();
            options.minTreeSupport = minTreeSlider.getValue();
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select AIM species tree file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });
        
//        JFileChooser geneFileChooser = new JFileChooser();
//        geneFileButton.addActionListener(e -> {
//        	geneFileChooser.setDialogTitle("Select file name for gene flow log.");
//            if (options.inFile != null)
//            	geneFileChooser.setCurrentDirectory(options.inFile);
//            else
//            	geneFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//
//            geneFileChooser.setSelectedFile(options.outFile);
//            int returnVal = geneFileChooser.showOpenDialog(dialog);
//
//            if (returnVal == JFileChooser.APPROVE_OPTION) {
//                options.outFile = geneFileChooser.getSelectedFile();
//                outFilename.setText(geneFileChooser.getSelectedFile().getName());
//            }
//        });

        
//        JFileChooser targetFileChooser = new JFileChooser();
//        targetFileButton.addActionListener(e -> {
//            targetFileChooser.setDialogTitle("Select target tree name.");
//            if (options.inFile != null)
//            	targetFileChooser.setCurrentDirectory(options.inFile);
//            else
//            	targetFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
//
//            targetFileChooser.setSelectedFile(options.outFile);
//            int returnVal = targetFileChooser.showOpenDialog(dialog);
//
//            if (returnVal == JFileChooser.APPROVE_OPTION) {
//                options.targetFile = targetFileChooser.getSelectedFile();
//                outFilename.setText(targetFileChooser.getSelectedFile().getName());
//            }
//        });
             			
        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Isolation with Migration Annotator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "Isolation with Migration Annotator - produces summaries of AIM species tree files.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-positions {mean,median} Choose position summary method.\n"
                    + "                         (default mean)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-minTreeSupport percentage    Choose the minimum support a.\n"
                    + "                         Tree has to have to be included (Default 0%)\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'summary.tree'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve ACGAnnotator options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, AIMAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;

                case "-positions":
                    if (args.length<=i+1) {
                        printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");
                    }

                    if (args[i+1].toLowerCase().equals("mean")) {
                        options.summaryStrategy = SummaryStrategy.MEAN;

                        i += 1;
                        break;
                    }

                    if (args[i+1].toLowerCase().equals("median")) {
                        options.summaryStrategy = SummaryStrategy.MEDIAN;

                        i += 1;
                        break;
                    }

                    printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");

                     
                case "-minTreeSupport":
                    if (args.length<=i+1)
                        printUsageAndError("-minTreeSupport must be followed by a number (percent)");

                    try {
                        options.minTreeSupport = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.minTreeSupport<0 || options.burninPercentage>100) {
                        printUsageAndError("minTreeSupport percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;
                    
                case "-userank":
                    if (args.length<=i+1)
                        printUsageAndError("-userank by true or false");

                    try {
                        options.useRank = Boolean.parseBoolean(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing useRank.");
                    }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	AIMAnnotatorOptions options = new AIMAnnotatorOptions();

    	if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new AIMannotator(options);
        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }

}