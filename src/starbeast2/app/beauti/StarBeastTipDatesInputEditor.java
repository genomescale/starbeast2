package starbeast2.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.TraitSet;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiPanel;
import beastfx.app.inputeditor.GuessPatternDialog;
import beastfx.app.util.FXUtils;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.Node;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import starbeast2.SpeciesTree;

import java.text.DateFormat;
import java.util.*;

public class StarBeastTipDatesInputEditor extends BEASTObjectInputEditor {

    public StarBeastTipDatesInputEditor(BeautiDoc doc) {
        super(doc);
    }
    private static final long serialVersionUID = 1L;

    DateFormat dateFormat = DateFormat.getDateInstance();

    @Override
    public Class<?> type() {
        return SpeciesTree.class;
    }
    SpeciesTree tree;
    TraitSet traitSet;
    ComboBox<TraitSet.Units> unitsComboBox;
    ComboBox<String> relativeToComboBox;
    List<String> taxa;
//    Object[][] tableData;
    ObservableList<TipDate> tableData;
    TableView<TipDate> table;
    String m_sPattern = ".*(\\d\\d\\d\\d).*";
    ScrollPane scrollPane;
    List<Taxon> taxonsets;

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
        m_bAddButtons = addButtons;
        this.itemNr = itemNr;

        pane.getChildren().clear();
        pane = FXUtils.newHBox();

        if (itemNr >= 0) {
            tree = (SpeciesTree) ((List<?>) input.get()).get(itemNr);
        } else {
            tree = (SpeciesTree) input.get();
        }
        if (tree != null) {
            try {
                m_input = ((BEASTInterface) tree).getInput("trait");
            } catch (Exception e1) {
                // TODO Auto-generated catch block
                e1.printStackTrace();
            }
            m_beastObject = tree;
            traitSet = tree.getDateTrait();

            VBox box = FXUtils.newVBox();

            CheckBox useTipDates = new CheckBox("Use tip dates");
            useTipDates.setSelected(traitSet != null);
            useTipDates.selectedProperty().addListener(new ChangeListener<>() {
                @Override
                public void changed(ObservableValue<? extends Boolean> observable, Boolean oldValue, Boolean newValue) {
                    try {
                        if (newValue) {
                            if (traitSet == null) {
                                traitSet = new TraitSet();
                                traitSet.initByName("traitname", "date",
                                        "taxa", tree.getTaxonset(),
                                        "value", "");
                                traitSet.setID("dateTrait.t:" + BeautiDoc.parsePartition(tree.getID()));
                            }
                            tree.setDateTrait(traitSet);
                        } else {
                            tree.setDateTrait(null);
                        }

                        refreshPanel();
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }

                }
            });
            HBox box2 = FXUtils.newHBox();
            box2.getChildren().add(useTipDates);
//            box2.add(Box.createHorizontalGlue());
            box.getChildren().add(box2);

            if (traitSet != null) {
                box.getChildren().add(createButtonBox());
//                box.getChildren().add(createListBox());
            }
            getChildren().add(box);
        }
    } // init

    public class TipDate {
        String taxon;
        String date;
        Double age;

        TipDate(String taxon, String date, Double age) {
            this.taxon = taxon;
            this.date = date;
            this.age = age;
        }

        public String getTaxon() {
            return taxon;
        }
        public void setTaxon(String taxon) {
            this.taxon = taxon;
        }
        public String getDate() {
            return date;
        }
        public void setDate(String date) {
            this.date = date;
        }
        public Double getAge() {
            return age;
        }
        public void setAge(Double age) {
            this.age = age;
        }
    } // class TipDate
    
    private Node createListBox() {
        taxa = traitSet.taxaInput.get().asStringList();
        List<TipDate> list = new ArrayList<>();
        for (String taxon : taxa) {
            list.add(new TipDate(taxon, "0", 0.0));
        }
        tableData = FXCollections.observableArrayList(list);

        // String[] columnData = new String[]{"Name", "Date (raw value)", "Height"};
        // tableData = new Object[taxa.size()][3];

        table = new TableView<>();
        table.setPrefWidth(800);
        table.setMinSize(doc.beauti.frame.getWidth()-12, doc.beauti.frame.getHeight()-230);
        BeautiPanel.resizeList.clear();
        BeautiPanel.resizeList.add(table);

        table.setEditable(true);

        TableColumn<TipDate,String> col1 = new TableColumn<>("Name");
        col1.setPrefWidth(500);
        col1.setEditable(false);
        col1.setCellValueFactory(
                new PropertyValueFactory<>("Taxon")
        );
        table.getColumns().add(col1);

        TableColumn<TipDate,String> col2 = new TableColumn<>("Date");
        col2.setPrefWidth(150);
        col2.setEditable(true);
        col2.setCellValueFactory(
                new PropertyValueFactory<>("Date")
        );
        table.getColumns().add(col2);

        col2.setCellFactory(TextFieldTableCell.forTableColumn());
        col2.setOnEditCommit(
                event -> {
                    String newValue = event.getNewValue();
                    TipDate tipDate = event.getRowValue();
                    tipDate.setDate(newValue);
                    convertTableDataToTrait();
                    convertTraitToTableData();
                }
        );

        TableColumn<TipDate,Double> col3 = new TableColumn<>("Height");
        col3.setPrefWidth(150);
        col3.setEditable(false);
        col3.setCellValueFactory(
                new PropertyValueFactory<TipDate,Double>("Age")
        );
        table.getColumns().add(col3);

        table.setItems(tableData);

        convertTraitToTableData();

        return table;

        // set up table.
        // special features: background shading of rows
        // custom editor allowing only Date column to be edited.
//        table = new Table(tableData, columnData) {
//            private static final long serialVersionUID = 1L;
//
//            // method that induces table row shading
//            @Override
//            public Component prepareRenderer(TableCellRenderer renderer, int Index_row, int Index_col) {
//                Component comp = super.prepareRenderer(renderer, Index_row, Index_col);
//                //even index, selected or not selected
//                if (isCellSelected(Index_row, Index_col)) {
//                    comp.setBackground(Color.lightGray);
//                } else if (Index_row % 2 == 0 && !isCellSelected(Index_row, Index_col)) {
//                    comp.setBackground(new Color(237, 243, 255));
//                } else {
//                    comp.setBackground(Color.white);
//                }
//                return comp;
//            }
//        };
//
//        // set up editor that makes sure only doubles are accepted as entry
//        // and only the Date column is editable.
//        table.setDefaultEditor(Object.class, new TableCellEditor() {
//            JTextField m_textField = new JTextField();
//            int m_iRow,
//                    m_iCol;
//
//            @Override
//            public boolean stopCellEditing() {
//                table.removeEditor();
//                String text = m_textField.getText();
////                try {
////                    Double.parseDouble(text);
////                } catch (Exception e) {
////                	try {
////                		Date.parse(text);
////                	} catch (Exception e2) {
////                        return false;
////					}
////                }
//                tableData[m_iRow][m_iCol] = text;
//                convertTableDataToTrait();
//                convertTraitToTableData();
//                return true;
//            }
//
//            @Override
//            public boolean isCellEditable(EventObject anEvent) {
//                return table.getSelectedColumn() == 1;
//            }
//
//            @Override
//            public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int rowNr, int colNr) {
//                if (!isSelected) {
//                    return null;
//                }
//                m_iRow = rowNr;
//                m_iCol = colNr;
//                m_textField.setText((String) value);
//                return m_textField;
//            }
//
//            @Override
//            public boolean shouldSelectCell(EventObject anEvent) {
//                return false;
//            }
//
//            @Override
//            public void removeCellEditorListener(CellEditorListener l) {
//            }
//
//            @Override
//            public Object getCellEditorValue() {
//                return null;
//            }
//
//            @Override
//            public void cancelCellEditing() {
//            }
//
//            @Override
//            public void addCellEditorListener(CellEditorListener l) {
//            }
//        });
//        int fontsize = table.getFont().getSize();
//        table.setRowHeight(24 * fontsize / 13);
//        scrollPane = new ScrollPane(table);

// AJD: This ComponentListener breaks the resizing of the tip dates table, so I have removed it.
//        scrollPane.addComponentListener(new ComponentListener() {
//            @Override
//            public void componentShown(ComponentEvent e) {}
//
//            @Override
//            public void componentResized(ComponentEvent e) {
//                Component c = (Component) e.getSource();
//                while (c.getParent() != null && !(c.getParent() instanceof JSplitPane)) {
//                    c = c.getParent();
//                }
//                if (c.getParent() != null) {
//                    Dimension preferredSize = c.getSize();
//                    preferredSize.height = Math.max(preferredSize.height - 170, 0);
//                    preferredSize.width = Math.max(preferredSize.width - 25, 0);
//                    scrollPane.setPreferredSize(preferredSize);
//                } else if (doc.getFrame() != null) {
//                    Dimension preferredSize = doc.getFrame().getSize();
//                    preferredSize.height = Math.max(preferredSize.height - 170, 0);
//                    preferredSize.width = Math.max(preferredSize.width - 25, 0);
//                    scrollPane.setPreferredSize(preferredSize);
//                }
//            }
//
//            @Override
//            public void componentMoved(ComponentEvent e) {}
//
//            @Override
//            public void componentHidden(ComponentEvent e) {}
//        });

    } // createListBox

    /* synchronise table with data from traitSet BEASTObject */
    private void convertTraitToTableData() {
        for (int i = 0; i < tableData.size(); i++) {
            tableData.get(i).setTaxon(taxa.get(i));
            tableData.get(i).setDate("0");
            tableData.get(i).setAge(0.0);
        }
        String[] traits = traitSet.traitsInput.get().split(",");
        for (String trait : traits) {
            trait = trait.replaceAll("\\s+", " ");
            String[] strs = trait.split("=");
            if (strs.length != 2) {
                break;
                //throw new Exception("could not parse trait: " + trait);
            }
            String taxonID = normalize(strs[0]);
            int taxonIndex = taxa.indexOf(taxonID);
//            if (taxonIndex < 0) {
//                throw new Exception("Trait (" + taxonID + ") is not a known taxon. Spelling error perhaps?");
//            }
            if (taxonIndex >= 0) {
                tableData.get(taxonIndex).setDate(normalize(strs[1]));
                tableData.get(taxonIndex).setTaxon(taxonID);
            } else {
                Log.warning.println("WARNING: File contains taxon " + taxonID + " that cannot be found in alignment");
            }
        }
        if (traitSet.traitNameInput.get().equals(TraitSet.DATE_BACKWARD_TRAIT)) {
            Double minDate = Double.MAX_VALUE;
            for (int i = 0; i < tableData.size(); i++) {
                minDate = Math.min(minDate, parseDate(tableData.get(i).getDate()));
            }
            for (int i = 0; i < tableData.size(); i++) {
                tableData.get(i).setAge( parseDate(tableData.get(i).getDate()) - minDate );
            }
        } else {
            Double maxDate = 0.0;
            for (int i = 0; i < tableData.size(); i++) {
                maxDate = Math.max(maxDate, parseDate(tableData.get(i).getDate()));
            }
            for (int i = 0; i < tableData.size(); i++) {
                tableData.get(i).setAge( maxDate - parseDate(tableData.get(i).getDate()) );
            }
        }
        table.refresh();

//        if (table != null) {
//            for (int i = 0; i < tableData.length; i++) {
//                table.setValueAt(tableData.get(i).[1], i, 1);
//                table.setValueAt(tableData.get(i).[2], i, 2);
//            }
//        }
    } // convertTraitToTableData

    private double parseDate(String str) {
        // default, try to interpret the string as a number
        try {
            return Double.parseDouble(str);
        } catch (NumberFormatException e) {
            // does not look like a number, try parsing it as a date
            if (str.matches(".*[a-zA-Z].*")) {
                str = str.replace('/', '-');
            }

            //try {

            // unfortunately this deprecated date parser is the most flexible around at the moment...
            long time = Date.parse(str);
            Date date = new Date(time);

            // AJD
            // Ideally we would use a non-deprecated method like this one instead but it seems to have
            // far less support for different date formats.
            // for example it fails on "12-Oct-2014"
            //dateFormat.setLenient(true);
            //Date date = dateFormat.parse(str);

            Calendar calendar = dateFormat.getCalendar();
            calendar.setTime(date);

            // full year (e.g 2015)
            int year = calendar.get(Calendar.YEAR);
            double days = calendar.get(Calendar.DAY_OF_YEAR);

            double daysInYear = 365.0;

            if (calendar instanceof GregorianCalendar &&(((GregorianCalendar) calendar).isLeapYear(year))) {
                daysInYear = 366.0;
            }

            double dateAsDecimal = year + days/daysInYear;

            return dateAsDecimal;
            //}
            //catch (ParseException e1) {
            //    System.err.println("*** WARNING: Failed to parse '" + str + "' as date using dateFormat " + dateFormat);
            //}
        }
        //return 0;
    } // parseStrings

    private String normalize(String str) {
        if (str.charAt(0) == ' ') {
            str = str.substring(1);
        }
        if (str.endsWith(" ")) {
            str = str.substring(0, str.length() - 1);
        }
        return str;
    }

    /**
     * synchronise traitSet BEAST object with table data
     */
    private void convertTableDataToTrait() {
        String trait = "";
        for (int i = 0; i < tableData.size(); i++) {
            trait += taxa.get(i) + "=" + tableData.get(i).getDate();
            if (i < tableData.size() - 1) {
                trait += ",\n";
            }
        }
        try {
            traitSet.traitsInput.setValue(trait, traitSet);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * create box with comboboxes for selectin units and trait name *
     */
    private HBox createButtonBox() {
        HBox buttonBox = FXUtils.newHBox();

        Label label = new Label("Dates specified as: ");
//        label.setMaximumSize(label.getPreferredSize());
        buttonBox.getChildren().add(label);
        unitsComboBox = new ComboBox<>(FXCollections.observableArrayList( TraitSet.Units.values() ));
        unitsComboBox.getSelectionModel().select(traitSet.unitsInput.get());
        unitsComboBox.setOnAction(e -> {
            String selected = unitsComboBox.getSelectionModel().getSelectedItem().toString();
            try {
                traitSet.unitsInput.setValue(selected, traitSet);
                //System.err.println("Traitset is now: " + m_traitSet.m_sUnits.get());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        });
//        Dimension d = unitsComboBox.getPreferredSize();
//        unitsComboBox.setMaximumSize(new Dimension(Integer.MAX_VALUE, unitsComboBox.getPreferredSize().height));
//        unitsComboBox.setSize(d);
        buttonBox.getChildren().add(unitsComboBox);

        relativeToComboBox = new ComboBox<>(
                FXCollections.observableArrayList("Since some time in the past", "Before the present") );
        relativeToComboBox.setTooltip(new Tooltip("Whether dates go forward or backward"));
        if (traitSet.traitNameInput.get().equals(TraitSet.DATE_BACKWARD_TRAIT)) {
            relativeToComboBox.getSelectionModel().select(1);
        } else {
            relativeToComboBox.getSelectionModel().select(0);
        }
        relativeToComboBox.setOnAction(e -> {
            String selected = TraitSet.DATE_BACKWARD_TRAIT;
            if (relativeToComboBox.getSelectionModel().getSelectedIndex() == 0) {
                selected = TraitSet.DATE_FORWARD_TRAIT;
            }
            try {
                traitSet.traitNameInput.setValue(selected, traitSet);
                Log.warning.println("Relative position is now: " + traitSet.traitNameInput.get());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
            convertTraitToTableData();
        });
//        relativeToComboBox.setMaximumSize(new Dimension(Integer.MAX_VALUE, relativeToComboBox.getPreferredSize().height));
        buttonBox.getChildren().add(relativeToComboBox);

//        buttonBox.getChildren().add(Box.createHorizontalGlue());
        Button guessButton = new Button("Auto-configure");
        guessButton.setTooltip(new Tooltip("Automatically configure dates based on taxon names"));
        guessButton.setId("Guess");
        guessButton.setOnAction(e -> {
            GuessPatternDialog dlg = new GuessPatternDialog(null, m_sPattern);
            dlg.allowAddingValues();
            String trait = "";
            switch (dlg.showDialog("Guess dates")) {
                case canceled:
                    return;
                case trait:
                    trait = dlg.getTrait();
                    break;
                case pattern:
                    for (String taxon : taxa) {
                        String match = dlg.match(taxon);
                        if (match == null) {
                            return;
                        }
                        double date = parseDate(match);
                        if (trait.length() > 0) {
                            trait += ",";
                        }
                        trait += taxon + "=" + date;
                    }
                    break;
            }
            try {
                traitSet.traitsInput.setValue(trait, traitSet);
                convertTraitToTableData();
                convertTableDataToTrait();
            } catch (Exception ex) {
                // TODO: handle exception
            }
            refreshPanel();
        });
        buttonBox.getChildren().add(guessButton);


        Button clearButton = new Button("Clear");
        clearButton.setTooltip(new Tooltip("Set all dates to zero"));
        clearButton.setOnAction(e -> {
            try {
                traitSet.traitsInput.setValue("", traitSet);
            } catch (Exception ex) {
                // TODO: handle exception
            }
            refreshPanel();
        });
        buttonBox.getChildren().add(clearButton);

        return buttonBox;
    } // createButtonBox
}
