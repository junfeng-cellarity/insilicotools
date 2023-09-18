package com.insilico.application.insilicotools.gui.modeling;

import ca.odell.glazedlists.*;
import ca.odell.glazedlists.matchers.MatcherEditor;
import ca.odell.glazedlists.swing.DefaultEventListModel;
import ca.odell.glazedlists.swing.TextComponentMatcherEditor;
import com.insilico.application.insilicotools.gui.SVGTableCellRenderer;
import com.insilico.application.insilicotools.gui.filter.GhostGlassPane;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.jidesoft.swing.JideSplitPane;
import org.RDKit.MolDraw2DSVG;
import org.RDKit.MolWriter;
import org.jdesktop.swingx.JXImageView;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;
import org.RDKit.RDKFuncs;

/**
 * Created by jfeng1 on 8/8/16.
 */
public class DockingInputDialog extends JDialog {
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    String ligandString;
    String pdbString;
    JXImageView imagePanel = new JXImageView();
    String receptorName;
    SortedList<String> pdbList;
    JTextField filterField;

    public DockingInputDialog(JFrame owner, Vector<String> receptorList) {
        super(owner);
        filterField = new JTextField(10);
        JPanel filterPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        filterPanel.add(new JLabel("Filter:"));
        filterPanel.add(filterField);
        setGlassPane(new GhostGlassPane());
        imagePanel.setEditable(true);

        JButton zoomInBtn = new JButton("+");
        zoomInBtn.addActionListener(imagePanel.getZoomInAction());
        JButton zoomOutBtn = new JButton("-");
        zoomOutBtn.addActionListener(imagePanel.getZoomOutAction());
        imagePanel.setBorder(new TitledBorder("Drag to move the image"));
        JPanel contentPanel = new JPanel(new BorderLayout());
        JideSplitPane splitPane = new JideSplitPane(JideSplitPane.VERTICAL_SPLIT);
        BasicEventList<String> eventList = new BasicEventList<>();
        for(String s:receptorList){
            eventList.add(s);
        }
        pdbList = new SortedList<String>(eventList);
        MatcherEditor<String> matcherEditor = new TextComponentMatcherEditor<String>(filterField, new TextFilterator<String>() {
            @Override
            public void getFilterStrings(List<String> list, String s) {
                String[] split = s.split("_");
                if(split!=null&&split.length>0){
                    for(String s1:split){
                        list.add(s1);
                    }
                }
            }
        });
        FilterList<String> filterList = new FilterList<String>(pdbList,matcherEditor);
        DefaultEventListModel<String> listModel = new DefaultEventListModel<String>(filterList);
        splitPane.setProportionalLayout(true);
        final JList<String> l = new JList<String>(listModel);
        splitPane.add(new JScrollPane(l));
        splitPane.add(new JScrollPane(imagePanel));
        splitPane.setProportions(new double[]{0.3});
        l.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        l.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()){
                    receptorName = (String)l.getSelectedValue();
                    final GhostGlassPane glassPane = (GhostGlassPane) DockingInputDialog.this.getGlassPane();
                    glassPane.start();
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<String> receptorAndLigand = ChemFunc.getReceptorAndLigand(receptorName);
                            pdbString = receptorAndLigand.get(0);
                            ligandString = receptorAndLigand.get(1);
                            org.RDKit.MolDraw2DSVG svg = new MolDraw2DSVG(600,400);
                            org.RDKit.ROMol mol = org.RDKit.RWMol.MolFromMolBlock(ligandString);
                            svg.drawMolecule(mol);
                            svg.finishDrawing();
                            String svgString = svg.getDrawingText();
                            return SVGTableCellRenderer.convertSVGStringToImage(svgString,new Rectangle(600,400));
                        }

                        @Override
                        protected void done() {
                            try {
                                BufferedImage image = (BufferedImage) get();
                                imagePanel.setImage(image);
                                okBtn.setEnabled(true);
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(DockingInputDialog.this,e1.getMessage());
                            }finally {
                                glassPane.stop();
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        contentPanel.add(filterPanel,BorderLayout.NORTH);
        contentPanel.add(splitPane,BorderLayout.CENTER);
        JPanel btnPanel = new JPanel();
        btnPanel.add(zoomInBtn);
        btnPanel.add(zoomOutBtn);
        okBtn = new JButton("OK");
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(pdbString!=null&&ligandString!=null){
                    isCommitted=true;
                }
                setVisible(false);
            }
        });
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                setVisible(false);
            }
        });
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        contentPanel.add(btnPanel,BorderLayout.SOUTH);
        setContentPane(contentPanel);
        setSize(new Dimension(1024,768));
        setModal(true);
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public String getReceptorName() {
        return receptorName;
    }

    public String getLigandString() {
        return ligandString;
    }

    public String getPdbString() {
        return pdbString;
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        String buffer = "MLKL_BIO-0918847new\n" +
                "CDK12_4NST\n" +
                "IRAK4_BIO-0736524\n" +
                "MLKL_BIO-0918847\n" +
                "MLKL_BIO-0722425\n" +
                "IRAK4_1kocr\n" +
                "FAK_BIO-0770701\n" +
                "TTBK1_4NFN\n" +
                "3VNG_A\n" +
                "5CGJ\n" +
                "4XMB\n" +
                "4N1B_C\n" +
                "4N1B_B\n" +
                "4N1B_A\n" +
                "4L7D_A\n" +
                "4IQK_A\n" +
                "4L7C_A\n" +
                "3VNH_A\n" +
                "4IFN_A\n" +
                "IRAK4_BIO-0918018\n" +
                "DAPK1_BIO-0774048\n" +
                "IRAK4_BIO-0919408";
        Vector<String> dockingList = new Vector<String>(Arrays.asList(buffer.split("\\n")));

        dockingList.add("IRAK4_BIO-0919408");
        dockingList.add("");
        DockingInputDialog dialog = new DockingInputDialog(f,dockingList);
        dialog.setVisible(true);
        f.dispose();
    }

}
