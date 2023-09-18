package com.insilico.application.insilicotools.chart;

import com.insilico.application.insilicotools.data.MolProperty;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.util.ImageUtil;
import com.google.common.base.Charsets;
import com.google.common.io.Resources;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.general.DatasetChangeEvent;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.Layer;
import org.jfree.ui.RefineryUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.util.Vector;

/**
 * Created by jfeng1 on 6/16/16.
 */
public class ClearancePlot extends JPanel{
    private final XYLineAnnotation clearanceAnnotation;

    private double minX = Double.MAX_VALUE;
    private double minY = Double.MAX_VALUE;
    private double maxX = Double.MIN_VALUE;
    private double maxY = Double.MIN_VALUE;

    private Range defaultDomainRange;
    private Range defaultRangeRange;

    private final XYPlot plot;
    private final JFreeChart chart;
    private final Scatterplot scatterplot;
    private final SelectionXYDataset data;
    private int mol_series = -1;

    public ClearancePlot() {
        this("<html><body><font color=\"#0000FF\">Clearance Route Prediction</font></body></html>","PC1","PC2",true,true);
    }

    public ClearancePlot(String title, String xAxis, String yAxis, boolean createLegend, boolean showToolTip) {
        super(new BorderLayout());


        JLabel titleLabel = new JLabel(title, JLabel.CENTER);
        titleLabel.setFont(titleLabel.getFont().deriveFont(18f));
        add(titleLabel, BorderLayout.NORTH);
        scatterplot = new Scatterplot(null, xAxis, yAxis, createLegend) {
            public void chartMouseClicked(ChartMouseEvent event) {
                super.chartMouseClicked(event);
//
//                Rectangle2D dataArea = getChartRenderingInfo().getPlotInfo().getDataArea();
//                Point pt = event.getTrigger().getPoint();
//
//                if (isPointOnEdge(PMIPct, pt, dataArea)) {
//                    selectPointsInShape(pmiShape);
//                }
            }
        };

        chart = scatterplot.getChart();
        plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        data = (SelectionXYDataset) plot.getDataset();
        loadRefPoints();
        if (showToolTip) {
//            for(int i=0;i<scatterplot.getPlot().getSeriesCount();i++){
//                scatterplot.getPlot().getRenderer().setSeriesToolTipGenerator(i, new XYToolTipGenerator() {
//                    @Override
//                    public String generateToolTip(XYDataset dataset, int series, int item) {
//                        try {
//                            System.out.println(series);
//                            Object key = ((SelectionXYDataset) dataset).getItem(series, item).getKey();
//                            String smiles = (key == null) ? null : key.toString();
//
////                        return "<html><body><img src=\"com.My.application.insilicotools.data:image/gif;base64,R0lGODlhUAAPAKIAAAsLav///88PD9WqsYmApmZmZtZfYmdakyH5BAQUAP8ALAA AAABQAA8AAAPbWLrc/jDKSVe4OOvNu/9gqARDSRBHegyGMahqO4R0bQcjIQ8E4BMCQc930JluyGRmdAAcdiigMLVrpTYWy5FKM1IQe+Mp+L4rphz+qIOBAUYeCY4p2tGrJZeH9y79mZsawFoaIRxF3JyiYxuHiMGb5KTkpFvZj4ZbYeCiXaiKBwnxh4fnt9e3ktgZyHhrChinONs3cFAShFF2JhvCZlG5uchYNun5eedRxMAF15XEFRXgZWWdciuM8GCmdSQ84lLQY5R14wDB5Lyon4ubwS7jx9NcV9/j5+g4JADs=\" alt=\"\" width=\"80\" height=\"15\" /></body></html>";
//                            return String.format("<html><body>%s</body></html>",smiles);
//                        }
//                        catch (Exception err) {
//                            err.printStackTrace();
//                            return null;
//                        }
//                    }
//                });
//            }
            scatterplot.setToolTipGenerator(new XYToolTipGenerator() {
                public String generateToolTip(XYDataset dataset, int series, int item) {
                    try {
                        Object key = ((SelectionXYDataset) dataset).getItem(series, item).getKey();
                        if(key!=null&&key instanceof PropertyMolecule) {
                            PropertyMolecule molecule = (PropertyMolecule)key;
                            MolProperty property = molecule.getProperty("Clearance Route");
                            if(property!=null) {
                                return ImageUtil.smi2ImgUrl(String.format("%s %s %s", molecule.getSmiles(), molecule.getName(), property.getProperty()));
                            }else{
                                return ImageUtil.smi2ImgUrl(String.format("%s %s", molecule.getSmiles(), molecule.getName()));
                            }
                        }else{
                            return null;
                        }
                    }
                    catch (Exception err) {
                        err.printStackTrace();
                        return null;
                    }
                }
            });
        }

        clearanceAnnotation = new XYLineAnnotation(minX,4*minX+1,maxX,4*maxX+1);
        plot.getRenderer().addAnnotation(clearanceAnnotation);

        data.addChangeListener(new DatasetChangeListener() {
            public void datasetChanged(DatasetChangeEvent event) {
                Range domainRange = plot.getDomainAxis().getRange();
                Range rangeRange = plot.getRangeAxis().getRange();
                defaultDomainRange = new Range(Math.min(minX, domainRange.getLowerBound()), Math.max(maxX, domainRange.getUpperBound()));
                defaultRangeRange = new Range(Math.min(minY, rangeRange.getLowerBound()), Math.max(maxY, rangeRange.getUpperBound()));
            }
        });

        plot.getDomainAxis().setAutoRange(false);
        plot.getRangeAxis().setAutoRange(false);

        reset();

        add(scatterplot, BorderLayout.CENTER);

    }

    public Scatterplot getScatterplot() {
        return scatterplot;
    }

    public void selectAll() {
        scatterplot.selectAll();
    }

    public void selectAll(int series) {
        scatterplot.selectAll(series);
    }

    public synchronized void addMouseListener(MouseListener l) {
        scatterplot.addMouseListener(l);
    }

    public synchronized void removeMouseListener(MouseListener l) {
        scatterplot.removeMouseListener(l);
    }

    public void setBackground(Color bg) {
        if (scatterplot != null) {
            scatterplot.setBackground(bg);
        } else {
            super.setBackground(bg);
        }
    }

    public Color getBackground() {
        return scatterplot != null ? scatterplot.getBackground() : super.getBackground();
    }

    public void addSelectionListener(SelectionListener l) {
        scatterplot.addSelectionListener(l);
    }

    public void removeSelectionListener(SelectionListener l) {
        scatterplot.removeSelectionListener(l);
    }

    public boolean selectPoint(int series, Object key) {
        return scatterplot.selectPoint(series, key);
    }

    public Object[] getSelectedPoints() {
        return scatterplot.getSelectedPoints();
    }

    public void reset() {
        plot.getDomainAxis().setRange(defaultDomainRange, false, false);
        plot.getRangeAxis().setRange(defaultRangeRange, false, false);
    }

    private Shape createPolygon(double[][] points) {
        Path2D polygon = new Path2D.Double();
        for (int i = 0; i < points.length; i++) {
            if (i==0)
                polygon.moveTo(points[i][0], points[i][1]);
            else
                polygon.lineTo(points[i][0], points[i][1]);
        }
        return polygon;
    }

    public int addSeries(String name) {
        int series = scatterplot.addSeries(name, Color.BLUE);
        if(name.equals("Molecule")){
            mol_series = series;
        }
        scatterplot.getPlot().getRenderer().setSeriesShape(series, new Ellipse2D.Double(-4,-4,8,8));
        return series;
    }

    public void addPoint(int series, double x, double y, PropertyMolecule mol) {
        scatterplot.addPoint(series, x, y, mol);
    }


    public void writeChartAsSVG(OutputStream out) throws IOException {
        scatterplot.writeChartAsSVG(out);
    }

    public void writeChartAsPDF(OutputStream out) throws IOException, ScatterplotException {
        scatterplot.writeChartAsPDF(out);
    }

    public void resetRanges(){
        defaultDomainRange = new Range(minX-1, maxX+1);
        defaultRangeRange = new Range(minY-1,maxY+1);

        reset();

    }

    public void clear() {
        if(mol_series!=-1) {
            scatterplot.removeSeries(mol_series);
        }
        resetRanges();
    }

    public void clearSelection() {
        scatterplot.clearSelection();
    }

    class AnnotationIcon implements Icon {

        private final Stroke stroke;
        private final Color color;

        public AnnotationIcon(Stroke stroke, Color color) {
            this.stroke = stroke;
            this.color = color;
        }

        public int getIconHeight() {
            return 5;
        }

        public int getIconWidth() {
            return 15;
        }

        public void paintIcon(Component component, Graphics graphics, int x, int y) {
            Graphics2D g = (Graphics2D) graphics;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
            g.setStroke(stroke);
            g.setColor(color);
            g.draw(new Line2D.Double(x, y + getIconHeight() / 2f, x + getIconWidth(), y + getIconHeight() / 2f));
        }
    }

    private void loadRefPoints(){
        int series_m = addSeries("Metabolic Clearance");
        int series_r = addSeries("Renal Clearance");
        int series_b = addSeries("Bile Clearance");
        scatterplot.getPlot().getRenderer().setSeriesPaint(series_m,new Color(255,0,0,127));
        scatterplot.getPlot().getRenderer().setSeriesPaint(series_b,new Color(0,255,0,127));
        scatterplot.getPlot().getRenderer().setSeriesPaint(series_r,new Color(0,0,255,127));
        URL resource = getClass().getClassLoader().getResource("clearance_ref.sdf");
        String text = null;
        try {
            text = Resources.toString(resource, Charsets.UTF_8);
        } catch (IOException e) {
            JOptionPane.showMessageDialog(ClearancePlot.this,e.getMessage());
            e.printStackTrace();
        }
        oemolistream ifs = new oemolistream();
        ifs.openstring(text);
        ifs.SetFormat(OEFormat.SDF);
        OEGraphMol mol = new OEGraphMol();
        Vector m_v = new Vector();
        Vector r_v = new Vector();
        Vector b_v = new Vector();
        while(oechem.OEReadMolecule(ifs,mol)){
            String route = oechem.OEGetSDData(mol,"Code (mechanism and clearance degree index)");
            if(route==null){
                continue;
            }

            double pc_1 = Double.parseDouble(oechem.OEGetSDData(mol,"PC1"));
            double pc_2 = Double.parseDouble(oechem.OEGetSDData(mol,"PC2"));
            Vector v = new Vector();
            v.add(pc_1);
            v.add(pc_2);
            v.add(new PropertyMolecule(mol));

            if(route.equals("M1")){
                m_v.add(v);
//                addPoint(series_m,pc_1,pc_2,oechem.OEMolToSmiles(mol));
            }else if(route.equals("R1")){
                r_v.add(v);
//                addPoint(series_r,pc_1,pc_2,oechem.OEMolToSmiles(mol));
            }else{
                b_v.add(v);
//                addPoint(series_b,pc_1,pc_2,oechem.OEMolToSmiles(mol));
            }
            if(pc_1<minX){
                minX = pc_1;
            }
            if(pc_2<minY){
                minY = pc_2;
            }
            if(pc_1>maxX){
                maxX = pc_1;
            }
            if(pc_2>maxY){
                maxY = pc_2;
            }

        }
        ifs.close();
        defaultDomainRange = new Range(minX-1, maxX+1);
        defaultRangeRange = new Range(minY-1,maxY+1);

        double y1 = minY-1;
        double x11 = (minY-2)/4.0-1;
        double x12 = (minY-2)/4.0+1;

        double y2 = maxY+1;
        double x21 = (maxY)/4.0-1;
        double x22 = (maxY)/4.0+1;

        double [] p11 = {0,1};
        double [] p21 = {-5.0/34.0,7.0/17.0};
        double [] p12 = {(2-minY)/1.33,minY-1};
        double [] p22 = {(minY-1.5)/0.6,minY-1};

        Shape shape_gray1 = createPolygon(new double[][]{{x11,y1},{x12,y1},{x22,y2},{x21,y2}});
        XYShapeAnnotation gray1 = new XYShapeAnnotation(shape_gray1,null,null, new Color(192,192,192,127));
        plot.getRenderer().addAnnotation(gray1,Layer.BACKGROUND);

        Shape shape_green = createPolygon(new double[][]{p11,p12,p22,p21,p11});
        XYShapeAnnotation green2 = new XYShapeAnnotation(shape_green,null,null, new Color(103,246,101,127));
        plot.getRenderer().addAnnotation(green2,Layer.BACKGROUND);

        for(Object v :m_v){
            Vector v1 = (Vector)v;
            double x = (Double)v1.get(0);
            double y = (Double)v1.get(1);
            PropertyMolecule molecule = (PropertyMolecule)v1.get(2);
            addPoint(series_m,x,y,molecule);
        }
        for(Object v :r_v){
            Vector v1 = (Vector)v;
            double x = (Double)v1.get(0);
            double y = (Double)v1.get(1);
            PropertyMolecule molecule = (PropertyMolecule)v1.get(2);
            addPoint(series_r,x,y,molecule);
        }
        for(Object v :b_v){
            Vector v1 = (Vector)v;
            double x = (Double)v1.get(0);
            double y = (Double)v1.get(1);
            PropertyMolecule molecule = (PropertyMolecule) v1.get(2);
            addPoint(series_b,x,y,molecule);
        }

    }

    public static void main(String[] args) {
        final ClearancePlot plot = new ClearancePlot();

//        int series = plot.addSeries("BIO Number");
//        for (int i = 0; i < 10; i++) {
//            double abs1 = Math.abs(Math.random());
//            double abs2 = 0.5+Math.abs(Math.random())/2;
//            plot.addPoint(series, abs1, abs2, "CCCCC");
//        }
//
        final JFrame f = new JFrame("Clearance Example");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(plot, BorderLayout.CENTER);
        JButton exportPdfButton = new JButton("Export PDF");
        f.getContentPane().add(exportPdfButton, BorderLayout.SOUTH);
        exportPdfButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JFileChooser fileChooser = new JFileChooser();
                if (fileChooser.showSaveDialog(f) == JFileChooser.APPROVE_OPTION) {
                    try {
                        plot.writeChartAsPDF(new FileOutputStream(fileChooser.getSelectedFile()));
                    }
                    catch (Exception err) {
                        err.printStackTrace();
                    }
                }
            }
        });
        JButton exportSVGButton = new JButton("Export SVG");
        //f.getContentPane().add(exportSVGButton, BorderLayout.SOUTH);
        exportSVGButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                JFileChooser fileChooser = new JFileChooser();
                if (fileChooser.showSaveDialog(f) == JFileChooser.APPROVE_OPTION) {
                    try {
                        plot.writeChartAsSVG(new FileOutputStream(fileChooser.getSelectedFile()));
                    }
                    catch (Exception err) {
                        err.printStackTrace();
                    }
                }
            }
        });
        f.pack();

        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                RefineryUtilities.centerFrameOnScreen(f);
                f.setVisible(true);
            }
        });
    }

}
