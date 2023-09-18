package com.insilico.application.insilicotools.chart;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.util.ImageUtil;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYShapeAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
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
import java.awt.geom.Rectangle2D;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Created by jfeng1 on 6/16/16.
 */
public class PrincipalMomentOfInertiaPlot extends JPanel{
    private final static double[][] PMIPct = {{0.0,1.0},{0.5,0.5},{1.0,1.0},{0.0,1.0}};

    private final Shape pmiShape;

    private final XYShapeAnnotation pmiAnnotation;
    private final XYTextAnnotation rodAnnotation =  new XYTextAnnotation("Rod-like",0.0,1.02);
    private final XYTextAnnotation discAnnotation = new XYTextAnnotation("Disc-like",0.5,0.48);
    private final XYTextAnnotation sphereAnnotation = new XYTextAnnotation("Sphere-like",1.0,1.02);

    private double minAnnotationX = Double.MAX_VALUE;
    private double minAnnotationY = Double.MAX_VALUE;
    private double maxAnnotationX = Double.MIN_VALUE;
    private double maxAnnotationY = Double.MIN_VALUE;

    private Range defaultDomainRange;
    private Range defaultRangeRange;

    private final XYPlot plot;
    private final JFreeChart chart;
    private final Scatterplot scatterplot;
    private final SelectionXYDataset data;


    public PrincipalMomentOfInertiaPlot(String title, String xAxis, String yAxis, boolean createLegend, boolean showToolTip) {
        super(new BorderLayout());

        JLabel titleLabel = new JLabel(title, JLabel.CENTER);
        titleLabel.setFont(titleLabel.getFont().deriveFont(18f));
        add(titleLabel, BorderLayout.NORTH);
        scatterplot = new Scatterplot(null, xAxis, yAxis, createLegend) {
            public void chartMouseClicked(ChartMouseEvent event) {
                super.chartMouseClicked(event);

                Rectangle2D dataArea = getChartRenderingInfo().getPlotInfo().getDataArea();
                Point pt = event.getTrigger().getPoint();

                if (isPointOnEdge(PMIPct, pt, dataArea)) {
                    selectPointsInShape(pmiShape);
                }
            }
        };

        chart = scatterplot.getChart();
        plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        data = (SelectionXYDataset) plot.getDataset();

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
                            PropertyMolecule propertyMolecule = (PropertyMolecule)key;
                            return ImageUtil.smi2ImgUrl(String.format("%s %s",propertyMolecule.getSmiles(),propertyMolecule.getName()));
                        }
                        return null;
                    }
                    catch (Exception err) {
                        err.printStackTrace();
                        return null;
                    }
                }
            });
        }

        final BasicStroke solid = new BasicStroke(1.0f);

        pmiShape = createPolygon(PMIPct);
        pmiAnnotation = new XYShapeAnnotation(pmiShape, solid, Color.black);
        plot.getRenderer().addAnnotation(pmiAnnotation, Layer.BACKGROUND);
        plot.getRenderer().addAnnotation(rodAnnotation, Layer.BACKGROUND);
        plot.getRenderer().addAnnotation(discAnnotation, Layer.BACKGROUND);
        plot.getRenderer().addAnnotation(sphereAnnotation, Layer.BACKGROUND);
        Rectangle2D bounds = pmiShape.getBounds2D();
        minAnnotationX = Math.min(minAnnotationX, bounds.getMinX());
        minAnnotationY = Math.min(minAnnotationY, bounds.getMinY());
        maxAnnotationX = Math.max(maxAnnotationX, bounds.getMaxX());
        maxAnnotationY = Math.max(maxAnnotationY, bounds.getMaxY());

//        Range domainRange = plot.getDomainAxis().getRange();
//        Range rangeRange = plot.getRangeAxis().getRange();

        defaultDomainRange = new Range(minAnnotationX-0.1, maxAnnotationX+0.1);
        defaultRangeRange = new Range(minAnnotationY-0.1, maxAnnotationY+0.1);


        data.addChangeListener(new DatasetChangeListener() {
            public void datasetChanged(DatasetChangeEvent event) {
                Range domainRange = plot.getDomainAxis().getRange();
                Range rangeRange = plot.getRangeAxis().getRange();
                defaultDomainRange = new Range(Math.min(minAnnotationX, domainRange.getLowerBound()), Math.max(maxAnnotationX, domainRange.getUpperBound()));
                defaultRangeRange = new Range(Math.min(minAnnotationY, rangeRange.getLowerBound()), Math.max(maxAnnotationY, rangeRange.getUpperBound()));
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
        scatterplot.getPlot().getRenderer().setSeriesShape(series, new Ellipse2D.Double(-4,-4,8,8));
        return series;
    }

    public void addPoint(int series, double x, double y, PropertyMolecule molecule) {
        scatterplot.addPoint(series, x, y, molecule);
    }


    public void writeChartAsSVG(OutputStream out) throws IOException {
        scatterplot.writeChartAsSVG(out);
    }

    public void writeChartAsPDF(OutputStream out) throws IOException, ScatterplotException {
        scatterplot.writeChartAsPDF(out);
    }

    public void resetRanges(){
        Rectangle2D bounds = pmiShape.getBounds2D();
        minAnnotationX = Math.min(minAnnotationX, bounds.getMinX());
        minAnnotationY = Math.min(minAnnotationY, bounds.getMinY());
        maxAnnotationX = Math.max(maxAnnotationX, bounds.getMaxX());
        maxAnnotationY = Math.max(maxAnnotationY, bounds.getMaxY());

        defaultDomainRange = new Range(minAnnotationX-0.1, maxAnnotationX+0.1);
        defaultRangeRange = new Range(minAnnotationY-0.1, maxAnnotationY+0.1);

        reset();

    }

    public void clear() {
        scatterplot.clear();
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

    public static void main(String[] args) {
        final PrincipalMomentOfInertiaPlot plot = new PrincipalMomentOfInertiaPlot("<html><body><font color=\"#0000FF\">Principal Moment of Inertia</font></body></html>", "NPR1", "NPR2", true, true);

        int series = plot.addSeries("BIO Number");
        for (int i = 0; i < 10; i++) {
            OEGraphMol mol = new OEGraphMol();
            oechem.OEParseSmiles(mol,"CCCCC");
            double abs1 = Math.abs(Math.random());
            double abs2 = 0.5+Math.abs(Math.random())/2;
            plot.addPoint(series, abs1, abs2, new PropertyMolecule(mol));
        }

        final JFrame f = new JFrame("PMI Example");
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
