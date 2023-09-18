package com.insilico.application.insilicotools.chart;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Point;
import java.awt.Shape;
import java.awt.event.MouseEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.event.EventListenerList;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.entity.ChartEntity;
import org.jfree.chart.entity.XYItemEntity;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RefineryUtilities;
import org.jfree.util.ShapeUtilities;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

public class Scatterplot extends ChartPanel implements ChartMouseListener {

    private Point2D selectPoint;
    protected SelectionXYDataset data;

    private Rectangle2D.Double selectRectangle;
    private Point2D zoomPoint = null;
    private int zoomTriggerDistance = 10; // default zoom trigger distance from 0.9
    private transient Rectangle2D zoomRectangle = null;

    protected EventListenerList listenerList = new EventListenerList();
    protected XYPlot plot;
    protected JFreeChart chart;
    protected StandardXYItemRenderer renderer;

    public Scatterplot(String title, String xAxis, String yAxis) {
        this(title, xAxis, yAxis, false, false);
    }

    public Scatterplot(String title, String xAxis, String yAxis, boolean createLegend) {
        this(title, xAxis, yAxis, false, false, createLegend, false);
    }

    public Scatterplot(String title, String xAxis, String yAxis, boolean isXLog, boolean isYLog) {
        this(title, xAxis, yAxis, isXLog, isYLog, true, false);
    }

    public Scatterplot(String title, String xAxis, String yAxis, boolean isXLog, boolean isYLog, boolean createLegend) {
        this(title, xAxis, yAxis, isXLog, isYLog, createLegend, false);
    }

    public Scatterplot(String title, String xAxis, String yAxis, boolean isXLog, boolean isYLog, boolean createLegend, boolean createLine) {
        super(ChartFactory.createScatterPlot(title, xAxis, yAxis, new SelectionXYDataset(), PlotOrientation.VERTICAL, createLine ? false : createLegend, false, false));

        chart = getChart();

        plot = chart.getXYPlot();
        data = (SelectionXYDataset) plot.getDataset();

        if (isXLog) {
            LogarithmicAxis xAxisLog = new LogarithmicAxis("Log(x)");
            xAxisLog.setAllowNegativesFlag(false);
            xAxisLog.setLog10TickLabelsFlag(false);
            plot.setDomainAxis(xAxisLog);
        }
        if (isYLog) {
            LogarithmicAxis yAxisLog = new LogarithmicAxis("Log(y)");
            yAxisLog.setAllowNegativesFlag(false);
            yAxisLog.setLog10TickLabelsFlag(false);
            plot.setRangeAxis(yAxisLog);
        }

        if (createLine) {
            plot.setDataset(1, data);
            XYLineAndShapeRenderer lineShapeRenderer = new XYLineAndShapeRenderer() {
                public Paint getItemPaint(int row, int column) {
                    if ((data.getItem(row, column)).isSelected()) {
                        return Color.red;
                    } else {
                        return super.getItemPaint(row, column);
                    }
                }
            };
            lineShapeRenderer.setSeriesPaint(0, Color.blue);
            lineShapeRenderer.setToolTipGenerator(new XYToolTipGenerator() {
                public String generateToolTip(XYDataset xyDataset, int i, int i1) {
                    return String.format("%8.3f,%8.3f", xyDataset.getX(i, i1), xyDataset.getY(i, i1));
                }
            });
            plot.setRenderer(1, lineShapeRenderer);
        }

        chart.setBackgroundPaint(getBackground());
        renderer = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES) {
            public Paint getItemPaint(int row, int column) {
                if ((data.getItem(row, column)).isSelected()) {
                    return Color.red;
                } else {
                    return super.getItemPaint(row, column);
                }
            }
        };
        renderer.setBaseShape(new Ellipse2D.Float(0, 0, 3, 3));
        plot.setRenderer(renderer);

        setPopupMenu(null);
        
        this.setDomainZoomable(false);
        this.setRangeZoomable(false);
        addChartMouseListener(this);
    }

    public void setToolTipGenerator(XYToolTipGenerator generator) {
        plot.getRenderer().setBaseToolTipGenerator(generator);
    }

    public void setBackground(Color bg) {
        super.setBackground(bg);
        if (chart != null) {
            chart.setBackgroundPaint(bg);
        }
    }

    //compatibility with jfreechart 0.9
    private Rectangle2D getScaledDataArea() {
    	return getScreenDataArea();
    }
    
    public void mouseReleased(MouseEvent event) {
        if (SwingUtilities.isRightMouseButton(event)) {
            if (zoomRectangle != null) {
                boolean zoomTrigger1 = Math.abs(event.getX() - zoomPoint.getX()) >= zoomTriggerDistance;
                boolean zoomTrigger2 = Math.abs(event.getY() - zoomPoint.getY()) >= zoomTriggerDistance;
                if (zoomTrigger1 || zoomTrigger2) {
                    double x, y, w, h;
                    Rectangle2D scaledDataArea = getScaledDataArea();

                    x = zoomPoint.getX();
                    y = zoomPoint.getY();
                    w = Math.min(zoomRectangle.getWidth(), scaledDataArea.getMaxX() - zoomPoint.getX());
                    h = Math.min(zoomRectangle.getHeight(), scaledDataArea.getMaxY() - zoomPoint.getY());

                    if (((event.getX() < zoomPoint.getX())) || ((event.getY() < zoomPoint.getY()))) {
                    	//plot.zoomHorizontalAxes(Math.abs(w / 10f));
                    	//plot.zoomVerticalAxes(Math.abs(h / 10f));
                    	plot.zoomRangeAxes(Math.abs(w / 10f), null, null);
                    	plot.zoomDomainAxes(Math.abs(h / 10f), null, null);
                    } else {
                        Rectangle2D zoomArea = new Rectangle2D.Double(x, y, w, h);
                        zoom(zoomArea);
                    }
                    zoomPoint = null;
                    zoomRectangle = null;
                } else {
                    Graphics2D g2 = (Graphics2D) getGraphics();
                    g2.setXORMode(java.awt.Color.gray);
                    g2.draw(zoomRectangle);
                    g2.dispose();
                    zoomRectangle = null;
                }
            }
        } else {
            if (SwingUtilities.isLeftMouseButton(event)) {
                if (selectRectangle != null) {
                    Graphics2D g2 = (Graphics2D) getGraphics();
                    g2.setXORMode(java.awt.Color.gray);
                    g2.draw(selectRectangle);
                    g2.dispose();

                    Point2D bottomLeftPoint = translateScreenToJava2D(new Point((int) selectRectangle.getX(), (int) selectRectangle.getY()));
                    Point2D topRightPoint = translateScreenToJava2D(new Point((int) (selectRectangle.getX() + selectRectangle.width), (int) (selectRectangle.getY() + selectRectangle.height)));

                    XYPlot plot = getChart().getXYPlot();
                    Rectangle2D dataArea = getChartRenderingInfo().getPlotInfo().getDataArea();
                    double dataSetLeftX = plot.getDomainAxis().java2DToValue(bottomLeftPoint.getX(), dataArea, plot.getDomainAxisEdge());
                    double dataSetBottomY = plot.getRangeAxis().java2DToValue(bottomLeftPoint.getY(), dataArea, plot.getRangeAxisEdge());

                    double dataSetRightX = plot.getDomainAxis().java2DToValue(topRightPoint.getX(), dataArea, plot.getDomainAxisEdge());
                    double dataSetTopY = plot.getRangeAxis().java2DToValue(topRightPoint.getY(), dataArea, plot.getRangeAxisEdge());

                    Rectangle2D.Double r = new Rectangle2D.Double(dataSetLeftX, dataSetTopY, Math.abs(dataSetRightX - dataSetLeftX), Math.abs(dataSetTopY - dataSetBottomY));

                    int seriesCount = data.getSeriesCount();
                    for (int i = 0; i < seriesCount; i++) {
                        int itemCount = data.getItemCount(i);
                        for (int j = 0; j < itemCount; j++) {
                            SelectionXYDataItem item = data.getItem(i, j);
                            if (r.contains(item.getX().doubleValue(), item.getY().doubleValue())) {
                                item.setSelected(true);
                            }
                        }
                    }

                    fireSelectionChanged();

                    repaint();
                    selectRectangle = null;
                }
            }
        }
    }

    public void mouseDragged(MouseEvent event) {
        if (SwingUtilities.isRightMouseButton(event)) {
            Graphics2D g2 = (Graphics2D) getGraphics();

            // use XOR to erase the previous zoom rectangle (if any)...
            g2.setXORMode(java.awt.Color.gray);
            if (zoomRectangle != null) {
                g2.draw(zoomRectangle);
            }

            Rectangle2D scaledDataArea = getScaledDataArea();
            // selected rectangle shouldn't extend outside the com.My.application.insilicotools.data area...
            double xmax = Math.min(event.getX(), scaledDataArea.getMaxX());
            double ymax = Math.min(event.getY(), scaledDataArea.getMaxY());
            zoomRectangle = new Rectangle2D.Double(zoomPoint.getX(), zoomPoint.getY(), xmax - zoomPoint.getX(), ymax - zoomPoint.getY());

            if (zoomRectangle != null) {
                g2.draw(zoomRectangle);
            }
            g2.dispose();
        } else {
            if (SwingUtilities.isLeftMouseButton(event)) {
                Graphics2D g2 = (Graphics2D) getGraphics();

                // use XOR to erase the previous zoom rectangle (if any)...
                g2.setXORMode(Color.gray);
                if (selectRectangle != null) {
                    g2.draw(selectRectangle);
                }

                Rectangle2D scaledDataArea = getScaledDataArea();
                double xmax = Math.min(Math.max(selectPoint.getX(), event.getX()), scaledDataArea.getMaxX());
                double ymax = Math.min(Math.max(selectPoint.getY(), event.getY()), scaledDataArea.getMaxY());
                double x = Math.max(Math.min(selectPoint.getX(), event.getX()), scaledDataArea.getMinX());
                double y = Math.max(Math.min(selectPoint.getY(), event.getY()), scaledDataArea.getMinY());
                selectRectangle = new Rectangle2D.Double(x, y, Math.abs(xmax - x), Math.abs(ymax - y));

                if (selectRectangle != null) {
                    g2.draw(selectRectangle);
                }
                g2.dispose();
            }
        }
    }

    public void mousePressed(MouseEvent event) {
        if (SwingUtilities.isRightMouseButton(event)) {
            if (zoomRectangle == null) {
                zoomPoint = ShapeUtilities.getPointInRectangle(event.getX(), event.getY(), getScaledDataArea());
            }
        } else {
            if (SwingUtilities.isLeftMouseButton(event)) {
                if (selectRectangle == null) {
                    if (!event.isShiftDown()) {
                        data.clearSelection();
                    }
                    selectPoint = ShapeUtilities.getPointInRectangle(event.getX(), event.getY(), getScaledDataArea());
                }
            }
        }
    }

    public void chartMouseClicked(ChartMouseEvent event) {
        if (!event.getTrigger().isShiftDown()) {
            data.clearSelection();
            fireSelectionChanged();
        }

        ChartEntity entity = event.getEntity();
        if (entity != null && entity instanceof XYItemEntity) {
            data.getItem(((XYItemEntity) entity).getSeriesIndex(), ((XYItemEntity) entity).getItem()).setSelected(true);
            fireSelectionChanged();
        }
    }

    public void chartMouseMoved(ChartMouseEvent event) {
    }

    public void addAnnotation(XYAnnotation xyAnnotation) {
        plot.addAnnotation(xyAnnotation);
    }

    public int addSeries(String name) {
        return addSeries(name, Color.blue);
    }

    public void clearAnnotations() {
        plot.clearAnnotations();
    }

    public Object[] getSelectedPoints() {
        SelectionXYDataItem[] items = data.getSelectedPoints();
        Object[] keys = new Object[items.length];
        for (int i = 0; i < items.length; i++) {
            keys[i] = items[i].getKey();
        }
        return keys;
    }

    public int addSeries(String name, Color color) {
        int i = data.addSeries(name);
        renderer.setSeriesPaint(i, color);
        return i;
    }

    public void removeSeries(int series) {
        data.deleteSeries(series);
    }

    public void addPoint(int series, double x, double y, Object key) {
        data.addItem(series, new SelectionXYDataItem(x, y, key));
    }

    public void addSelectionListener(SelectionListener l) {
        listenerList.add(SelectionListener.class, l);
    }

    public void removeSelectionListener(SelectionListener l) {
        listenerList.remove(SelectionListener.class, l);
    }

    public SelectionListener[] getSelectionListeners() {
        return (SelectionListener[]) listenerList.getListeners(SelectionListener.class);
    }

    public boolean selectPoint(int series, Object key) {
        return data.selectPoint(series, key);
    }

    public void selectAll() {
        data.selectAll();
    }

    public void deselectAll() {
        data.deselectAll();
    }

    public void selectAll(int series) {
        data.selectAll(series);
    }

    public void clearSelection() {
        data.clearSelection();
    }

    public void zoomIn() {
        zoomInBoth(0, 0);
    }

    public void zoomOut() {
        zoomOutBoth(0, 0);
    }

    public void reset() {
    	if (this.getChart() != null && this.getChart().getXYPlot() != null) {
    		this.getChart().getXYPlot().getRangeAxis().setAutoRange(true);
    		this.getChart().getXYPlot().getDomainAxis().setAutoRange(true);
    	}
    }

    public void fireSelectionChanged() {
        // Guaranteed to return a non-null array
        Object[] listeners = listenerList.getListenerList();
        // Process the listeners last to first, notifying
        // those that are interested in this event
        SelectionEvent evt = null;
        for (int i = listeners.length - 2; i >= 0; i -= 2) {
            if (listeners[i] == SelectionListener.class) {
                if (evt == null) {
                    evt = new SelectionEvent(this, data.getSelectedPoints());
                }
                ((SelectionListener) listeners[i + 1]).selectionChanged(evt);
            }
        }
        this.getChart().fireChartChanged();
    }

    public void writeChartAsJpeg(OutputStream out, int width, int height) throws IOException {
        Dimension d = getSize();
        ChartUtilities.writeChartAsJPEG(out, chart, d.width, d.height);
    }

    public void writeChartAsJpeg(OutputStream out) throws IOException {
        Dimension d = getSize();
        writeChartAsJpeg(out, d.width, d.height);
    }

    public void writeChartAsPNG(OutputStream out, int width, int height) throws IOException {
        Dimension d = getSize();
        writeChartAsPNG(out, d.width, d.height);
    }

    public void writeChartAsPNG(OutputStream out) throws IOException {
        Dimension d = getSize();
        ChartUtilities.writeChartAsPNG(out, chart, d.width, d.height);
    }

    public void writeChartAsSVG(OutputStream out, int width, int height) throws IOException {
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        Document document = domImpl.createDocument(null, "svg", null);
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
        svgGenerator.getGeneratorContext().setPrecision(6);
        Dimension d = getSize();
        chart.draw(svgGenerator, new Rectangle2D.Double(0, 0, d.width, d.height), null);
        boolean useCSS = true;
        Writer w = new OutputStreamWriter(out, "UTF-8");
        svgGenerator.stream(w, useCSS);
        w.flush();
        w.close();
    }

    public void writeChartAsSVG(OutputStream out) throws IOException {
        Dimension d = getSize();
        writeChartAsSVG(out, d.width, d.height);
    }

    public void writeChartAsPDF(OutputStream out, int width, int height) throws IOException, ScatterplotException {
        Rectangle pagesize = new Rectangle(width, height);
        com.itextpdf.text.Document document = new com.itextpdf.text.Document(pagesize, 50, 50, 50, 50);
        PdfWriter writer = null;
        try {
            writer = PdfWriter.getInstance(document, out);
        }
        catch (DocumentException err) {
            throw new ScatterplotException(err);
        }
        document.addAuthor("Cellarity");
        document.open();
        PdfContentByte cb = writer.getDirectContent();
        PdfTemplate tp = cb.createTemplate(width, height);
        Graphics2D g2 = tp.createGraphics(width, height, new DefaultFontMapper());
        Rectangle2D r2D = new Rectangle2D.Double(0, 0, width, height);
        chart.draw(g2, r2D);
        g2.dispose();
        cb.addTemplate(tp, 0, 0);
        document.close();
        out.flush();
        out.close();
    }

    public void writeChartAsPDF(OutputStream out) throws IOException, ScatterplotException {
        Dimension d = getSize();
        writeChartAsPDF(out, d.width, d.height);
    }

    public void print() {
        createChartPrintJob();
    }

    public void clear() {
        data.clear();
    }

    public XYPlot getPlot() {
        return plot;
    }

    public void selectPointsInShape(Shape s) {
        int seriesCount = data.getSeriesCount();
        for (int i = 0; i < seriesCount; i++) {
            int itemCount = data.getItemCount(i);
            for (int j = 0; j < itemCount; j++) {
                SelectionXYDataItem item = data.getItem(i, j);
                if (s.contains(item.getX().doubleValue(), item.getY().doubleValue())) {
                    item.setSelected(true);
                }
            }
        }

        fireSelectionChanged();

        repaint();
    }

    public boolean isPointOnEdge(double[][] shape, Point p, Rectangle2D dataArea) {
        for (int i = 0; i < shape.length - 1; i++) {
            if (new Line2D.Double(convertPoint(shape[i][0], shape[i][1], dataArea), convertPoint(shape[i + 1][0], shape[i + 1][1], dataArea)).ptSegDist(p) < 7) {
                return true;
            }
        }
        if (new Line2D.Double(convertPoint(shape[0][0], shape[0][1], dataArea), convertPoint(shape[shape.length - 1][0], shape[shape.length - 1][1], dataArea)).ptSegDist(p) < 7) {
            return true;
        } else {
            return false;
        }
    }

    private Point2D convertPoint(double x, double y, Rectangle2D dataArea) {
        return new Point2D.Double(plot.getDomainAxis().valueToJava2D(x, dataArea, plot.getDomainAxisEdge()), plot.getRangeAxis().valueToJava2D(y, dataArea, plot.getRangeAxisEdge()));
    }

    public static void main(String[] args) {
        final Scatterplot plot = new Scatterplot("BBB-ABS", "PSA", "ClogP", false, false);

        int series = plot.addSeries("Bios");
        for (int i = 0; i < 5; i++) {
            plot.addPoint(series, 10*i, 5*i, "CCCCC");
        }

        series = plot.addSeries("Other", Color.yellow);
        for (int i = 0; i < 5; i++) {
            plot.addPoint(series, i * 100, i * 100, "CCCCC");
        }

        final JFrame f = new JFrame("Scatterplot Example");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setContentPane(plot);
        f.pack();

        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                RefineryUtilities.centerFrameOnScreen(f);
                f.setVisible(true);
                //plot.format();
            }
        });
    }
}