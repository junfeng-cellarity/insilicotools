package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.data.MolProperty;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.util.ExcelPixelUtil;
import com.insilico.application.insilicotools.util.ChemFunc;
import org.apache.poi.hssf.usermodel.*;
import org.apache.poi.ss.usermodel.ClientAnchor;
import org.apache.poi.ss.usermodel.Workbook;

import javax.swing.table.DefaultTableModel;
import java.util.Vector;

/**
 * Created by jfeng1 on 9/8/15.
 */
public class PropertyMolTableModel extends DefaultTableModel {
    Vector<PropertyMolecule> propertyMolecules;
    Vector<String> selectedProperties;
    Vector<String> existingTags;
    boolean showMarked = false;

    public PropertyMolTableModel(Vector<PropertyMolecule> molecules, Vector<String> tags){
        this(molecules,tags,false);
    }

    public PropertyMolTableModel(Vector<PropertyMolecule> molecules, Vector<String> tags, boolean showMarked) {
        propertyMolecules = molecules;
        selectedProperties = new Vector<String>();
        existingTags = tags;
        if(molecules==null){
            propertyMolecules = new Vector<PropertyMolecule>();
            existingTags = new Vector<String>();
        }
        this.showMarked = showMarked;
    }

    public void setSelectedProperties(Vector<String> selectedProperties) {
        if(selectedProperties!=this.selectedProperties) {
            this.selectedProperties.clear();
            for(String p:selectedProperties){
                if(!this.selectedProperties.contains(p)){
                    this.selectedProperties.add(p);
                }
            }
        }
        fireTableStructureChanged();
        fireTableDataChanged();
    }

    public Vector<String> getSelectedProperties() {
        return selectedProperties;
    }

    public Vector<PropertyMolecule> getPropertyMolecules() {
        return propertyMolecules;
    }

    public Vector<String> getExistingTags() {
        return existingTags;
    }

    @Override
    public int getRowCount() {
        return propertyMolecules==null?0:propertyMolecules.size();
    }

    @Override
    public int getColumnCount() {
        if(propertyMolecules!=null) {
            if (propertyMolecules.size() > 0) {
                if(showMarked){
                    return selectedProperties.size() + 4;
                }else {
                    return selectedProperties.size() + 3;
                }
            }
        }
        return 0;
    }

    @Override
    public String getColumnName(int column) {
        if(propertyMolecules!=null) {
            if (propertyMolecules.size() > 0) {
                if (column > 2) {
                    if(showMarked){
                        if(column<getColumnCount()-1){
                            if (selectedProperties.size() > 0) {
                                return selectedProperties.get(column - 3);
                            }
                        }else{
                            return "Marked";
                        }
                    }else {
                        if (selectedProperties.size() > 0) {
                            return selectedProperties.get(column - 3);
                        }
                    }
                } else {
                    if (column == 0) {
                        return "";
                    }
                    if (column == 1) {
                        return "Structure";
                    }
                    if (column == 2) {
                        return "Name";
                    }
                }
            }
        }
        return "";
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        if(column==0){
            return true;
        }else {
            if(showMarked){
                if(column==getColumnCount()-1){
                    return true;
                }else{
                    return false;
                }
            }else {
                return false;
            }
        }
    }



    @Override
    public void setValueAt(Object aValue, int row, int column) {
        if(column==0 && aValue!=null&& aValue instanceof Boolean){
            propertyMolecules.get(row).setIsSelected((Boolean)aValue);
//            fireTableDataChanged();
        }else if(showMarked&&column==getColumnCount()-1&& aValue!=null&& aValue instanceof Boolean){
            propertyMolecules.get(row).setMarked((Boolean)aValue);
//            fireTableDataChanged();
        }
    }

    @Override
    public Object getValueAt(int row, int column) {
        if(propertyMolecules.size()==0){
            return null;
        }
        if(propertyMolecules.size()<=row){
            return null;
        }
        if(column>2) {
            if(!showMarked) {
                return propertyMolecules.get(row).getProperty(getColumnName(column));
            }else{
                if(column==getColumnCount()-1){
                    return propertyMolecules.get(row).isMarked();
                }else{
                    return propertyMolecules.get(row).getProperty(getColumnName(column));
                }
            }
        }else{
            if(column==0){
                return propertyMolecules.get(row).isSelected();
            }
            else if(column==1) {
                return propertyMolecules.get(row);
            }else{
                return propertyMolecules.get(row).getName();
            }
        }
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        if(columnIndex>2){
            if(showMarked){
                if(columnIndex==getColumnCount()-1){
                    return Boolean.class;
                }else{
                    return MolProperty.class;
                }
            }else {
                return MolProperty.class;
            }
        }else{
            if(columnIndex==0){
                return Boolean.class;
            }else if(columnIndex==1){
                return PropertyMolecule.class;
            }else {
                return String.class;
            }
        }
    }

    public void clearProperties(){
        selectedProperties.clear();
        fireTableStructureChanged();
        fireTableDataChanged();
    }

    public HSSFWorkbook saveAsExcelSheet(ProgressReporter progressReporter, boolean saveImage){
        int rowNum = this.getRowCount();
        int colNum = this.getColumnCount();

        HSSFWorkbook excelWorkBook = new HSSFWorkbook();
        HSSFSheet sheet1 = excelWorkBook.createSheet();

        if(saveImage) {
//            sheet1.setColumnWidth(0, ExcelPixelUtil.pixel2WidthUnits(200));
            sheet1.setColumnWidth(0, ExcelPixelUtil.pixel2WidthUnits(350));
        }

        HSSFRow headerRow = sheet1.createRow(0);
        for (int i = 1; i < colNum; i++) {
            HSSFCell headerCell = headerRow.createCell(i-1);
            headerCell.setCellValue(this.getColumnName(i));
        }

        for (int i = 0; i < rowNum; i++) {
            if(progressReporter!=null){
                progressReporter.reportProgress("Progress", 100*i/rowNum);
            }
            HSSFRow row = sheet1.createRow(i + 1);
            row.setHeight(ExcelPixelUtil.pixel2HeightUnits(300));
//            row.setHeight(ExcelPixelUtil.pixel2HeightUnits(150));
            for (int j = 1; j < colNum; j++) {
                HSSFCell cell = row.createCell(j-1);
                Object cellValue = this.getValueAt(i, j);
                if (cellValue == null) {
                    cell.setCellValue("");
                } else {
                    if(cellValue instanceof PropertyMolecule){
                        PropertyMolecule mol = (PropertyMolecule) cellValue;
                        if(saveImage) {
                            byte[] pngBytes = ChemFunc.mol2png(mol,500,500);
                            ClientAnchor imageAnchor = new HSSFClientAnchor();
                            imageAnchor.setCol1(cell.getColumnIndex());
                            imageAnchor.setRow1(cell.getRowIndex());
                            int picture_id = excelWorkBook.addPicture(pngBytes, Workbook.PICTURE_TYPE_PNG);
                            HSSFPatriarch drawing = sheet1.createDrawingPatriarch();
                            HSSFPicture picture = drawing.createPicture(imageAnchor, picture_id);
                            double pic_width = picture.getImageDimension().getWidth();
                            double pic_height = picture.getImageDimension().getHeight();
                            double column_width = ExcelPixelUtil.heightUnits2Pixel(row.getHeight());
                            double row_height = ExcelPixelUtil.widthUnits2Pixel((short)cell.getSheet().getColumnWidth(cell.getColumnIndex()));
                            double pic_dim = Math.max(pic_width, pic_height);
                            double cell_dim = Math.min(column_width,row_height);
                            System.out.println(cell_dim+" "+pic_dim);
                            picture.resize(cell_dim/pic_dim);
                        }else{
                            cell.setCellValue(mol.getSmiles());
                        }
                    }else if(cellValue instanceof MolProperty) {
                        MolProperty prop = (MolProperty)cellValue;
                        if(prop.isNumerical()){
                            cell.setCellValue(prop.getValue());
                        }else{
                            cell.setCellValue(prop.getProperty());
                        }
                    }
                    else{
                        cell.setCellValue(cellValue.toString());
                    }
                }
            }
        }
        return excelWorkBook;
    }



}
