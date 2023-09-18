package com.insilico.application.insilicotools.util;

import com.google.common.collect.Lists;
import com.google.common.io.CharStreams;
import openeye.oechem.oechem;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.List;

public final class OEChemWebLicenseInstaller {
    private OEChemWebLicenseInstaller() {
    }

    public static void loadOELicenseFromWeb() throws IOException {
//      List<String> lines = CharStreams.readLines(new FileReader(new File("/Users/junfeng/oe_license.txt")));
        List<String> lines = CharStreams.readLines(new InputStreamReader(new URL("http://10.74.2.128:8080/insilico_tools/oe_license.txt").openStream()));
        List<OEProduct> products = Lists.newArrayList();
        for(String line : lines) {
            if(line.trim().isEmpty()) {
                continue;
            }
            String[] tokens = line.split(":");
            if(tokens.length < 2) {
                continue;
            }
            String value = tokens[1].trim();

            if(line.startsWith("#PRODUCT:")) {
                OEProduct product = new OEProduct(value);
                products.add(product);
                continue;
            }
            if(products.isEmpty()) {
                continue;
            }
            OEProduct currentProduct = products.get(products.size() - 1);
            if(line.startsWith("#LICENSEE")) {
                currentProduct.setLicensee(value);
            }

            if(line.startsWith("#SITE")) {
                currentProduct.setSite(value);
            }

            if(line.startsWith("LICENSE_KEY:")) {
                currentProduct.setLicenseKey(value);
            }
        }

        for(OEProduct product : products) {
            product.addLicense();
        }
    }

    static final class OEProduct {
        private String productName;
        private String licenseKey;
        private String site;
        private String licensee;

        public boolean addLicense() {
            if(licenseKey != null && site != null && licensee != null) {
                oechem.OEAddLicenseKey(licenseKey, licensee, site);
                return true;
            }
            return false;
        }

        OEProduct(String productName) {
            this.productName = productName;
        }

        public String getProductName() {
            return productName;
        }

        public void setProductName(String productName) {
            this.productName = productName;
        }

        public String getLicenseKey() {
            return licenseKey;
        }

        public void setLicenseKey(String licenseKey) {
            this.licenseKey = licenseKey;
        }

        public String getSite() {
            return site;
        }

        public void setSite(String site) {
            this.site = site;
        }

        public String getLicensee() {
            return licensee;
        }

        public void setLicensee(String licensee) {
            this.licensee = licensee;
        }
    }
}
