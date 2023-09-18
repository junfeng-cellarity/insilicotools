package com.insilico.application.insilicotools.gui.util;

import com.schrodinger.jymol.JyMol;

import java.awt.*;

public class JyMolUtilities {
	public static final String LINE_STYLE = "lines";
	public static final String STICK_STYLE = "sticks";
	public static final String CPK_STYLE = "spheres";
	public static final String MESH_STYLE = "mesh";
	public static final String SURFACE_STYLE = "surface";
	public static final String CARTOON_STYLE = "cartoon";
	public static final String PUTTY_STYLE = "putty";
	public static final String COLOR_B_FACTOR = "B Factor";
	public static final String COLOR_SECONDARY_STRUCTURE = "Secondary Structure";
	public static final String COLOR_EP = "EP";
	public static String[] colors = new String[]{"green","aquamarine","bluewhite","brightorange","deepolive","deepsalmon","hotpink","magenta","yellow","violet"};


	public static void setHighQuality(JyMol jymol) {
		if (jymol != null) {
			jymol.cmd.set("valence", "1");
			jymol.cmd.set("auto_show_nonbonded", "1");
			jymol.cmd.set("h_bond_from_proton", "0");
//			jymol.cmd.set("antialias", "1", "");
//			jymol.cmd.set("stick_quality", "5", "");
//			jymol.cmd.set("stick_nub", "0.7", "");
//			jymol.cmd.set("stick_ball", "on", "");
//			jymol.cmd.set("stick_ball_ratio", "1.0", "");
//			jymol.cmd.set("stick_overlap", "0.2", "");
//			jymol.cmd.set("surface_best", "0.25", "");
//			jymol.cmd.set("surface_poor", "0.85", "");
			jymol.cmd.set("surface_quality", "2.0", "");
//			jymol.cmd.set("surface_debug", "0", "");
//			jymol.cmd.set("surface_trim_cutoff", "0.2", "");
//			jymol.cmd.set("surface_trim_factor", "2.0", "");
//			jymol.cmd.set("sphere_quality", "5.0", "");
//			jymol.cmd.set("cgo_sphere_quality", "10.0", "");
//			jymol.cmd.set("cgo_transparency", "0.2");
//			jymol.cmd.set("two_sided_lighting", "1");
//            jymol.cmd.set("nonbonded_size","1");
			//jymol.cmd.set("selection_width","17");
		}
	}

	public static void setUltraHighQuality(JyMol jymol) {
		if (jymol != null) {
			jymol.cmd.set("valence", "1");
			jymol.cmd.set("antialias", "1", "");
			jymol.cmd.set("stick_quality", "15", "");
			jymol.cmd.set("stick_nub", "0.7", "");
			jymol.cmd.set("stick_ball", "on", "");
			jymol.cmd.set("stick_ball_ratio", "1.0", "");
			jymol.cmd.set("stick_overlap", "0.2", "");
			jymol.cmd.set("surface_best", "0.25", "");
			jymol.cmd.set("surface_poor", "0.85", "");
			jymol.cmd.set("surface_quality", "3.0", "");
			jymol.cmd.set("surface_debug", "0", "");
			jymol.cmd.set("surface_trim_cutoff", "0.2", "");
			jymol.cmd.set("surface_trim_factor", "2.0", "");
			jymol.cmd.set("sphere_quality", "5.0", "");
			jymol.cmd.set("cgo_sphere_quality", "10.0", "");
			jymol.cmd.set("cgo_transparency", "0.2");
			jymol.cmd.set("two_sided_lighting", "1");
			jymol.cmd.set("h_bond_from_proton", "0");
//            jymol.cmd.set("nonbonded_size","1");
			jymol.cmd.set("auto_show_nonbonded", "1");
			//jymol.cmd.set("selection_width","17");
		}
	}

	public static void label(JyMol jymol, String selection, String label) {
		if (selection != null && label != null) {
			if (label.equals("")) {
				jymol.cmd.label(selection, "");
			} else {
				String formattedLabel = String.format("%s", label);
				jymol.cmd.label(selection, formattedLabel);
			}
		}
	}

	public static String convertColorRGB(float[] rgbf) {
		String separator = ",";
		if (rgbf == null || rgbf.length != 3) {
			rgbf = Color.WHITE.getRGBColorComponents(null);
		}
		return String.valueOf(rgbf[0]) + separator + rgbf[1] + separator + rgbf[2];
	}

	public static void setMolColor(String colorName, String selectionName, JyMol jymol) {
		jymol.cmd.color(colorName, String.format("(%s and e. c)", selectionName));
	}

	public static void hide_all_hydrogens(JyMol jymol) {
		jymol.cmd.hide("everything", "h.");
	}

	public static void hide_polar_hydrogens(JyMol jymol) {
		jymol.cmd.hide("everything", "h. and (elem c extend 1)");
	}
}
