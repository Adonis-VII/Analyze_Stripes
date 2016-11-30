//Copyright 2016 Ryan P. Sullivan.
//Copyright 2013 Justin R. Bickford.
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”),
//copy, modify, merge, publish, or otherwise alter this software for educational or academic purposes subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//The copyright holders of other software modified and included in the Software retain their rights and the licenses on that software should not be removed.
//Cite the authors (above) of this plugin in any publication that relies on the Software. Also cite those projects on which the Software relies when applicable.
//
//THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
//WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.macro.Interpreter;
import ij.plugin.*;
import ij.plugin.frame.RoiManager;
import ij.plugin.ContrastEnhancer;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.util.Tools;
import ij.gui.PolygonRoi;
import ij.gui.PointRoi;
import ij.gui.WaitForUserDialog;
import ij.gui.Roi;

public class Analyze_Stripes implements PlugIn {

protected ImagePlus Image, ActiveImage;
ImageProcessor ActiveImageProcessor;
ImageConverter ActiveImageConverter;
ResultsTable rt;
double averageangle, linewidth, RMS_edge_roughness, pkpk_edge_roughness;
int nResults;
int similarity = 5; // controls how wide to cast the net of similarity
int options, measurements;
List<Double> group1, group2;
double[][] Data;
String unit = "pixels";
RoiManager ActiveROIManager;
Roi[] RoiArray, RoiArray1, RoiArray2;
Roi Roi1, Roi2;
boolean userInput = false;

	public void run(String arg) {
		Image=IJ.getImage();//Opened Image
		ActiveImage=Image.duplicate();//Image on which all the processing is done
		ActiveImageProcessor = ActiveImage.getProcessor();
		ActiveImageConverter = new ImageConverter(ActiveImage); 
		unit = ActiveImage.getCalibration().getUnits();
		ActiveROIManager = new RoiManager(true);
		rt = new ResultsTable();
		options = ParticleAnalyzer.CLEAR_WORKSHEET + ParticleAnalyzer.ADD_TO_MANAGER;
		measurements = Measurements.AREA + Measurements.MEAN + Measurements.STD_DEV +  Measurements.MIN_MAX + Measurements.CENTROID + Measurements.ELLIPSE;
		IJ.run("Set Measurements...", "area mean standard min centroid fit redirect=None decimal=3");
		//Convert image to binary edges
		findEdgesAndThreshold();
		//measure edges
		ParticleAnalyzer ActiveParticleAnalyzer = new ParticleAnalyzer(options,measurements,rt,0,Double.POSITIVE_INFINITY,0.00,0.10);
		ActiveImageConverter.convertToGray8();
		ActiveParticleAnalyzer.setRoiManager(ActiveROIManager);
		ActiveParticleAnalyzer.analyze(ActiveImage);
		//Apply rotation to edges
		calculateAngleAndRotate();
		//collect Data
		ActiveROIManager.runCommand("Measure");
		setMeasurementArray();
		sortDataByMean();
		mergeROIs();
		analyzeData();
		ActiveROIManager.moveRoisToOverlay(Image);//Overlays Edges on original image. Stop prompt from next line and keep the Edges
		if (linewidth==0){
			IJ.log("Failed to find distinct edges.");
			return;
		}
		outputToResults();
		rt.show("Results");
		return;
	}
	//Finds the gradient^4 of the image
	//converts to pixels > mean*8 to binary and skeletonizes
	public void findEdgesAndThreshold() {
		ActiveImageConverter.convertToGray16();
		ContrastEnhancer enhancer = new ContrastEnhancer();
		enhancer.equalize(ActiveImage);
		ActiveImageProcessor = ActiveImage.getProcessor();
		ActiveImageProcessor.findEdges();
		ActiveImageProcessor.sqr(); // turns the gradient into gradient squared
		ActiveImageProcessor.sqr(); // further enhances the good edges
		enhancer.equalize(ActiveImage);
		double mean = ActiveImage.getStatistics().mean;
		double max = ActiveImage.getStatistics().max;
		ActiveImageProcessor.setThreshold(mean*8,max,ActiveImageProcessor.OVER_UNDER_LUT);
		if(userInput)
		{
			IJ.run(ActiveImage,"Threshold...","");
			new WaitForUserDialog("OK").show();
			if (WindowManager.getWindow("Threshold")!= null){
				IJ.selectWindow("Threshold");
				IJ.run("Close");  
			}
		}
		IJ.run(ActiveImage,"Convert to Mask","");
		IJ.run(ActiveImage,"Skeletonize","");
		return;
	}
	//Calculates the area weighted angle and rotates the data using v = x*sin(angle) + y*cos(angle)
	public void calculateAngleAndRotate(){
		double cumarea = 0;
		double sina = 0;
		double cosa = 0;
		double area, angle;
		for (int n = 0; n < rt.size();n++) {
			area = rt.getValueAsDouble(rt.getColumnIndex("Area"),n);
			angle = 2*rt.getValueAsDouble(rt.getColumnIndex("Angle"),n);
			sina = sina + area*Math.sin(angle*Math.PI/180);
			cosa = cosa + area*Math.cos(angle*Math.PI/180);
			cumarea = cumarea+area;
		}
		averageangle = Math.abs(0.5*(180/Math.PI)*Math.atan2(sina/cumarea,cosa/cumarea)); // this is the area weighted average angle
		// rotate the data 
		IJ.run(ActiveImage,"Select All","");
		ActiveImageConverter.convertToGray32();
		IJ.run(ActiveImage,"Macro...", "code=[v= x*sin(PI/180*"+averageangle+")+y*cos(PI/180*"+averageangle+")]");
		return;
	}
	//Grabs the mean, standard deviation, max, and min values of the edges and puts them in Data[][]
	public void setMeasurementArray(){
		nResults = rt.size();
		Data = new double[4][nResults];
		for(int i =1; i<nResults; i++) {
			Data[0][i]=i;
			Data[1][i]=rt.getValueAsDouble(rt.getColumnIndex("Mean"),i); 
			Data[2][i]=rt.getValueAsDouble(rt.getColumnIndex("StdDev"),i);
			Data[3][i]=rt.getValueAsDouble(rt.getColumnIndex("Max"),i)-rt.getValueAsDouble(rt.getColumnIndex("Min"),i);
		}
		return;
	}
	//Sort Data[][] by mean
	public void sortDataByMean() {
		for(int i = 1; i<nResults-1; i++){
			double tempIndex = Data[0][i];
			double tempMean = Data[1][i];
			double tempStdDev = Data[2][i];
			double tempextremea = Data[3][i];
			int j = i - 1;
			while (j>=0 && Data[1][j] > tempMean){
				Data[0][j+1]=Data[0][j];
				Data[1][j+1]=Data[1][j];
				Data[2][j+1]=Data[2][j];
				Data[3][j+1]=Data[3][j];
			}
			Data[0][j+1]=tempIndex;
			Data[1][j+1]=tempMean;
			Data[2][j+1]=tempStdDev;
			Data[3][j+1]=tempextremea;
		}
		return;
	}
	//Filter and Merge ROIs (edges) to only 2 outer most edges
	public void mergeROIs() {
		group1 = new ArrayList<Double>();
		group1.add(Data[0][0]);
		for (int i = 1; i < (nResults-1);i++) {
			if ((Data[1][i]-similarity*Data[2][i]) < (Data[1][0]+similarity*Data[2][0]) ) { // if the lines are similar to the minimum line, add them to group1
				group1.add(Data[0][i]);//addes new element to group 1 with value of roiindex[i]
			}
		}
		int[] group1a = new int[group1.size()];
		for(int i =0; i < group1.size(); i++){
			double temp = group1.get(i);
			group1a[i] = (int)temp;
		}
		group2 = new ArrayList<Double>();
		group2.add(Data[0][nResults-1]);
		for (int i=(nResults-2);i>1;i--) {
			if ((Data[1][i]+similarity*Data[2][i]) > (Data[1][nResults-1]-similarity*Data[2][nResults-1]) ) { 
			// if the lines are similar to the maximum line, add them to group2
			group2.add(0,Data[0][i]);
			}
		}
		int[] group2a = new int[group2.size()];
		for(int i =0; i < group2.size(); i++){
			double temp = group2.get(i);
			group2a[i] = (int)temp;
		}
		int count;
		//IJ.run("Select None");
		ActiveROIManager.deselect();
		ActiveROIManager.setSelectedIndexes(group1a);
		if(group1a.length > 1) {
			ActiveROIManager.runCommand(ActiveImage,"Combine"); 
			//Creates new ROI that is combination of group 2
			count = ActiveROIManager.getCount();
			Roi1 = ActiveROIManager.getRoi(count-1);//Selects the combined group1 ROI
		}
		else{
			Roi1 = ActiveROIManager.getRoi(group1a[0]);
		}
		ActiveROIManager.setSelectedIndexes(group2a);
		if(group2a.length > 1) {
			ActiveROIManager.runCommand(ActiveImage,"Combine"); 
			//Creates new ROI that is combination of group 2
			count = ActiveROIManager.getCount();
			Roi2 = ActiveROIManager.getRoi(count-1); //Selects the combined group2 ROI
		}
		else{
			Roi2 = ActiveROIManager.getRoi(group2a[0]);
		}
		ActiveROIManager.reset();
		rt.reset();
		ActiveROIManager.add(ActiveImage,Roi1,0);
		ActiveROIManager.select(ActiveImage,0);
		Analyzer ActiveAnalyzer = new Analyzer(ActiveImage,measurements,rt);
		ActiveAnalyzer.measure();
		ActiveROIManager.add(ActiveImage,Roi2,1);
		ActiveROIManager.select(ActiveImage,1);
		ActiveAnalyzer.measure();
		ActiveROIManager.runCommand(Image,"Show All without labels");// removes the labels, which a user may find confusing in the context of this macro
		return;
	}
	//Calculate Linewidth, RMS edge roughness, and PkPk edge roughness
	public void analyzeData(){
		linewidth = rt.getValue("Mean", 1) - rt.getValue("Mean", 0);
		RMS_edge_roughness = Math.sqrt((Math.pow(rt.getValue("StdDev", 1), 2) + Math.pow(rt.getValue("StdDev", 0), 2))/2);  // this is the correct form
		if(rt.getValue("Max",1) - rt.getValue("Min", 1) > rt.getValue("Max", 0) - rt.getValue("Min", 0)){
			pkpk_edge_roughness = rt.getValue("Max", 1) – rt.getValue("Min", 1);
		}
		else{
			pkpk_edge_roughness = (rt.getValue("Max", 0) – rt.getValue("Min", 0));
		}
		return;
	}
	public void setUserInput(boolean userInput){
		this.userInput = userInput;
	}
	public void outputToResults(){
		rt = new ResultsTable();
		rt.incrementCounter();
		rt.addLabel("Average Angle");
		rt.setValue(1,0,averageangle);
		rt.setValue(2,0,"degrees");
		rt.incrementCounter();
		rt.addLabel("Line width");
		rt.setValue(1,1,linewidth);
		rt.setValue(2,1,unit);
		rt.incrementCounter();
		rt.addLabel("RMS Edge Roughness (Rq)");
		rt.setValue(1,2,RMS_edge_roughness);
		rt.setValue(2,2,unit);
		rt.incrementCounter();
		rt.addLabel("Peak-to-Peak Edge Roughness (Rt)");
		rt.setValue(1,3,pkpk_edge_roughness);
		rt.setValue(2,3,unit);
		return;
	}