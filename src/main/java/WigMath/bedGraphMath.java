package WigMath;
/*
 * Copyright (C) 2019  Evan Tarbell and Tao Liu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Max;

import Node.TagNode;

import static HMMR_ATAC.Utils.toMap;

public class bedGraphMath {
	private final HashMap<String,ArrayList<TagNode>> bedgraph;
	private double mean;
	private double std;
	
	public bedGraphMath(ArrayList<TagNode> bdg){
		bedgraph = toMap(bdg);
		setMeanAndStd();
	}


	public double getMean(){return mean;}
	public double getSTD(){return std;}
	public HashMap<String,ArrayList<TagNode>> getMappedBedgraph(){return bedgraph;}

	public static double[] createFilter(int stdev){
		// Use a window size equal to +/- 3 SD's
		double[] filter = new double[6 * stdev + 1];
		double sum = 0;
		for (int i = 0; i < filter.length; i++) {
			double x = (i - 3);
			double value = Math.exp((x * x) / -2);
			filter[i] = value;
			sum += value;
		}
		// Normalize so that the filter is area-preserving (has total area = 1)
		for (int i = 0; i < filter.length; i++) {
			filter[i] /= sum;
		}
		return filter;
	}

	public static void setSmooth(int stdev, TagNode tag, List<TagNode> overlaps, double[] filter) {
		int start = tag.getStart();
		int stop = tag.getStop();

		if (stop - start <= 0){
			return;
		}

		int paddedStart = start - (3 * stdev);
		int paddedStop = stop + (3 * stdev);
		double[] data = new double[paddedStop - paddedStart];

		for (TagNode node2 : overlaps) {
			double value = node2.getScore2();
			for (int i = Math.max(node2.getStart(), paddedStart); i < Math.min(node2.getStop(), paddedStop); i++) {
				data[i - paddedStart] = value;
			}
		}

		// Convolve the data with the filter and find summit
		int largestIndex = -1;
		double tempMax = 0.0;
		for (int i = 0; i < stop - start; i++) {
			double smoothed = 0.0;
			for (int j = 0; j < filter.length; j++) {
				smoothed += data[i + j] * filter[j];
			}
			if (smoothed > tempMax) {
				tempMax = smoothed;
				largestIndex = i;
			}
		}
		if (largestIndex >= 0){
			tag.setSummit(new TagNode(tag.getChrom(), start + largestIndex, start + largestIndex + 1));
		}
	}
	public static void set(TagNode tag, List<TagNode> overlaps, double[] result){
		Max m = new Max();
		Mean mu = new Mean();
		ArrayList<Double> values = new ArrayList<>();
		for (TagNode node2 : overlaps) {
			double value = node2.getScore2();
			m.increment(value);
			for (int i = node2.getStart(); i < node2.getStop(); i++) {
				if (i >= tag.getStart() && i < tag.getStop()) {
					values.add(value);
					mu.increment(value);

				}
			}
		}

		if (values.size() > 0) {
			result[0] = mu.getResult();
			result[1] = values.get(values.size() / 2);
			result[2] = m.getResult();
		} else {
			result[0] = 0;
			result[1] = 0;
			result[2] = 0;
		}
	}
	public ArrayList<TagNode> getAboveZscore(double z){
		ArrayList<TagNode> results = new ArrayList<>();
		for (String chr : bedgraph.keySet()){
			ArrayList<TagNode> inTemp = bedgraph.get(chr);
			inTemp.sort(TagNode.basepairComparator);
			for (TagNode tagNode : inTemp) {
				double value = tagNode.getScore2();
				if (((value - mean) / std) >= z) {
					results.add(tagNode);
				}
			}
		}
		return results;
	}
	
	public ArrayList<TagNode> getBetweenRanges(double upper,double lower){
		ArrayList<TagNode> results = new ArrayList<>();
		for (String chr : bedgraph.keySet()){
			
			ArrayList<TagNode> inTemp = bedgraph.get(chr);
			inTemp.sort(TagNode.basepairComparator);
			for (TagNode tagNode : inTemp) {
				double value = tagNode.getScore2();
				if ((value / mean) >= lower && (value / mean) <= upper) {
					results.add(tagNode);
				}
			}
		}
		return results;
	}
	
	private void setMeanAndStd(){
		Mean mu = new Mean();
		StandardDeviation dev = new StandardDeviation();
		for (String chr : bedgraph.keySet()){
			
			ArrayList<TagNode> inTemp = bedgraph.get(chr);
			inTemp.sort(TagNode.basepairComparator);
			for (TagNode tagNode : inTemp) {
				int length = tagNode.getLength();
				double value = tagNode.getScore2();
				for (int a = 0; a < length; a++) {
					mu.increment(value);
					dev.increment(value);
				}
			}
		}
		mean = mu.getResult();
		std = dev.getResult();
	}
}


