package HMMR_ATAC;

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

import Node.TagNode;

public class ExtendBed {
	
	private final ArrayList<TagNode> input;
	private final int extSize;
	private final ArrayList<TagNode> output;
	
	/**
	 * Constructor for creating ExtendBed object and performing extension
	 * @param i an ArrayList of TagNode representing the BED data to be extended
	 * @param ext an integer representing the upstream and downstream extension size
	 */
	public ExtendBed(ArrayList<TagNode> i,int ext){
		input = i;
		extSize = ext;
		output = new ArrayList<>();
		set();
	}
	/**
	 * Access the extended data
	 * @return an ArrayList of TagNode representing the extended data
	 */
	public ArrayList<TagNode> getResults(){
		return output;
	}
	/**
	 * Extend the data by the extension size
	 */
	private void set(){
		for (TagNode tagNode : input) {
			String chr = tagNode.getChrom();
			int start = tagNode.getStart();
			if (start - extSize > 0) {
				start = start - extSize;
			}
			int stop = tagNode.getStop() + extSize;
			TagNode temp = new TagNode(chr, start, stop);
			output.add(temp);
		}
	}
	

}
