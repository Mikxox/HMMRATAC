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
import java.util.HashMap;

import Node.OverlapNode;
import Node.TagNode;

import static HMMR_ATAC.Utils.toMap;

public class SubtractBed {
	
	private final ArrayList<TagNode> input;
	private final ArrayList<TagNode> exclude;
	private final ArrayList<TagNode> output;
	
	/**
	 * Constructor for creating a SubtractBed object and subtracting the data
	 * @param i an ArrayList of TagNode representing the data to be subtracted from
	 * @param e an ArrayList of TagNode representing the data to be subtracted by
	 */
	public SubtractBed(ArrayList<TagNode> i, ArrayList<TagNode> e){
		input = i;
		exclude = e;
		output = new ArrayList<>();
		subtract();
	}
	/**
	 * Access the subtracted data
	 * @return an ArrayList of TagNode representing the subtracted data
	 */
	public ArrayList<TagNode> getResults(){return output;}
	/**
	 * Subtract the data
	 */
	private void subtract(){
		HashMap<String,ArrayList<TagNode>> in = toMap(input);
		HashMap<String,ArrayList<TagNode>> ex = toMap(exclude);
		for (String chr : in.keySet()){
			if (ex.containsKey(chr)){
				ArrayList<TagNode> inTemp = in.get(chr);
				ArrayList<TagNode> exTemp = ex.get(chr);
				for (TagNode tagNode : inTemp) {
					ArrayList<OverlapNode> results = new ArrayList<>();
					for (TagNode value : exTemp) {
						OverlapNode node = overlap(tagNode, value);
						if (node.hasHit()) {
							// If A is completely consumed by B, dont add as there wont be a subtraction.
							if (!node.isConsumed()) {
								results.add(node);
							}
						}
					}
					//No hits, therefore add the whole A entry into output
					if (results.size() == 0) {
						output.add(tagNode);
					}
					//Only one hit recorded. Only need to find one set of outputs
					else if (results.size() == 1) {
						TagNode hit = results.get(0).getHit();
						/*
						 * +++++++++++
						 *    ----
						 * ===    ====
						 */
						if (hit.getStart() > tagNode.getStart() && hit.getStop() < tagNode.getStop()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStart()));
							output.add(new TagNode(hit.getChrom(), hit.getStop(), tagNode.getStop()));
						}
						/*
						 * ++++++++++
						 * -------
						 *        ===
						 */
						else if (hit.getStart() == tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), hit.getStop(), tagNode.getStop()));
						}
						/*
						 * +++++++++++
						 *    --------
						 * ===
						 */
						else if (hit.getStop() == tagNode.getStop()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStop()));
						}
						/*
						 *     +++++++++
						 * ---------
						 *     =====
						 */
						else if (hit.getStart() < tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), tagNode.getStart(), hit.getStop()));
						}
						/*
						 * ++++++++
						 *     -------
						 *     ====
						 */
						else if (hit.getStart() > tagNode.getStart()) {
							output.add(new TagNode(hit.getChrom(), hit.getStart(), tagNode.getStop()));
						}
					}
					//More than one overlap
					else {
						//Scan the hits to look for which bases in A survive. Then report contiguous intervals that survive
						ArrayList<TagNode> res = new ArrayList<>();
						for (OverlapNode result : results) {
							res.add(result.getHit());
						}
						res.sort(TagNode.basepairComparator);

						int index;
						if (res.get(0).getStart() < tagNode.getStart()) {
							output.add(new TagNode(tagNode.getChrom(), res.get(0).getStop(),
									res.get(1).getStart()));
							index = 1;
						} else {
							output.add(new TagNode(tagNode.getChrom(), tagNode.getStart(),
									res.get(0).getStart()));
							index = 0;
						}
						for (int x = index; x < res.size() - 1; x++) {
							output.add(new TagNode(tagNode.getChrom(), res.get(x).getStop()
									, res.get(x + 1).getStart()));
						}
						if (res.get(res.size() - 1).getStop() < tagNode.getStop()) {
							output.add(new TagNode(tagNode.getChrom(), res.get(res.size() - 1).getStop(),
									tagNode.getStop()));
						}
					}

				}
			
			} else {
				output.addAll(in.get(chr));
			}
		}
		
	}
	/**
	 * Determine if two entries overlap each other
	 * @param node1 a TagNode representing one entry
	 * @param node2 a TagNode representing a seconfd entry
	 * @return a OverlapNode representing the overlap between the two TagNode
	 */
	public static OverlapNode overlap(TagNode node1,TagNode node2){
		if(node1.getChrom().equals(node2.getChrom())){
			int start1 = node1.getStart();
			int start2 = node2.getStart();
			int stop1 = node1.getStop();
			int stop2 = node2.getStop();

			/*
			 * ++++++++++++
			 *    -----
			 *    =====
			 */
			if (start1 <= start2 && stop1 >= stop2){
				return new OverlapNode(node2,true,false);
			}
			
			/*
			 *    ++++++
			 *  -----  
			 *    ===
			 */
			else if(start1 >= start2 && start1 <= stop2){
				return new OverlapNode(new TagNode(node1.getChrom(),node1.getStart(),node2.getStop()),true,false);
			}
			
			/*
			 * ++++++++++
			 *      --------
			 *      =====
			 */
			else if (stop1 >= start2 && stop1 <= stop2){
				return new OverlapNode(new TagNode(node1.getChrom(),node2.getStart(),node1.getStop()),true,false);
			}
			
			/*
			 *    +++++++
			 * --------------   
			 */
			else if (start1 >= start2 && stop1 <= stop2){
				return new OverlapNode(node1,true,true);
			}
			else{
				return new OverlapNode(null,false,false);
			}
		} else{
			return new OverlapNode(null,false,false);
		}
	}

	public static boolean hasOverlap(TagNode node1, TagNode node2){
		if(node1.getChrom().equals(node2.getChrom())){
			int start1 = node1.getStart();
			int start2 = node2.getStart();
			int stop1 = node1.getStop();
			int stop2 = node2.getStop();

			return stop1 >= start2 && stop2 >= start1;
		}
		return false;
	}
}

