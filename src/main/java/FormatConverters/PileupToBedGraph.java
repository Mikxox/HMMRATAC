package FormatConverters;

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

import Node.PileupNode;
import Node.TagNode;

public class PileupToBedGraph {

	
	private ArrayList<TagNode> _bedGraph;
	private final int step;
	
	/**
	 * Constructor for creating new PileupToBedGraph object and generating the data
	 * @param pile an ArrayList of PileupNode2's to convert into a bedgraph
	 * @param x an integer representing the step in the pileup
	 */
	public PileupToBedGraph(ArrayList<PileupNode> pile, int x){
		HashMap<String,ArrayList<PileupNode>> map = toMap(pile);
		step = x;
		run(map);
	}

	/**
	 * WARNING This constructor only works with step size 1!
	 */
	public PileupToBedGraph(double[] pile, TagNode t){
		step = 1;
		run(t, pile);
	}

	public PileupToBedGraph(int[] pile, int step, TagNode t){
		this.step = step;
		run(t, pile);
	}
	/**
	 * Convert the entire pileup into bedgraph
	 */
	private void run(HashMap<String,ArrayList<PileupNode>> map){
		_bedGraph = new ArrayList<>();
		for (String chr : map.keySet()){
			_bedGraph.addAll(generateGraph2(map.get(chr)));
		}
	}

	private void run(TagNode t, double[] array){
		_bedGraph = new ArrayList<>();
		_bedGraph = generateGraph2(t, array);
	}
	private void run(TagNode t, int[] array){
		_bedGraph = new ArrayList<>();
		_bedGraph = generateGraph2(t, array);
	}
	/**
	 * Convert one portion of the pileup into bedgraph
	 * @param states an ArrayList of PileupNode2's to convert
	 * @return an ArrayList of TagNode's representing the bedgraph data
	 */
	private ArrayList<TagNode> generateGraph2(ArrayList<PileupNode> states){
		ArrayList<TagNode> bedGraph = new ArrayList<>();
		
		String chr = states.get(0).getChrom();
		double value = states.get(0).getScore();
		int start = states.get(0).getBase();
		
		for (int i = 0; i < states.size()-1; i++){
			double newValue = states.get(i).getScore();
			int currentBase = states.get(i).getBase();
			int nextBase = states.get(i+1).getBase();

			if (newValue != value || currentBase != (nextBase-step)){
				bedGraph.add(new TagNode(chr, start, currentBase, value));
				start = currentBase;
				value = newValue;
				
			}
		}
		return bedGraph;
	}

	/**
	 * This should only be called when step size is 1
	 */
	private ArrayList<TagNode> generateGraph2(TagNode t, double[] states){
		int offset = t.getStart();
		int start = t.getStart();
		String chr = t.getChrom();
		ArrayList<TagNode> bedGraph = new ArrayList<>();
		double value = states[0];

		for (int i = 0; i < states.length-1; i++){
			double newValue = states[i];
			if (newValue != value){
				bedGraph.add(new TagNode(chr, start, offset+i, value));
				start = offset+i;
				value = newValue;

			}
		}
		return bedGraph;
	}

	private ArrayList<TagNode> generateGraph2(TagNode t, int[] states){
		int offset = t.getStart();
		int start = t.getStart();
		String chr = t.getChrom();
		ArrayList<TagNode> bedGraph = new ArrayList<>();
		int value = states[0];

		for (int i = 0; i < states.length-1; i++){
			int newValue = states[i];
			if (newValue != value){
				bedGraph.add(new TagNode(chr, start, offset+(i*step), value));
				start = offset+(i*step);
				value = newValue;
			}
		}
		int remainder = t.getLength() % 10;
		bedGraph.add(new TagNode(chr, start, offset+((states.length-1)*step)-remainder, value));
		return bedGraph;
	}

	/**
	 * Access the finished bedgraph data
	 * @return an ArrayList of TagNode's representing the completed bedgraph
	 */
	public ArrayList<TagNode> getBedGraph(){
		return _bedGraph;
	}
	/**
	 * Split the pileup into individual chromosomes for ease of computation
	 * @return a HashMap of ArrayList's of PileupNode2's. Each entry in the HasHMap represents one chromosome
	 */
	private HashMap<String,ArrayList<PileupNode>> toMap(ArrayList<PileupNode> _pileup){
		HashMap<String,ArrayList<PileupNode>> map = new HashMap<>();

		for (PileupNode pileupNode : _pileup) {
			String key = pileupNode.getChrom();
			ArrayList<PileupNode> temp;
			if (map.containsKey(key)) {
				temp = map.get(key);
			} else {
				temp = new ArrayList<>();
			}
			temp.add(pileupNode);
			map.put(key, temp);
		}
		
		return map;
	}
}
