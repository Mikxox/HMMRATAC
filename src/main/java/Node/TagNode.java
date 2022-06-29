package Node;

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

import java.util.Comparator;

public class TagNode {
	private String uniqID = "";
	private final String CHROM;
	private final int BP_START;
	private int BP_STOP;
	private double score2 = 0.0;
	private String score3 = "";
	private TagNode summit=null;
	private TagNode upstream=null;
	private TagNode downstream=null;
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 * @param uniqId A String representing a unique ID for the entry
	 */
	public TagNode (String chr, int start, int stop,String uniqId) {
		CHROM = chr;
		BP_START = start;
		BP_STOP = stop;
		uniqID = uniqId;

	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 */
	public TagNode(String chr, int start, int stop) {
		
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
	}
	/**
	 * Constructor
	 * @param chr a String representing the name of the chromosome
	 * @param start an integer representing the start position
	 * @param stop an integer representing the stop position
	 * @param score a double representing the score for the entry
	 */
	public TagNode(String chr,int start,int stop,double score){
		CHROM=chr;
		BP_START = start;
		BP_STOP = stop;
		score2 = score;
	}
	/**
	 * Set a summit as another TagNode
	 * @param t a TagNode representing the summit
	 */
	public void setSummit(TagNode t){
		summit=t;
	}
	/**
	 * Access the summit 
	 * @return a TagNode representing the summit
	 */
	public TagNode getSummit(){return summit;}
	
	/**
	 * Set the third score variable
	 * @param s a double representing a score
	 */
	public void setScore3(String s){score3 = s;}
	/**
	 * Access third score variable
	 * @return a double representing a score
	 */
	public String getScore3(){return score3;}
	/**
	 * Access the score
	 * @return a double representing the entry's score
	 */
	public double getScore2(){
		return score2;
	}
	/**
	 * Access the length of the entry
	 * @return an integer representing the length of the entry
	 */
	public int getLength(){
		return BP_STOP - BP_START;
	}
	/**
	 * Access textual description of the entry
	 * @return a String representing textual description of the entry
	 */
	public String toString(){
		return CHROM+"\t"+BP_START+"\t"+BP_STOP;
	}
	/**
	 * Access textual description of the entry
	 * @return a String representing textual description of the entry
	 */
	public String toString2(){
		return CHROM+"\t"+BP_START+"\t"+BP_STOP+"\t"+"E"+(int)score2;
	}
	/**
	 * Access textual description of the entry, as a scored bedgraph 
	 * @return a String representing the scored bedgraph entry
	 */
	public String toString_ScoredBdg(){
		return CHROM+"\t"+BP_START+"\t"+BP_STOP+"\t"+"E"+(int)score2+"\t"+score3;
	}
	/**
	 * Access textual description of the entry, as a scored summit 
	 * @return a String representing the scored summit entry
	 */
	public String toString_ScoredSummit(){

		return CHROM+"\t"+ this.getSummit().getStart()+"\t"+ this.getSummit().getStop()+"\t"+uniqID+"\t"+score3;
	}
	/**
	 * Access textual description of the entry, as a HMMR gappedPeak 
	 * @return a String representing the scored summit entry
	 */
	public String toString_gappedPeak(){
		
		if (this.upstream == null){
			upstream = this;
		}
		
		if (this.downstream == null){
			downstream=this;
		}
		
		String middleValues = "3"+"\t"+
				"1"+","+getLength()+","+
					"1"+"\t"+"0,"+
				(getStart()-upstream.getStart())+","+
					((downstream.getStop()-upstream.getStart())-1);

		return CHROM+"\t"+upstream.getStart()+"\t"+
				downstream.getStop()+"\t"+uniqID+"\t"+".\t."+"\t"+BP_START+"\t"+
				BP_STOP+"\t"+"255,0,0"+"\t"+middleValues+"\t"+score3+"\t"+"-1\t-1";
	}
	/**
	 * Set the unique ID
	 * @param newid a String representing the entry's unique ID
	 */
	public void setID(String newid) {
		uniqID = newid;		
	}
	/**
	 * Access the chromosome name
	 * @return a String representing the chromosome name
	 */
	public String getChrom() {
		return CHROM;		
	}
	/**
	 * Access the stop position	
	 * @return an integer representing the stop position
	 */
	public int getStop() {
		return BP_STOP;		
	}
	/**
	 * Set the stop position
	 * @param newbp an integer representing the stop position
	 */
	public void setStop(int newbp) {
		BP_STOP = newbp;		
	}
	/**
	 * Access the start position
	 * @return an integer representing the start position
	 */
	public int getStart() {
		return BP_START;		
	}
	/**
	 * Set an upstream tagnode
	 * @param u A TagNode representing the upstream feature
	 */
	public void setUpstream(TagNode u){upstream =u;}
	/**
	 * Set an upstream tagnode
	 * @param d A TagNode representing the upstream feature
	 */
	public void setDownstream(TagNode d){downstream =d;}
	/**
	 * Method for comparing and sorting TagNode
	 */
	public static final Comparator<TagNode> basepairComparator = Comparator.comparingInt(node -> node.BP_START);
}
