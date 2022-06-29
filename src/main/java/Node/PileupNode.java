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
public class PileupNode {
	
	
	private final int _base;
	private final double _score;
	private final String _chrom;
	/**
	 * Constructor for creating new PileupNode
	 * @param base an integer representing the base position
	 * @param score a double representing the base's score
	 * @param chrom a String representing the chromosome name
	 */
	public PileupNode(int base, double score, String chrom){
		_base = base;
		_score = score;
		_chrom = chrom;
	}

	/**
	 * Access the base position
	 * @return an integer representing the base position
	 */
	public int getBase(){
		return _base;
	}
	/**
	 * Access the score
	 * @return a double representing the score at the position
	 */
	public double getScore(){
		return _score;
	}
	/**
	 * Access the chromosome name
	 * @return a String representing the chromosome na,e
	 */
	public String getChrom(){
		return _chrom;
	}
}
