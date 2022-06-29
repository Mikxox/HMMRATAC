package RobustHMM;
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
import java.util.List;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;

public class RobustHMM {
	private final int[] states;
	/**
	 * Constructor
	 * @param o a List of Observations representing data to be decoded
	 * @param h an initial hmm
	 */
	public RobustHMM(List<double[]> o, Hmm<ObservationVector> h){
		ViterbiCalculator vit = new ViterbiCalculator(o, h);
		states = vit.stateSequence();
	}
	/**
	 * Access the state annotations
	 * @return an Array of integers representing the state assignments after decoding
	 */
	public int[] getStates(){
		return states;
	}
}
