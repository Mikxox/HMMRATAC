/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm;

import java.util.Iterator;
import java.util.List;


/**
 * This class can be used to compute the most probable state sequence matching
 * a given observation sequence (given an HMM).
 */
public class ViterbiCalculator
{	
	/*
	 * The psy and delta values, as described in Rabiner and Juand classical
	 * papers.
	 */
	private double[][] delta;
	private int[][] psy;
	private final int[] stateSequence;
	private final double[] array;


	/**
	 * Computes the most likely state sequence matching an observation
	 * sequence given an HMM.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 */
	public <O extends Observation>
	ViterbiCalculator(List<double[]> oseq, Hmm<O> hmm)
	{
		if (oseq.isEmpty()) {
			throw new IllegalArgumentException("Invalid empty sequence");
		}

		array = new double[oseq.get(0).length];

		delta = new double[2][hmm.nbStates()];
		psy = new int[oseq.size()][hmm.nbStates()];
		stateSequence = new int[oseq.size()];

		for (int i = 0; i < hmm.nbStates(); i++) {
			delta[0][i] = -Math.log(hmm.getPi(i)) - Math.log(hmm.getOpdf(i).probability(oseq.get(0), array));
			psy[0][i] = 0;
		}

		Iterator<double[]> oseqIterator = oseq.iterator();
		if (oseqIterator.hasNext()) {
			oseqIterator.next();
		}

		int t = 1;
		while (oseqIterator.hasNext()) {
			double[] observation = oseqIterator.next();

			for (int i = 0; i < hmm.nbStates(); i++) {
				computeStep(hmm, observation, t, i);
			}
			System.arraycopy(delta[1], 0, delta[0], 0, hmm.nbStates());

			t++;
		}

		double lnProbability = Double.MAX_VALUE;
		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisProbability = delta[0][i];

			if (lnProbability > thisProbability) {
				lnProbability = thisProbability;
				stateSequence[oseq.size() - 1] = i;
			}
		}

		for (int t2 = oseq.size() - 2; t2 >= 0; t2--) {
			stateSequence[t2] = psy[t2 + 1][stateSequence[t2 + 1]];
		}

		psy = null;
		delta = null;
	}

	/**
	 * Computes delta and psy[t][j] (t > 0)
	 */
	private <O extends Observation> void
	computeStep(Hmm<O> hmm, double[] o, int t, int j)
	{
		double minDelta = Double.MAX_VALUE;
		int min_psy = 0;

		for (int i = 0; i < hmm.nbStates(); i++) {
			double thisDelta = delta[0][i] - Math.log(hmm.getAij(i, j));

			if (minDelta > thisDelta) {
				minDelta = thisDelta;
				min_psy = i;
			}
		}

		delta[1][j] = minDelta - Math.log(hmm.getOpdf(j).probability(o, array));
		psy[t][j] = min_psy;
	}


	/**
	 * Returns the array containing the computed most likely
	 * state sequence.
	 *
	 * @return The state sequence; the i-th value of the array is the index
	 *         of the i-th state of the state sequence.
	 */
	public int[] stateSequence() 
	{
		return stateSequence;
	}
}
