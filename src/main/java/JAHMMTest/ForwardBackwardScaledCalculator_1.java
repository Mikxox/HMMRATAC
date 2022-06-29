/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;


/**
 * This class can be used to compute the probability of a given observations
 * sequence for a given HMM.
 * <p>
 * This class implements the scaling method explained in <i>Rabiner</i> and 
 * <i>Juang</i>, thus the {@link #alphaElement(int,int) alphaElement} and
 * {@link #betaElement(int,int) betaElement} return the scaled alpha and
 * beta elements.  The <code>alpha</code> array must always be computed
 * because the scaling factors are computed together with it.
 * <p>
 * For more information on the scaling procedure, read <i>Rabiner</i> and 
 * <i>Juang</i>'s <i>Fundamentals of speech recognition</i> (Prentice Hall,
 * 1993).
 */
public class ForwardBackwardScaledCalculator_1
extends ForwardBackwardCalculator_1
{
	/*
	 * Warning, the semantic of the alpha and beta elements are changed;
	 * in this class, they have their value scaled.
	 */
	// Scaling factors
	private final double[] ctFactors;
	
	//Scaling constant
	private final int constant;
	
	/**
	 * Computes the probability of occurence of an observation sequence
	 * given a Hidden Markov Model.  The algorithms implemented use scaling
	 * to avoid underflows.
	 *
	 * @param hmm A Hidden Markov Model;
	 * @param oseq An observations sequence.
	 * @param flags How the computation should be done. See the
	 *              {@link ForwardBackwardCalculator_1.Computation}.
	 *              The alpha array is always computed.
	 */
	public <O extends Observation> 
	ForwardBackwardScaledCalculator_1(List<? extends O> oseq,
			Hmm<O> hmm, EnumSet<Computation> flags, int c, double[][] alpha, double[][] beta, int threadAmount)
	{
		super(hmm, oseq.size());
		if (oseq.isEmpty()) {
			throw new IllegalArgumentException();
		}

		ctFactors = new double[oseq.size()];
		constant = c;
		int dim = ((ObservationVector) oseq.get(0)).dimension();
		int threadSize = oseq.size() / threadAmount;

		ArrayList<Thread> threads = new ArrayList<>();
		for (int threadCounter = 0; threadCounter < threadAmount; threadCounter++){
			int start = threadSize * threadCounter;
			int end;
			if (threadCounter == threadAmount - 1){
				end = oseq.size(); // Do this to avoid rounding error
			} else {
				end = threadSize * (threadCounter + 1);
			}

			Thread thread = new Thread(() -> {
				double[] matrix1 = new double[dim];
				for (int t = start; t < end; t++) {
					O observation = oseq.get(t);
					computeProbability(hmm, observation, t, matrix1);
				}
			});
			thread.start();
			threads.add(thread);
		}

		for (Thread thread : threads){
			try {
				thread.join();
			} catch (InterruptedException ignored){}
		}

		this.alpha = alpha;
		computeAlpha(hmm);

		if (flags.contains(Computation.BETA)) {
			this.beta = beta;
			computeBeta(hmm);
		}
	}
	
	
	/* Computes the content of the scaled alpha array */
	protected <O extends Observation> void
	computeAlpha(Hmm<? super O> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			computeAlphaInit(hmm, i);
		}
		scale(ctFactors, alpha, 0);

		for (int t = 1; t < size; t++) {
			Arrays.fill(alpha[t], .0);
			for (int j = 0; j < hmm.nbStates(); j++) {
				for (int i = 0; i < hmm.nbStates(); i++) {
					alpha[t][j] += alpha[t - 1][i] * hmm.getAij(i, j);
				}
				alpha[t][j] *= probabilities[j][t];
			}
			scale(ctFactors, alpha, t);
		}
	}
	
	
	/* Computes the content of the scaled beta array.  The scaling factors are
	 those computed for alpha. */
	protected <O extends Observation> void 
	computeBeta(Hmm<? super O> hmm)
	{
		for (int i = 0; i < hmm.nbStates(); i++) {
			beta[size - 1][i] = probabilities[i][size - 1] / ctFactors[size - 1];
		}
		double temp;
		for (int t = size - 2; t >= 0; t--) {
			Arrays.fill(beta[t], .0);

			for (int j = 0; j < hmm.nbStates(); j++) {
				temp = beta[t + 1][j] / ctFactors[t];
				for (int i = 0; i < hmm.nbStates(); i++) {
					beta[t][i] += hmm.getAij(i, j) * temp * probabilities[i][t];
				}
			}
		}
	}
	
	
	/* Normalize alpha[t] and put the normalization factor in ctFactors[t] */
	private void scale(double[] ctFactors, double[][] array, int t)
	{
		double[] table = array[t];
		double sum = 0.;

		for (double v : table) {
			sum += v;
		}
		
		//added: multiply sum by constant. original code did not multiply sum by constant
		//to force original method, set constant to one
		ctFactors[t] = sum * constant;
		for (int i = 0; i < table.length; i++) {
			table[i] /= sum * constant;
		}
	}

	public void clean(){
		super.clean();
	}
}
