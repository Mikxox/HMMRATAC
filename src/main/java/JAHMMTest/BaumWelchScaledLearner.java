/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the Baum-Welch learning algorithm.  It uses a
 * scaling mechanism so as to avoid underflows.
 * <p>
 * For more information on the scaling procedure, read <i>Rabiner</i> and 
 * <i>Juang</i>'s <i>Fundamentals of speech recognition</i> (Prentice Hall,
 * 1993).
 */
public class BaumWelchScaledLearner
extends BaumWelchLearner
{
	/**
	 * Initializes a Baum-Welch Scaler Learner instance.
	 * @param i length of first observation
	 * @param j number of nb states of hmm model
	 */
	public BaumWelchScaledLearner(int i, int j, int threadAmount)
	{
		super(i, j, threadAmount);
	}

	public void clean(){
		super.clean();
	}
	
	protected <O extends Observation> ForwardBackwardCalculator_1
	generateForwardBackwardCalculator(List<? extends O> sequence,
			Hmm<O> hmm, double[][] alpha, double[][] beta, int threadAmount)
	{
		/*
		 * Initializes a Baum-Welch algorithm implementation.
		 */
		int scaledConstant = 1;
		return new ForwardBackwardScaledCalculator_1(sequence, hmm,
				EnumSet.allOf(ForwardBackwardCalculator_1.Computation.class), scaledConstant, alpha, beta, threadAmount);
	}
}