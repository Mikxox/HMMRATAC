/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm;

import java.text.NumberFormat;
import java.util.ArrayList;

import be.ac.ulg.montefiore.run.distributions.MultiGaussianDistribution;


/**
 * This class represents a multivariate gaussian distribution function.
 */
public class OpdfMultiGaussian
implements Opdf<ObservationVector>
{	
	private MultiGaussianDistribution distribution;
	private final int dimension;
	
	/**
	 * Builds a new gaussian probability distribution with a given mean and
	 * covariance matrix.
	 *
	 * @param mean The distribution's mean.
	 * @param covariance The distribution's covariance matrix.
	 */
	public OpdfMultiGaussian(double[] mean, double[][] covariance)
	{		
		if (covariance.length == 0 || mean.length != covariance.length ||
				covariance.length != covariance[0].length)
			throw new IllegalArgumentException();
		
		distribution = new MultiGaussianDistribution(mean, covariance);
		dimension = covariance.length;
	}
	
	
	/**
	 * Returns (a copy of) this distribution's mean vector.
	 *
	 * @return The mean vector.
	 */
	public double[] mean()
	{
		return distribution.mean();
	}
	
	
	/**
	 * Returns this distribution's covariance matrix.
	 *
	 * @return The covariance matrix.
	 */
	public double[][] covariance()
	{
		return distribution.covariance();
	}

	public double probability(double[] o, double[] array)
	{
		return distribution.probability(o, array);
	}

	public double probability(ObservationVector observationVector, double[] array) {
		return distribution.probability(observationVector.value, array);
	}

	public void fit(ArrayList<? extends ObservationVector> co,
					double[][] weights, int state)
	{
		if (co.isEmpty() || co.size() != weights.length)
			throw new IllegalArgumentException();

		double[] mean = new double[dimension];
		int t = 0;
		// Even though observation vectors never change, the weights do, so we can't precalculate this
		// Compute mean
		for (ObservationVector o : co) {
			for (int r = 0; r < dimension; r++) {
				mean[r] += o.value[r] * weights[t][state];
			}
			t++;
		}

		// Compute covariance
		double[][] covariance = new double[dimension][dimension];
		double[] obs = new double[dimension];
		int i = 0;
		for (ObservationVector o : co) {
			for (int j = 0; j < dimension; j++) {
				obs[j] = o.value[j] - mean[j];
			}

			// using colt matrix zMult to onestep this is not faster due to dimension being so small
			// obsmeaned.zMult(obsmeaned, covariance, weights[i], 1, false, true);
			for (int r = 0; r < dimension; r++) {
				double temp = obs[r] * weights[i][state];
				for (int c = 0; c < dimension; c++) {
					covariance[r][c] += temp * obs[c];
				}
			}
			i++;
		}
		distribution = new MultiGaussianDistribution(mean, covariance);
	}

	public OpdfMultiGaussian clone()
	{
		try {
			return (OpdfMultiGaussian) super.clone();
		} catch(CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
	}
	
	
	public String toString()
	{
		return toString(NumberFormat.getInstance());
	}
	
	
	public String toString(NumberFormat numberFormat)
	{
		StringBuilder s = new StringBuilder("Multi-variate Gaussian distribution --- Mean: [ ");
		double[] mean = distribution.mean();

		for (double v : mean) {
			s.append(numberFormat.format(v)).append(" ");
		}
		
		return s + "]";
	}
}
