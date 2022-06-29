/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.distributions;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

import java.io.Serializable;


/**
 * This class implements a multi-variate Gaussian distribution.
 */
public class MultiGaussianDistribution implements Serializable
{
	final private int dimension;
	final private double[] mean;
	final private double[][] covariance;
	private double[][] covarianceInv;
	private double covarianceDet;
	private final double divider;
	private static final Algebra algebra = new Algebra(1.0E-20);
	private boolean inf;


	/**
	 * Creates a new pseudo-random, multivariate gaussian distribution.
	 *
	 * @param mean The mean vector of the generated numbers.  This array is
	 *             copied.
	 * @param covariance The covariance of the generated numbers.  This array
	 *                   is copied.  <code>covariance[r][c]</code> is the
	 *                   element at row <code>r</code> and column
	 *                   <code>c</code>.
	 */
	public MultiGaussianDistribution(double[] mean, double[][] covariance)
	{
		if (covariance.length != covariance[0].length) {
			throw new IllegalArgumentException("Covariance must be a square matrix");
		}

		dimension = covariance.length;
		if (mean.length != dimension) {
			throw new IllegalArgumentException("mean and covariance dimensions don't match");
		}

		this.inf = false;
		this.mean = mean.clone();
		this.covariance = covariance.clone();
		getDetAndInv();
		this.divider = Math.pow(2. * Math.PI, ((double) dimension) / 2.) * Math.pow(covarianceDet, .5);
	}

	private void getDetAndInv() {
		DoubleMatrix2D matrix = new DenseDoubleMatrix2D(covariance);
		covarianceDet = algebra.det(matrix);
		// If covarianceDet is zero, the probability will always be Infinity
		if (covarianceDet == 0) {
			this.inf = true;
			return;
		}
		covarianceInv = algebra.inverse(matrix).toArray();
		// Since covarianceInv will always be a symmetric matrix, we premultiply the top part of the matrix by a factor -1
		// We premultiply the diagonal by a factor -0.5
		// This way, we only need to loop over the diagonal and this top part to calculate the probability
		for (int j = 0; j<dimension-1; j++){
			for (int k = j+1; k<dimension; k++){
				covarianceInv[j][k] *= -1.0;
			}
		}
		for (int k = 0; k<dimension; k++){
			covarianceInv[k][k] *= -0.5;
		}
	}


	@SuppressWarnings("unused")
	public int dimension()
	{
		return dimension;
	}


	/**
	 * Returns (a copy of) this distribution's mean vector.
	 *
	 * @return This distribution's mean vector.
	 */
	public double[] mean()
	{
		return mean.clone();
	}


	/**
	 * Returns (a copy of) this distribution's covariance matrix.
	 *
	 * @return This distribution's covariance matrix.
	 */
	public double[][] covariance()
	{
		return covariance.clone();
	}
	public double probability(double[] v, double[] between)
	{
		if (this.inf){
			return Double.POSITIVE_INFINITY;
		}
		double res = 0;
		for (int j = 0; j<dimension; j++){
			between[j] = v[j] - mean[j];
			res += between[j] * between[j] * covarianceInv[j][j];
		}

		for (int j = 0; j<dimension-1; j++){
			double t1 = 0;
			for (int k = j+1; k<dimension; k++){
				t1 += between[k] * covarianceInv[j][k];
			}
			res += between[j] * t1;
		}

		return Math.exp(res) / divider;
	}

}


