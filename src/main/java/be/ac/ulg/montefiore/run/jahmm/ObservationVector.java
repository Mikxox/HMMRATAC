/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm;

import java.text.NumberFormat;


/**
 * This class holds an Observation described by a vector of reals.
 */
public class ObservationVector extends Observation
implements Cloneable
{	
	final double[] value;
	
	
	/**
	 * An observation whose components are 0.
	 *
	 * @param dimension The dimension of the resulting vector.
	 */
	public ObservationVector(int dimension)
	{
		if (dimension <= 0)
			throw new IllegalArgumentException("Dimension must be strictly " +
					"positive");
		
		this.value = new double[dimension];
	}
	
	
	/**
	 * An observation that can be described by a vector of reals.
	 *
	 * @param value The value of this observation.  This array is copied.
	 */
	public ObservationVector(double[] value)
	{
		this(value.length);

		System.arraycopy(value, 0, this.value, 0, value.length);
	}
	
	
	/**
	 * Returns the dimension of this vector.
	 */
	public int dimension()
	{
		return value.length;
	}
	
	
	/**
	 * Returns the values composing this observation.
	 *
	 * @return The values of this observation. The array is copied.
	 */
	@SuppressWarnings("unused")
	public double[] values()
	{
		return value.clone();
	}
	
	
	/**
	 * Returns one of the values composing the observation.
	 *
	 * @param i The dimension of interest (0 &le; i &lt; dimension).
	 * @return The value of the (i+1)-th dimension of this observation.
	 */
	public double value(int i)
	{
		return value[i];
	}

	public String toString(NumberFormat numberFormat)
	{
		StringBuilder s = new StringBuilder("[");

		for (double v : value) s.append(" ").append(numberFormat.format(v));
		
		return s + " ]";
	}
	
	
	public ObservationVector clone()
	{
		return new ObservationVector(value);
	}
}
