/*
 * Copyright (c) 2004-2009, Jean-Marc François. All Rights Reserved.
 * Licensed under the New BSD license.  See the LICENSE file.
 */

package be.ac.ulg.montefiore.run.jahmm;

import java.util.*;
import java.io.*;
import java.text.NumberFormat;


/**
 * Objects implementing this interface represent an observation probability
 * (distribution) function.
 * <p>
 * An <code>Opdf</code> can represent a probability function (if the 
 * observations can take discrete values) or a probability distribution (if
 * the observations are continous).
 */
public interface Opdf<O extends Observation> 
extends Cloneable, Serializable
{
    
    /**
     * Returns the probability (density) of an observation given a distribution.
     *
     * @param o An observation.
     * @return The probability (density, if <code>o</code> takes continuous
     *         values) of <code>o</code> for this function.
     */
    double probability(O o, double[] array);
    double probability(double[] o, double[] array);


    /**
     * Fits this observation probability (distribution) function to a
     * weighted (non empty) set of observations.  Equations (53) and (54)
     * of Rabiner's <i>A Tutorial on Hidden Markov Models and Selected 
     * Applications in Speech Recognition</i> explain how the weights can be
     * used.
     *
     * @param co A set of observations compatible with this factory.
     * @param weights The weight associated to each observation (such that
     *                <code>weight.length == o.length</code> and the sum of
     *                all the elements equals 1).
     * @param state the state to use for weights[loop][state]
     */
    void fit(ArrayList<? extends O> co, double[][] weights, int state);
    
    
    /**
     * Returns a {@link java.lang.String String} describing this distribution.
     * 
     * @param numberFormat A formatter used to convert the numbers (<i>e.g.</i>
     *      probabilities) to strings.
     * @return A {@link java.lang.String String} describing this distribution.
     */
    String toString(NumberFormat numberFormat);
    
    
    Opdf<O> clone();
}
