package HMMR_ATAC;

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
import java.util.ArrayList;
import java.util.Objects;

import JAHMMTest.BaumWelchScaledLearner;
import JAHMMTest.FitRobust;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;

public class BaumWelch {
	
	private final Hmm<ObservationVector> h;
	private final ArrayList<ArrayList<ObservationVector>> obs;
	private final int maxIter;
	private final double epsilon;
	private final int threadAmount;
	/**
	 * Constructor for creating new BaumWelch object
	 * @param H a hidden markov model
	 * @param o a List of a List of ObservationVector
	 * @param i an integer representing the maximum iterations to perform
	 */
	public BaumWelch(Hmm<ObservationVector> H, ArrayList<ArrayList<ObservationVector>> o, int i, int threadAmount){
		h=H;
		obs=o;
		maxIter=i;
		epsilon = 0.001;
		this.threadAmount = threadAmount;
	}
	
	/**
	 * Build the model
	 * @return a refined model after Baum Welch training
	 */
	public Hmm<ObservationVector> build(){
		Hmm<ObservationVector> firstHmm = h;
		checkModel(firstHmm);
		BaumWelchScaledLearner sbw = new BaumWelchScaledLearner(obs.get(0).size(), firstHmm.nbStates(), threadAmount);
		Hmm<ObservationVector> scaled = null;
		int iter = 0;
		while (iter < maxIter  ){
			try {
				scaled = sbw.iterate(firstHmm, obs);
				checkModel(scaled);

			} catch(IllegalArgumentException e){
				// This happens due to covariance matrix being NaN, can have multiple causes like determinant being so small Java rounds to zero etc.
				System.out.println("Rounding error during training! Model might not be fully converged! You may want to retrain with a different seed or more training data.");
				break;
			}
			if (converged(scaled,firstHmm)){
				break;
			}
			iter += 1;
			firstHmm = scaled;
		}
		sbw.clean();
		//Set proportional initial probabilities
		for (int i = 0; i < Objects.requireNonNull(scaled).nbStates(); i++){
			scaled.setPi(i,(double) 1/scaled.nbStates());
		}
		if (!Double.isNaN(scaled.getAij(0, 0))){
			return scaled;
		} else if (!Double.isNaN(firstHmm.getAij(0, 0))){
			for (int i = 0; i < Objects.requireNonNull(scaled).nbStates(); i++) {
				firstHmm.setPi(i,(double) 1/firstHmm.nbStates());
			}
			return firstHmm;
		} else {
			return h;
		}
	}
	/**
	 * Check the model
	 *
	 * @param hmm HMM to check
	 */
	private void checkModel(Hmm<ObservationVector> hmm){
		for (int i = 0; i < hmm.nbStates(); i++){
			OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
			double[][] cov = pdf.covariance();
			FitRobust fitter = new FitRobust(cov);
			double[][] temp = fitter.getCovariance();
			OpdfMultiGaussian t = new OpdfMultiGaussian(pdf.mean(),temp);
			hmm.setOpdf(i, t);
		}
	}
	/**
	 * Access whether the model converged
	 * @param h1 HMM before BW iteration
	 * @param h2 HMM after BW iteration
	 * @return	a boolean representing whether the model converged within epsilon value
	 */
	private boolean converged(Hmm<ObservationVector> h1, Hmm<ObservationVector> h2){
		int counter = 0;
		for (int i = 0;i < h1.nbStates();i++){
			OpdfMultiGaussian pdf1 = (OpdfMultiGaussian) h1.getOpdf(i);
			OpdfMultiGaussian pdf2 = (OpdfMultiGaussian) h2.getOpdf(i);
			double[] value1 = pdf1.mean();
			double[] value2 = pdf2.mean();
			for (int a = 0; a < value1.length; a++){
				if (Math.abs(value1[a] - value2[a]) > epsilon){
					counter += 1;
				}
			}
		}

		return counter == 0;
	}
}
