package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;


/**
 * This class can be used to compute the probability of a given observations
 * sequence for a given HMM.  Once the probability has been computed, the
 * object holds various information such as the <i>alpha</i> (and possibly
 * <i>beta</i>) array, as described in <i>Rabiner</i> and <i>Juang</i>.
 * <p>
 * Computing the <i>beta</i> array requires a O(1) access time to the 
 * observation sequence to get a theoretically optimal performance.
 */
public class ForwardBackwardCalculator_1
{       
        /**
         * Flags used to explain how the observation sequence probability
         * should be computed (either forward, using the alpha array, or backward,
         * using the beta array).
         */
        public enum Computation { ALPHA, BETA }
        
        
        /* alpha[t][i] = P(O(1), O(2),..., O(t+1), i(t+1) = i+1 | hmm), that is the
         probability of the beginning of the state sequence (up to time t+1)
         with the (t+1)th state being i+1. */
        protected double[][] alpha = null;
        protected double[][] beta = null;
        protected int size;
        // Keep track of all probabilities, Array here is much faster than HashMap in distribution
        protected double[][] probabilities;


        protected ForwardBackwardCalculator_1(Hmm<?> hmm, int osecSize)
        {
                this.size = osecSize;
                probabilities = new double[hmm.nbStates()][osecSize];
        }

        public <O extends Observation> void computeProbability(Hmm<? super O> hmm, O o, int t, double [] array){
                for (int i = 0; i < hmm.nbStates(); i++) {
                        probabilities[i][t] = hmm.getOpdf(i).probability(o, array);
                }
        }

        public void clean(){
                probabilities = null;
                alpha = null;
                beta = null;
        }
        
        
        /**
         * Computes the probability of occurence of an observation sequence
         * given a Hidden Markov Model.
         *
         * @param hmm A Hidden Markov Model;
         * @param flags How the computation should be done. See the
         *              {@link Computation Computation} enum.
         */
        public <O extends Observation>
        ForwardBackwardCalculator_1(Hmm<O> hmm, EnumSet<Computation> flags)
        {
            if (flags.contains(Computation.ALPHA)) {
                    computeAlpha(hmm);
            }
                
            if (flags.contains(Computation.BETA)) {
                    computeBeta(hmm);
            }
        }
        
        /* Computes the content of the alpha array */
        protected <O extends Observation> void
        computeAlpha(Hmm<? super O> hmm)
        {
                alpha = new double[size][hmm.nbStates()];
                
                for (int i = 0; i < hmm.nbStates(); i++) {
                        computeAlphaInit(hmm, i);
                }
                
                for (int t = 1; t < size; t++) {
                        for (int i = 0; i < hmm.nbStates(); i++) {
                                computeAlphaStep(hmm, t, i);
                        }
                }
        }
        
        
        /* Computes alpha[0][i] */
        protected <O extends Observation> void
        computeAlphaInit(Hmm<? super O> hmm, int i)
        {
                alpha[0][i] = hmm.getPi(i) * probabilities[i][0];
        }
        
        
        /* Computes alpha[t][j] (t > 0) */
        protected <O extends Observation> void 
        computeAlphaStep(Hmm<? super O> hmm, int t, int j)
        {
                double sum = 0.;
                for (int i = 0; i < hmm.nbStates(); i++) {
                        sum += alpha[t - 1][i] * hmm.getAij(i, j);
                }
                alpha[t][j] = sum * probabilities[j][t];
        }
        
        
        /* Computes the content of the beta array.  Needs an O(1) access time
         to the elements of oseq to get a theoretically optimal algorithm. */
        protected <O extends Observation> void 
        computeBeta(Hmm<? super O> hmm)
        {
                beta = new double[size][hmm.nbStates()];
                
                for (int i = 0; i < hmm.nbStates(); i++) {
                        beta[size - 1][i] = 1.;
                }
                
                for (int t = size-2; t >= 0; t--) {
                        for (int i = 0; i < hmm.nbStates(); i++) {
                                computeBetaStep(hmm, t, i);
                        }
                }
        }
        
        
        /* Computes beta[t][i] (t < obs. seq.le length - 1) */
        protected <O extends Observation> void 
        computeBetaStep(Hmm<? super O> hmm, int t, int i)
        {
                beta[t][i] = 0;
                for (int j = 0; j < hmm.nbStates(); j++) {
                        beta[t][i] += beta[t + 1][j] * hmm.getAij(i, j) * probabilities[j][t+1];
                }
        }
        
        
        /**
         * Returns an element of the <i>alpha</i> array.
         * 
         * @param t The temporal argument of the array (positive but strictly
         *          smaller than the length of the sequence that helped generating
         *          the array).
         * @param i A state index of the HMM that helped generating the array.
         * @throws UnsupportedOperationException if alpha array has not been computed.
         * @return The <i>alpha</i> array (t, i) element.
         */ 
        public double alphaElement(int t, int i)
        {
                if (alpha == null) {
                        throw new UnsupportedOperationException("Alpha array has not been computed");
                }
                
                return alpha[t][i];
        }
        
        
        /**
         * Returns an element of the <i>beta</i> array.
         * 
         * @param t The temporal argument of the array (positive but smaller than
         *          the length of the sequence that helped generating the array).
         * @param i A state index of the HMM that helped generating the array.
         * @throws UnsupportedOperationException if beta array has not been
         *          computed.
         * @return The <i>beta</i> beta (t, i) element.
         */ 
        public double betaElement(int t, int i)
        {
                if (beta == null) {
                        throw new UnsupportedOperationException("Beta array has not been computed");
                }
                
                return beta[t][i];
        }
}