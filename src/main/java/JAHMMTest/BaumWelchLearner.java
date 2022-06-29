package JAHMMTest;

import java.util.*;

import be.ac.ulg.montefiore.run.jahmm.*;


/**
 * An implementation of the Baum-Welch learning algorithm.  This algorithm
 * finds a HMM that models a set of observation sequences.
 */
public class BaumWelchLearner
{
        private double[][] gamma;
        private double[][] alpha;
        private double[][] beta;
        private int currentI;
        private int currentJ;
        private final int threadAmount;
        /**
         * Initializes a Baum-Welch instance.
         * @param i length of first observation
         * @param j number of nb states of hmm model
         */
        public BaumWelchLearner(int i, int j, int threadAmount){
                gamma = new double[i][j];
                alpha = new double[i][j];
                beta = new double[i][j];
                currentI = i;
                currentJ = j;
                this.threadAmount = threadAmount;
        }

        public void clean(){
//                xi = null;
                gamma = null;
                alpha = null;
                beta = null;
                currentJ = 0;
                currentI = 0;
        }
        
        /**
         * Performs one iteration of the Baum-Welch algorithm.
         * In one iteration, a new HMM is computed using a previously estimated
         * HMM.
         *
         * @param hmm A previously estimated HMM.
         * @param sequences The observation sequences on which the learning is
         *         based.  Each sequence must have a length higher or equal to
         *         2.
         * @return A new, updated HMM.
         */
        public Hmm<ObservationVector>
        iterate(Hmm<ObservationVector> hmm, ArrayList<ArrayList<ObservationVector>> sequences) {
                Hmm<ObservationVector> nhmm;
                nhmm = hmm.clone();

                /* gamma and xi arrays are those defined by Rabiner and Juang */
                /* gamma = gamma array associated to observation sequence n */

                /* a[i][j] = aijNum[i][j] / aijDen[i]
                 * aijDen[i] = expected number of transitions from state i
                 * aijNum[i][j] = expected number of transitions from state i to j
                 */
                double[][] aijNum = new double[hmm.nbStates()][hmm.nbStates()];
                double[] aijDen = new double[hmm.nbStates()];

                Arrays.fill(aijDen, 0.);
                for (int i = 0; i < hmm.nbStates(); i++) {
                        Arrays.fill(aijNum[i], 0.);
                }

                for (ArrayList<ObservationVector> obsSeq : sequences) {
                        ForwardBackwardCalculator_1 fbc = generateForwardBackwardCalculator(obsSeq, hmm, alpha, beta, threadAmount);


                        double[][] xiSum = new double[hmm.nbStates()][hmm.nbStates()];
                        double[] gammaSum = new double[hmm.nbStates()];
                        estimateXi(obsSeq, fbc, hmm, xiSum, gammaSum);

                        aijDen = gammaSum;
                        aijNum = xiSum;

                        // clear all information stored inside fbc to save memory
                        fbc.clean();
                }

                for (int i = 0; i < hmm.nbStates(); i++) {
                        if (aijDen[i] == 0.) { // State i is not reachable
                                for (int j = 0; j < hmm.nbStates(); j++) {
                                        nhmm.setAij(i, j, hmm.getAij(i, j));
                                }
                        }
                        else {
                                for (int j = 0; j < hmm.nbStates(); j++) {
                                        nhmm.setAij(i, j, aijNum[i][j] / aijDen[i]);
                                }
                        }
                }

                /* pdfs computation */
                ArrayList<Thread> threads = new ArrayList<>();
                for (int i = 0; i < hmm.nbStates(); i++) {
                        int finalI = i;
                        Thread thread = new Thread(() -> {
                                double sum = 0.;
                                int j;
                                for (j = 0; j < sequences.get(0).size(); j++) {
                                        sum += gamma[j][finalI];
                                }

                                if (sum != 0.) {
                                        for (j--; j >= 0; j--) {
                                                gamma[j][finalI] /= sum;
                                        }
                                }
                                Opdf<ObservationVector> opdf = nhmm.getOpdf(finalI);

                                opdf.fit(sequences.get(0), gamma, finalI);
                        });
                        thread.start();
                        threads.add(thread);
                }

                for (Thread thread : threads){
                        try {
                                thread.join();
                        } catch (InterruptedException ignored){}
                }

                return nhmm;
        }
        
        protected <O extends Observation> ForwardBackwardCalculator_1
        generateForwardBackwardCalculator(List<? extends O> sequence, Hmm<O> hmm, double[][] alpha, double[][] beta, int threadAmount)
        {
                return new ForwardBackwardCalculator_1(hmm, EnumSet.allOf(ForwardBackwardCalculator_1.Computation.class));
        }
        
        
        protected <O extends Observation> void
        estimateXi(List<? extends O> sequence, ForwardBackwardCalculator_1 fbc,
                        Hmm<O> hmm, double[][] xiSum, double[] gammaSum)
        {       
                if (sequence.size() <= 1) {
                        throw new IllegalArgumentException("Observation sequence too short");
                }

                if (currentI != sequence.size()-1 || currentJ != hmm.nbStates()){
                        currentI = sequence.size()-1;
                        currentJ = hmm.nbStates();
                        gamma = new double[sequence.size()][hmm.nbStates()];
                }

                for (int t = 0; t < sequence.size() - 1; t++) {
                        Arrays.fill(gamma[t], .0);
                        for (int i = 0; i < hmm.nbStates(); i++) {

                                for (int j = 0; j < hmm.nbStates(); j++) {
                                        double temp = fbc.alphaElement(t, j) * hmm.getAij(j, i) * fbc.betaElement(t + 1, i);
                                        gamma[t][j] += temp;
                                        xiSum[j][i] += temp;
                                        gammaSum[j] += temp;
                                }
                        }
                }

                Arrays.fill(gamma[sequence.size()-1], .0);
                for (int j = 0; j < hmm.nbStates(); j++) {
                        for (int i = 0; i < hmm.nbStates(); i++) {
                                double temp = fbc.alphaElement(sequence.size() - 1, j) * hmm.getAij(j, i) * fbc.betaElement(sequence.size()-1, i);
                                gamma[sequence.size()-1][j] += temp;
                                gammaSum[j] += temp;
                        }
                }

        }

}
