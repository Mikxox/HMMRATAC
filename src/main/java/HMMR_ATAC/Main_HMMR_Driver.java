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

/*
 * Written by: Evan Tarbell 
 * evantarb@buffalo.edu
 *
 * Edited by Michael Meuleman
 * https://github.com/Mikxox
 */

import java.io.*;
import java.util.*;

import ATACFragments.FragPileupGen;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussian;
import be.ac.ulg.montefiore.run.jahmm.io.*;
import FormatConverters.PileupToBedGraph;
import GEMM.HMMR_EM;
import GenomeFileReaders.GenomeFileReader;
import GenomeFileReaders.bedFileReader;
import Node.TagNode;
import RobustHMM.KMeansToHMM;
import RobustHMM.RobustHMM;

import WigMath.bedGraphMath;
import WigMath.pileup;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import static HMMR_ATAC.Utils.toMapAndClear;

public class Main_HMMR_Driver {


	@SuppressWarnings("unchecked") //Binary reader will never be type safe
	public static void main(String[] args) throws IOException, FileFormatException {

		/*
		 * Version number. Change as needed
		 */
		String versionNum = "2.0.0";
		ArgParser p = new ArgParser(args, versionNum);

		//Required inputs
		File bamFile = p.getBam();
		File indexFile = p.getIndex();
		String genomeFile = p.getGenome();

		// Optional Inputs
		// comma separated list of initial mean values for frag dist
		String means = p.getMeans();
		// comma separated list of initial standard deviations for frag dist
		String stddevs = p.getStd();
		// whether to perform fragment dist em
		boolean fragEM = p.getEM();
		// minimum mapping quality of reads to keep
		int minMapQ = p.getMinQ();
		// lower bound for fold change range for choosing training sites
		int lower = p.getLower();
		// upper bound for fold change range for choosing training sites
		int upper = p.getUpper();
		// zscored read coverage to exclude from viterbi decoding
		int zscore = p.getZscore();
		// output name
		String output = p.getOutput();
		// whether to print peaks
		boolean peaks = p.getPeaks();
		// whether to print bedgraph
		boolean bg = p.getBedgraph();
		String blacklist = p.getBlacklist();
		int minLength = p.getMinLength();
		String scoreSys = p.getScore();
		boolean BGScore = p.getBGScore();
		int k = p.getK();
		int trim = p.getTrim();
		String trainingRegions = p.getTrainingRegions();
		int vitWindow = p.getWindow();
		File modelFile = p.getModelFile();
		int maxTrain = p.getMaxTrain();
		boolean rmDup = p.getRemoveDuplicates();
		boolean printExclude = p.getPrintExclude();
		boolean printTrain = p.getPrintTrain();
		long randomTrainSeed = p.getRandomTrainSeed();
		double threshold = p.getThreshold();
		boolean printHMMRTracks = p.getPrintHMMRTracks();
		boolean stopAfterModel = p.getStopAfterModel();
		// If we should use HMMBinaryReader or HMMReader for modelFile
		boolean binary = p.getBinary();
		// Code uses int peak
		int peakState = p.getPeakState();
		int coresTraining = p.getCoresTraining();
		int coresViterbi = p.getCoresViterbi();
		// For run time calculation
		long startTime = System.currentTimeMillis();

		//Declare output name
		if (output == null){
			output = "NA";
		}
		PrintStream log = new PrintStream(output +".log");
		
		//Exit program if BAM file or Index file or Genome file not given
		if (bamFile == null || indexFile == null || genomeFile == null ){
			p.printUsage();
			System.exit(1);
		}
		
		/*
		 * Report version into log file
		 * 
		 */
		log.println("Version:"+"\t"+ versionNum);
		
		/*
		 * Report all input arguments into log file
		 */
		log.println("Arguments Used:");
		for (int i = 0; i < args.length-1;i+=2){
			log.println(args[i]+"\t"+args[i+1]);
		}
		
		//Read in genome size stats
		GenomeFileReader gReader = new GenomeFileReader(genomeFile);
		ArrayList<TagNode> genomeStats = gReader.getMap();
		
		//Read in blacklisted if inputted
		ArrayList<TagNode> black = null;
		if (blacklist != null){
			black = new bedFileReader(blacklist).getData();
		}
		
		/*
		 * Set fragment length distribution parameters. 
		 * Use inputted values to set initial values, if provided.
		 * Else use defaults
		 */
		double[] fragMeans = new double[4];
		double[] fragStddevs = new double[4];
		double[] mode = new double[4];
		mode[1] = mode[2] = mode[3] = 2;
		mode[0]=0.5;
		
		if (means != null){
			String[] mu = means.split(",");
			for (int i = 0; i < mu.length;i++){
				fragMeans[i] = Double.parseDouble(mu[i]);
			}
		} else{
			fragMeans[0] = 50.0;
			fragMeans[1] = 200.0; 
			fragMeans[2] = 400.0;
			fragMeans[3] = 600.0;
		}
		if (stddevs != null){
			String[] std = stddevs.split(",");
			for(int i = 0;i < std.length;i++){
				fragStddevs[i] = Double.parseDouble(std[i]);
			}
		} else{
			Arrays.fill(fragStddevs, 20.0);
		}
		
		
		/*
		 * Pull the lengths from the read data. Fragment distribution uses fragments with length > 100 to train the 
		 nucleosome distributions. the short distribution is set to 50 and remains unchanged
		 Only occurs if EM training occurs, else use the default settings
		*/
		
		
		if (fragEM){
			pullLargeLengths puller = new pullLargeLengths(bamFile, indexFile, minMapQ, genomeStats,java.util.Arrays.copyOf(fragMeans, 4), randomTrainSeed);
			double[] lengths = puller.getSampledLengths(10);
			double[] weights = puller.getWeights();

			//Perform EM training
			
			HMMR_EM em = new HMMR_EM(weights,java.util.Arrays.copyOfRange(fragMeans, 1,4),
					java.util.Arrays.copyOfRange(fragStddevs, 1, 4),lengths);
			em.learn();
			double[] tempMeans = em.getMeans();
			double[] tempLam = em.getLamda();

			for (int i = 0;i < tempMeans.length;i++){
				//This will update the parameters IFF they were updated. If they become NaN, leave as default
				if(!Double.isNaN(tempMeans[i]) && !Double.isNaN(tempLam[i])){
					fragMeans[i+1] = tempMeans[i];
					fragStddevs[i+1] = tempLam[i];
				}
			}
		}
		
		log.println("Fragment Expectation Maximum Done");
		for (int i = 0;i < fragMeans.length;i++){
			log.println("Mean\t"+fragMeans[i]+"\tStdDevs\t"+fragStddevs[i]);
		}

		Hmm<ObservationVector> hmm=null;
		pileup pileupData = new pileup(new SplitBed(genomeStats, vitWindow).getResult(), bamFile, indexFile, 0, rmDup);
		bedGraphMath fc = new bedGraphMath(pileupData.getBedGraph());
		
		//calculate the cpm scaling factor for input into FragPileupGen. Use cpmScale=1 for no scaling
		double cpmScale = pileupData.getCPMScale()/1000000;
		log.println("ScalingFactor\t"+cpmScale);
		
		
		double genomeMean = fc.getMean();
		double genomeStd = fc.getSTD();
		
		ArrayList<TagNode> train = new MergeBed(fc.getBetweenRanges(upper, lower)).getResults();
		
		
		
		ArrayList<TagNode> newTrain = new ArrayList<>();
		if (train.size() < maxTrain){
			maxTrain = train.size();
		}
		
		//Shuffle training list before choosing.
		Collections.shuffle(train, new Random(randomTrainSeed));
		for (int i = 0; i < maxTrain; i++){
			newTrain.add(train.get(i));
		}
		
		
		
		train = newTrain;
		train = new ExtendBed(train,5000).getResults();
		ArrayList<TagNode> exclude = new MergeBed(fc.getAboveZscore(zscore)).getResults();
		ArrayList<TagNode> addBack = exclude;
		
		if (blacklist != null){
			exclude.addAll(black);
			exclude = new MergeBed(exclude).getResults();
		}
			
		if(printExclude){
			PrintStream ex = new PrintStream(output +"_excluded.bed");
			for (TagNode tagNode : exclude) {
				ex.println(tagNode.toString() + "\t" + "exclude");
			}
			ex.close();
		}
		
		
		newTrain = new ArrayList<>();
		for (TagNode tagNode : train) {
			int counter = 0;
			for (TagNode node2 : exclude) {
				if (SubtractBed.hasOverlap(tagNode, node2)) {
					counter++;
				}
			}
			if (counter == 0) {
				newTrain.add(tagNode);

			}
		}
		
		train = newTrain;

		// Allows user to use training set for model generation
		if (trainingRegions != null){
			bedFileReader trainReader = new bedFileReader(trainingRegions);
			train = trainReader.getData();
		}
		
		log.println("Training Regions found and Zscore regions for exclusion found");
		
		
		/*
		 * Create the fragment pileup tracks using the training set and the fragment distribution parameters
		 */
		if (modelFile == null){
		
			
			if (printTrain) {
				PrintStream tr = new PrintStream(output + "_training.bed");
				for (TagNode tagNode : train) {
					tr.println(tagNode.toString() + "\t" + "training");
				}
				tr.close();
			}
			FragPileupGen gen = new FragPileupGen(train, mode, fragMeans, fragStddevs, minMapQ, rmDup,cpmScale);
			TrackHolder holder = new TrackHolder(gen.buildTracks(bamFile, indexFile), trim);
			
			log.println("Training Fragment Pileup completed");
			
			/*
			 * Create the initial model using KMeans and then refine it using Baum-Welch
			 */
			
			KMeansToHMM kmeans = new KMeansToHMM(holder.getDataSet(), k,Integer.MAX_VALUE,true,true,true);
			log.println("Kmeans Model:\n"+kmeans.getHMM().toString());
			
			hmm = new BaumWelch((Hmm<ObservationVector>) kmeans.getHMM(), holder.getBWObs(), 1000, coresTraining).build();
		}

		/*
		 * Use input model if available
		 */
		
		if (modelFile != null){
			if (binary){
				hmm = (Hmm<ObservationVector>) HmmBinaryReader.read(new FileInputStream(modelFile));
			} else {
				BufferedReader inReader = new BufferedReader(new FileReader(modelFile));
				hmm = HmmReader.read(inReader, new OpdfMultiGaussianReader());
				inReader.close();
			}
		}

		
		/*
		 * Identify peak state as the state with the highest total read signals.
		 */
		int peak = peakState;
		if (peak < 0) {
			double max = 0.0;
			for (int i = 0; i < hmm.nbStates(); i++) {
				OpdfMultiGaussian pdf = (OpdfMultiGaussian) hmm.getOpdf(i);
				double sh = Arrays.stream(pdf.mean()).sum();
				if (sh > max) {
					peak = i;
					max = sh;
				}
			}
		}
		
		
		
		/*
		 * Output binary model file
		 */
		if (modelFile == null) {
			File outputModel = new File(output + ".model");
			FileOutputStream outModel = new FileOutputStream(outputModel);
			HmmBinaryWriter.write(outModel, hmm);
			outModel.close();
			log.println("Model created and refined. See " + output + ".model for binary model and " + output + "_readable.model for human readable model.");
			log.println("Model:\n" + hmm);
		}

		/*
		 * Output human-readable model file
		 */
		if (modelFile == null) {
			File outputModel = new File(output + "_readable.model");
			FileOutputStream outModel = new FileOutputStream(outputModel);
			Writer outWriter = new BufferedWriter(new OutputStreamWriter(outModel));
			HmmWriter.write(outWriter, new OpdfMultiGaussianWriter(), hmm);
			outWriter.flush();
			outWriter.close();
			outModel.close();
		}

		long endTime2 = System.currentTimeMillis();
		long total2 = (endTime2 - startTime) / 1000;
		log.println("Model written (seconds): \t"+total2);

		/*
		 * Stop program if only model is desired
		 */
		if(stopAfterModel){
			System.exit(0);
		}
		
		/*
		 * Split the genome file into smaller 25MB chunks 
		 * Can also split into whatever sized chunks the users prefers
		 * May be necessary to split into smaller chunks for machines with less memory
		 */
		
		ArrayList<TagNode> split = new SplitBed(genomeStats, vitWindow).getResult();
		
		/*
		 * Subtract excluded regions from the split genome for Viterbi
		 */
		
		ArrayList<TagNode> vitBed = new SubtractBed(split, exclude).getResults();
		
		log.println("Genome split and subtracted masked regions");
		
		/*
		 * Run viterbi on the whole genome
		 */
		PrintStream NFR = null;
		PrintStream MONO = null;
		PrintStream DI = null;
		PrintStream TRI = null;
		if(printHMMRTracks){
			 NFR = new PrintStream(output +"_nfr.bedgraph");
			 MONO = new PrintStream(output +"_mono.bedgraph");
			 DI = new PrintStream(output +"_di.bedgraph");
			 TRI = new PrintStream(output +"_tri.bedgraph");
		}
		ArrayList<Thread> threads = new ArrayList<>();
		ArrayList<TagNode> genomeAnnotation = new ArrayList<>();
		Hmm<ObservationVector> finalHmm = hmm;
		SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamFile).index(indexFile));
		for (int i = 0; i < vitBed.size(); i++){
			if (vitBed.get(i).getLength() >= 10){
				ArrayList<TagNode> tempBed = new ArrayList<>();
				tempBed.add(vitBed.get(i));
				// Creating these tracks eats up a lot of memory so multithreading this is only feasible on very high-end machines and the time gained is very minor
				FragPileupGen vGen = new FragPileupGen(tempBed, mode, fragMeans, fragStddevs, minMapQ, rmDup, cpmScale);
				ArrayList<double[]> vGenTracks = vGen.buildTracks(reader);
				// If size is zero, this means the current chromosome in the genome file was not present in the input bam file
				// This of course means that no further processing is needed for this part
				if (vGenTracks.size() == 0){
					continue;
				}

				// Make sure we don't create too many threads to avoid memory overflow
				if (threads.size() >= coresViterbi){
					Thread wait = threads.remove(0);
					try {
						wait.join();
					} catch (Exception ignored){}
				}

				TrackHolder vHolder = new TrackHolder(vGenTracks, trim);

				if (printHMMRTracks){
					HMMRTracksToBedgraph tracks = new HMMRTracksToBedgraph(vHolder.getRawData(),vitBed.get(i),10);
					ArrayList<TagNode> nfr = tracks.getShort();
					ArrayList<TagNode> mono = tracks.getMono();
					ArrayList<TagNode> di = tracks.getDi();
					ArrayList<TagNode> tri = tracks.getTri();

					if (nfr != null) {
						for (TagNode tagNode : nfr) {
							NFR.println(tagNode.toString2());
						}
					}
					if (mono != null) {
						for (TagNode tagNode : mono) {
							MONO.println(tagNode.toString2());
						}
					}
					if (di != null) {
						for (TagNode tagNode : di) {
							DI.println(tagNode.toString2());
						}
					}
					if (tri != null) {
						for (TagNode tagNode : tri) {
							TRI.println(tagNode.toString2());
						}
					}

				}

				// Multithread the Viterbi calculations called by RobustHMM HMM
				// The eventual output in genomeAnnotation needs to be sorted before final processing anyway
				// so results being added out of order does not matter.
				int finalI = i;
				TagNode nextVit = vitBed.get(finalI);
				Thread thread = new Thread(() -> {
					RobustHMM HMM = new RobustHMM(vHolder.getTracks(), finalHmm);
					ArrayList<TagNode> temp = new PileupToBedGraph(HMM.getStates(), 10, nextVit).getBedGraph();
					synchronized (genomeAnnotation){
						genomeAnnotation.addAll(temp);
					}

					if (finalI % 50 == 0 || finalI == vitBed.size() - 1) {
						log.println(finalI + " round viterbi done");
					}
				});
				thread.start();
				threads.add(thread);
			}
		}
		try {
			reader.close();
		} catch (IOException ignored) {}

		for (Thread thread : threads){
			try {
				thread.join();
			} catch (Exception ignored){}
		}

		if(printHMMRTracks){
			NFR.close();MONO.close();DI.close();TRI.close();
		}

		/*
		 * Report the final results as peaks, bedgraphs and summits, if desired
		 */
		PrintStream bedgraph=null;
		if (bg){
			 bedgraph = new PrintStream(output +".bedgraph");
		}
		PrintStream pks=null;
		PrintStream summits=null;
		if (peaks){
			 pks = new PrintStream(output +"_peaks.gappedPeak");
			 summits = new PrintStream(output +"_summits.bed");
		}
		HashMap<String,ArrayList<TagNode>> bdg = fc.getMappedBedgraph();
		HashMap<String,ArrayList<TagNode>> hmmrBdg = toMapAndClear(genomeAnnotation);
		double[] filter = bedGraphMath.createFilter(20);
		int counter=1;
		for (String chr : hmmrBdg.keySet()){
			ArrayList<TagNode> hmmr = hmmrBdg.get(chr);
			ArrayList<TagNode> signal = bdg.get(chr);
			List<TagNode> overlaps = new LinkedList<>();
			double[] scores = new double[3];
			if (signal != null) {
				hmmr.sort(TagNode.basepairComparator);
				
				signal.sort(TagNode.basepairComparator);
				int index = 0;
				for (int i = 0; i < hmmr.size(); i++) {
					TagNode temp = hmmr.get(i);

					/*
					 * Execute the scoring commands if the state is a peak or if bedgraph scoring is on
					 */
					if ((int) temp.getScore2() == peak || BGScore) {
						boolean hasHadOverlap = false;
						overlaps.clear();
						int a;
						for (a = index; a < signal.size(); a++) {
							if (SubtractBed.hasOverlap(temp, signal.get(a))) {
								overlaps.add(signal.get(a));
								hasHadOverlap = true;
							} else if (hasHadOverlap) {
								index = a;
								break;
							} else if (signal.get(a).getStart() > temp.getStop()){
								// signal is sorted by start values, so once signal is passed temp, we can never have an overlap
								break;
							}
						}
						double MAX, MEAN, MEDIAN, ZSCORE, FOLDCHANGE;
						if (!hasHadOverlap){
							MAX = 0;
							MEAN = 0;
							MEDIAN = 0;
							ZSCORE = (0 - genomeMean) / genomeStd;
							FOLDCHANGE = 0 / genomeMean;
						} else {
							bedGraphMath.set(temp, overlaps, scores);
							MAX = scores[2];
							MEAN = scores[0];
							MEDIAN = scores[1];
							ZSCORE = (scores[0] - genomeMean) / genomeStd;
							FOLDCHANGE = scores[0] / genomeMean;
						}

						switch (scoreSys) {
							case "ave" -> temp.setScore3(Double.toString(MEAN));
							case "fc" -> temp.setScore3(Double.toString(FOLDCHANGE));
							case "zscore" -> temp.setScore3(Double.toString(ZSCORE));
							case "med" -> temp.setScore3(Double.toString(MEDIAN));
							case "all" -> {
								String ANSWER = MAX + "_" + MEAN + "_" + MEDIAN + "_" + ZSCORE + "_" + FOLDCHANGE;
								temp.setScore3(ANSWER);
							}
							default -> temp.setScore3(Double.toString(MAX));
						}
						if ((int) temp.getScore2() == peak) {
							bedGraphMath.setSmooth(20, temp, overlaps, filter);
							temp.setID("Peak_" + counter);
							if (i > 0) {
								temp.setUpstream(hmmr.get(i - 1));
							} else {
								temp.setUpstream(hmmr.get(i));
							}
							if (i < hmmr.size() - 1) {
								temp.setDownstream(hmmr.get(i + 1));
							} else {
								temp.setDownstream(hmmr.get(i));
							}
							counter++;
						}

					}
					/*
					 * report the bedgraph, is desired
					 */
					if (bg) {
						if (!BGScore) {
							bedgraph.println(temp.toString2());
						} else {
							bedgraph.println(temp.toString_ScoredBdg());
						}
					}
					/*
					 * report the peaks and summits, if desired
					 */
					if (peaks && (int) temp.getScore2() == peak
							&& temp.getLength() >= minLength &&
							Double.parseDouble(temp.getScore3()) >= threshold) {
						if (temp.getSummit() != null) {
							summits.println(temp.toString_ScoredSummit());
						}

						pks.println(temp.toString_gappedPeak());
					}

				}
			}
			
		}
		if (bg){
			bedgraph.close();
		}
		
		
		if (peaks){
			counter=0;

			for (TagNode tagNode : addBack) {
				String chrom = tagNode.getChrom();
				int start = tagNode.getStart();
				int stop = tagNode.getStop();

				pks.println(chrom + "\t" + start + "\t" + stop + "\t" + "HighCoveragePeak_" + counter + "\t.\t.\t0\t0\t255,0,0\t1\t" +
						tagNode.getLength() + "\t0\t-1\t-1\t-1");
			}
			
			pks.close();
			summits.close();
		}

		long endTime = System.currentTimeMillis();
		long total = (endTime - startTime) / 1000;
		log.println("Total time (seconds): \t"+total);
		log.close();
	}
	
	
}

