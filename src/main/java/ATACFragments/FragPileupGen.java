package ATACFragments;

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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import Node.TagNode;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class FragPileupGen {
    private final ArrayList<TagNode> genome;
    private final int minMapQ;
    final AbstractRealDistribution shortDist;
    final AbstractRealDistribution monoDist;
    final AbstractRealDistribution diDist;
    final AbstractRealDistribution triDist;
    private final boolean rmDup;
    private final double scale;

    /**
     * Constructor for creating a new FragPileupGen object
     * @param g an ArrayList of TagNodes representing the genomic regions over which the data will be created
     * @param mode an array of doubles representing the distributions to use for building the data tracks
     * @param means an array of doubles representing the means of the distributions
     * @param lamda an array of doubles representing the standard deviations of the distributions
     * @param q an integer representing the minimum mapping quality score of the reads to use for data track generation
     */
    public FragPileupGen(ArrayList<TagNode> g,double[] mode, double[] means,double[] lamda,
                         int q,boolean r,double s) {
        genome = g;
        shortDist = getDist(mode[0],means[0],lamda[0]);
        monoDist = getDist(mode[1],means[1],lamda[1]);
        diDist = getDist(mode[2],means[2],lamda[2]);
        triDist = getDist(mode[3],means[3],lamda[3]);
        minMapQ=q;
        rmDup = r;
        scale=s;
    }

    /**
     * Builds the data tracks.
     * Then averages the data in the tracks across a 10bp window
     * Then perform a scaling of the tracks based on inputted scaling factor.
     * Then perform a square root transformation of the tracks.
     * @param input a BAM file containing the ATAC-seq paired end reads.
     * @param index a BAM index file containing the index for the BAM file.
     * @return a new ArrayList containing the data
     */
    public ArrayList<double[]> buildTracks(File input, File index) {
        final SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(input).index(index));
        ArrayList<double[]> res = buildTracksReader(reader);
        try {
            reader.close();
        } catch (IOException ignored) {}
        return res;
    }

    /**
     * Builds the data tracks.
     * Then averages the data in the tracks across a 10bp window
     * Then perform a scaling of the tracks based on inputted scaling factor.
     * Then perform a square root transformation of the tracks.
     * @param reader a SamReader made from a bam and corresponding index file
     * @return a new ArrayList containing the data
     */
    public ArrayList<double[]> buildTracks(SamReader reader) {
        return buildTracksReader(reader);
    }

    public ArrayList<double[]> buildTracksReader(SamReader reader) {
        ArrayList<double[]> tracks = new ArrayList<>();

        for (TagNode tagNode : genome) {
            String chr = tagNode.getChrom();
            int bedStart = tagNode.getStart();
            int bedStop = tagNode.getStop();
            int offset = tracks.size();

            CloseableIterator<SAMRecord> iter = reader.query(chr, bedStart, bedStop, false);

            if (!iter.hasNext()){
                iter.close();
                continue;
            }

            tracks.ensureCapacity(offset + bedStop - bedStart + 10);
            for (int z = bedStart; z < bedStop; z++) {
                tracks.add(null);
            }

            while (iter.hasNext()) {
                SAMRecord record = iter.next();
                if (!record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && record.getFirstOfPairFlag()
                        && record.getMappingQuality() >= minMapQ && !(record.getDuplicateReadFlag() && rmDup && record.getInferredInsertSize() != 0)) {

                    int start = 0;
                    int stop = 0;
                    if (record.getInferredInsertSize() > 0) {
                        start = record.getAlignmentStart();
                        stop = record.getAlignmentStart() + record.getInferredInsertSize() - 1;
                    } else if (record.getInferredInsertSize() < 0) {
                        start = record.getAlignmentEnd() + record.getInferredInsertSize() + 1;
                        stop = record.getAlignmentEnd();

                    }

                    int length = stop - start;
                    double sh = shortDist.density(length);
                    double mono = monoDist.density(length);
                    double di = diDist.density(length);
                    double tri = triDist.density(length);

                    double total = sh + mono + di + tri;

                    if (total == 0) {
                        total = 1;
                    }

                    for (int x = max(start, bedStart) - bedStart; x < min(stop, bedStop) - bedStart; x++) {
                        double[] tr = tracks.get(x+offset);
                        if (tr == null) {
                            tracks.set(x+offset, new double[]{sh / total, mono / total, di / total, tri / total});
                        } else {
                            tr[0] += sh / total;
                            tr[1] += mono / total;
                            tr[2] += di / total;
                            tr[3] += tri / total;
                        }
                    }
                }
            }
            iter.close();
        }

        return averageAndScale(tracks);
    }

    /**
     * Function that will average the tracks in batches of 10
     * as well as scale them by given scaling factor
     * and then transforming them by taking the square root
     * @param tracks ArrayList containing the tracks to be averaged and scaled
     */
    private ArrayList<double[]> averageAndScale(ArrayList<double[]> tracks){
        // Reverse the ArrayList so we can loop it from the back to the front
        // This way we can remove each element in the loop before creating new elements in newTracks
        // This saves memory & removing from the back of an ArrayList is O(1) instead of O(n)
        Collections.reverse(tracks);
        ArrayList<double[]> newTracks = new ArrayList<>(tracks.size()/10);
        for (int i = tracks.size()-1; i >= 0; i-=10){
            double[] temp = new double[4];
            int end = Math.max(i-9, 0);
            double counter = Math.min(10, i+1);
            for (int a = i; a >= end; a--){
                double[] tr = tracks.remove(a);
                if (tr != null) {
                    temp[0] += tr[0];
                    temp[1] += tr[1];
                    temp[2] += tr[2];
                    temp[3] += tr[3];
                }
            }
            counter *= scale;
            for (int a = 0; a < temp.length; a++){
                temp[a] = Math.sqrt(temp[a] / counter);
            }
            newTracks.add(temp);
        }
        return newTracks;
    }

    /**
     * Returns a new Distribution
     * @param p a double specifying which distribution to use
     * @param m a double representing the mean for the desired distribution
     * @param l a double representing the standard deviation, or more generally, the lambda of the desired distribution
     * @return A new AbstractRealDistribution
     */
    private AbstractRealDistribution getDist(double p, double m, double l){
        if (p == 2){
            return new NormalDistribution(null, m, l);
        } else {
            return new ExponentialDistribution(null, m);
        }
    }
}
