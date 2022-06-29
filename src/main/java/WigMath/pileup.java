package WigMath;

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
import java.util.HashMap;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import FormatConverters.PileupToBedGraph;
import Node.TagNode;

import static HMMR_ATAC.Utils.toMap;

public class pileup {

	private final ArrayList<TagNode> intervals;
	private final File input;
	private final File index;
	private final int minMapQ;
	private final boolean rmDup;
	
	private final ArrayList<TagNode> bdg;
	private double cpmScale=0;
	
	public pileup(ArrayList<TagNode> t, File b, File i,int q,boolean r){
		intervals = t;
		input = b;
		index = i;
		minMapQ = q;
		rmDup=r;
		bdg = new ArrayList<>();
		build();
	}

	public double getCPMScale(){return cpmScale;}
	public ArrayList<TagNode> getBedGraph(){return bdg;}

	public void build(){

		HashMap<String,ArrayList<TagNode>> temp = toMap(intervals);
		for (String chr : temp.keySet()){
			ArrayList<TagNode> temp1 = temp.get(chr);
			temp1.sort(TagNode.basepairComparator);
			for (TagNode tagNode : temp1) {
				double[] block = makeBlock(tagNode);
				if (block.length > 0){
					bdg.addAll(new PileupToBedGraph(block, tagNode).getBedGraph());
				}
			}
		}
	}
	
	public double[] makeBlock(TagNode t){


		final SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(input).index(index));

		CloseableIterator<SAMRecord> iter = reader.query(t.getChrom(), t.getStart(), t.getStop(), false);

		if (!iter.hasNext()){
			return new double[0];
		}

		double[] temp = new double[t.getLength()];

		SAMRecord record;
		while (iter.hasNext()){
			record = iter.next();
			if(record != null){
				if (!record.getReadUnmappedFlag()
						&& !record.getMateUnmappedFlag()
						&& record.getMappingQuality() >= minMapQ
						&& !(record.getDuplicateReadFlag() && rmDup)
						&& Math.abs(record.getInferredInsertSize()) <= 1000
						&& record.getInferredInsertSize() != 0) {
					int readStart = record.getAlignmentStart();
					int readStop = record.getAlignmentEnd();

					cpmScale++;
					if (readStart < t.getStart()) {
						readStart = t.getStart();
					}
					if (readStop < t.getStart()) {
						readStop = t.getStart();
					}
					if (readStop >= t.getStop()) {
						readStop = t.getStop() - 1;
					}
					if (readStart >= t.getStop()) {
						readStart = t.getStop() - 1;
					}
					temp[readStart - t.getStart()]++;
					temp[readStop - t.getStart()]--;

				}
			}

		}
		iter.close();
		try {
			reader.close();
		} catch (IOException ignored) {}
		for (int i = 1; i < temp.length; i++){
			temp[i] = temp[i] + temp[i-1];
		}

		return temp;
	}
}
