/*
 * Copyright 2017-2018 Brian L. Browning
 *
 * This file is part of the sprime program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package sprime;

import blbutil.Const;
import blbutil.FileUtil;
import blbutil.IntList;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * <p>Class {@code SAnalyzer} performs an SPrime analysis.
 * </p>
 * <p>Instances of class {@code SAnalyzer} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SAnalyzer {

    private final SPar par;
    private final File outFile;
    private final SVariants variants;

    private int nOutgroupSamples = 0;
    private int nTargetSamples = 0;

    private int segmentCnt = 0;
    private int variantCnt = 0;

    /**
     * Constructs a new {@code SAnalyzer} class for the specified analysis
     * parameters.
     * @param par the SPrime analysis parameters
     * @throws NullPointerException if {@code par == null}
     */
    public SAnalyzer(SPar par) {
        this.par = par;
        this.variants = new SVariants(par);
        this.outFile = new File(par.out() + ".score");
    }

    /**
     * Runs an SPrime analysis.
     */
    public void run() {
        try (PrintWriter out = FileUtil.printWriter(outFile)) {
            out.println(SOutRec.header());
            List<SOutRec> outRecs = new ArrayList<>(1000);
            try (SWindow window = new SWindow(par)) {
                this.nOutgroupSamples = window.nOutgroupSamples();
                this.nTargetSamples = window.nTargetSamples();
                while (window.hasNext()) {
                    DoseRec[] inRecs = window.next();
                    variantCnt += inRecs.length;
                    analyzeWindow(inRecs, outRecs);
                }
            }
            outRecs.stream()
                    .sorted()
                    .forEach(rec -> out.println(rec.toString()));
            outRecs.clear();
        }
    }

    private void analyzeWindow(DoseRec[] inRecs,  List<SOutRec> outRecs) {
        int nRecs = inRecs.length;
        double[] scores = new double[nRecs];
        int[] prevMarker = new int[nRecs];
        BitSet changed = new BitSet(nRecs);
        SPairScorer scorer = new SPairScorer(inRecs);
        for (int j=0; j<nRecs; ++j) {
            setScore(scorer, j, scores, prevMarker);
        }
        int[] indices = storeMaxSeg(inRecs, scores, prevMarker, outRecs);
        while (indices.length>0) {
            excludeIndices(indices, scores, prevMarker, changed);
            int lastChangedPos = inRecs[indices[indices.length-1]].pos();
            boolean finished = false;
            for (int i=indices[0] + 1; i<nRecs && finished==false; ++i) {
                if (rescore(i, scores, prevMarker, changed))  {
                    setScore(scorer, i, scores, prevMarker);
                    changed.set(i);
                    if (inRecs[i].pos() > lastChangedPos) {
                        lastChangedPos = inRecs[i].pos();
                    }
                }
                finished = (inRecs[i].pos() - lastChangedPos) > SPairScorer.maxDist();
            }
            indices = storeMaxSeg(inRecs, scores, prevMarker, outRecs);
        }
    }

    private void setScore(SPairScorer scorer, int index, double[] scores,
            int[] prevMarker) {
        int start = scorer.start(index);
        int inclEnd = scorer.inclEnd(index);
        scores[index] = 0.0;
        prevMarker[index] = -1;
        if (start != -1) {
            for (int k=start; k<=inclEnd; ++k) {
                if (scores[k] >= 0.0) {
                    int chrom = scorer.rec(k).marker().chromIndex();
                    int p1 = scorer.rec(k).marker().pos();
                    int p2 = scorer.rec(index).marker().pos();
                    double mutPerCm = variants.mutPerCmPerGen(chrom, p1, p2);
                    double score = scores[k] + scorer.score(k, index, mutPerCm);
                    if (score > scores[index]) {
                        scores[index] = score;
                        prevMarker[index] = k;
                    }
                }
            }
        }
    }

    private int[] storeMaxSeg(DoseRec[] inRecs, double[] scores,
            int[] prevMarker, List<SOutRec> outRecs) {
        IntList indexList = new IntList(100);
        int index = maxIndex(scores);
        double score = scores[index];
        if (score >= par.minscore()) {
            outRecs.add(new SOutRec(inRecs[index], segmentCnt, score));
            indexList.add(index);
            while (prevMarker[index] != -1) {
                index = prevMarker[index];
                outRecs.add(new SOutRec(inRecs[index], segmentCnt, score));
                indexList.add(index);
            }
            ++segmentCnt;
        }
        int[] indices = indexList.toArray();
        Arrays.sort(indices);
        return indices;
    }

    private void excludeIndices(int[] indices, double[] scores, int[] prevMarker,
            BitSet changed) {
        changed.clear();
        for (int i : indices) {
            scores[i] = Double.NEGATIVE_INFINITY;
            prevMarker[i] = -1;
            changed.set(i);
        }
    }

    private static boolean rescore(int index, double[] scores, int[] prevMarker,
            BitSet modified) {
        int prev = prevMarker[index];
        return scores[index]>=0 && (prev != -1) && modified.get(prev);
    }

    private static int maxIndex(double[] scores) {
        double max = Double.NEGATIVE_INFINITY;
        int maxIndex = 0;
        for (int j=0; j<scores.length; ++j) {
            if (scores[j]>=max) {
                max = scores[j];
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public SPar par() {
        return par;
    }

    /**
     * Returns the number of chromosomes used to determine the global
     * variant density, which is the number of ALT variants per base pair.
     * Returns 0 if {@code this.run()} has not been invoked.
     * @return the number of chromosomes used to determine the global
     * variant density
     */
    public int nChrom() {
        return variants.nChrom();
    }

    /**
     * Returns the number of outgroup samples. Returns 0 if {@code this.run()}
     * has not been invoked.
     * @return the number of outgroup samples
     */
    public int nOutgroupSamples() {
        return nOutgroupSamples;
    }

    /**
     * Returns the number of target samples. Returns 0 if {@code this.run()}
     * has not been invoked.
     * @return the number of target samples
     */
    public int nTargetSamples() {
        return nTargetSamples;
    }

    /**
     * Returns the number of analyzed VCF variants, which are variants with
     * frequency less than or equal to {@code this.par().maxfreq()}.
     * Returns 0 if {@code this.run()} has not been invoked.
     *
     * @return the number of analyzed VCF variants
     */
    public int variantCnt() {
        return variantCnt;
    }

    /**
     * Returns the number of number of introgressed segments identified.
     * Returns 0 if {@code this.run()} has not been invoked.
     *
     * @return the number of number of introgressed segments identified
     */
    public int segmentCnt() {
        return segmentCnt;
    }

    private void dbgPrint(File dbgFile, DoseRec[] recs, double[] scores,
            int[] prevMarker) {
        try (PrintWriter dbg = FileUtil.printWriter(dbgFile)) {
            dbg.print("Index");
            dbg.print(Const.tab);
            dbg.print("Pos");
            dbg.print(Const.tab);
            dbg.print("Allele");
            dbg.print(Const.tab);
            dbg.print("Score");
            dbg.print(Const.tab);
            dbg.println("Previous");
            for (int j=0; j<recs.length; ++j) {
                dbg.print(j);
                dbg.print(Const.tab);
                dbg.print(recs[j].pos());
                dbg.print(Const.tab);
                dbg.print(recs[j].allele());
                dbg.print(Const.tab);
                dbg.print((int) Math.round(scores[j]));
                dbg.print(Const.tab);
                dbg.println(prevMarker[j]);
            }
        }
    }
}
