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

import java.util.Arrays;
import vcf.Marker;

/**
 * <p>Class {@code SPairScorer} stores allele dose data for a chromosome
 * and computes pairwise scores for all pairs of dose records.
 * </p>
 * Instances of class {@code SPairScorer} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SPairScorer {

    private static final int MIN_DIST = 10;
    private static final int MAX_DIST = 20_000;

    private final DoseRec[] recs;
    private final int[] start;      // inclusive start
    private final int[] inclEnd;    // inclusive end

    /**
     * Constructs a new {@code SPairScores} instance from the specified data.
     * The list of {@code DoseRec} instances must be derived from the same
     * chromosome.
     * @param recs a sorted list of {@code DoseRec} instances from a single
     * chromosome
     * @throws IllegalArgumentException if
     * {@code (recs[j].marker().chromIndex != recs[k].marker().chromIndex())}
     * for any {@code j, k} satisfying
     * {@code (0 <= j && j < k && k < recs.length)}
     * @throws NullPointerException if {@code recs == null}
     * @throws NullPointerException if {@code (recs[j] == null)} for some
     * {@code j} satisfying {@code (0 <= j && j < recs.length)}
     */
    public SPairScorer(DoseRec[] recs) {
        boolean isSorted = checkData(recs);
        DoseRec[] clone = recs.clone();
        if (isSorted==false) {
            Arrays.sort(clone);
        }
        this.recs = clone;
        this.start = new int[recs.length];
        this.inclEnd = new int[recs.length];
        for (int i=0; i<recs.length; ++i) {
            int pos = recs[i].marker().pos();
            int maxPos = pos - MIN_DIST;
            int minPos = pos - MAX_DIST;
            int j = i;
            while (j>0 && (recs[j-1].pos() > maxPos)) {
                --j;
            }
            if (j==0 || recs[j-1].pos() < minPos) {
                start[i] = -1;
                inclEnd[i] = -1;
            }
            else {
                --j;
                assert (minPos <= recs[j].pos() && recs[j].pos() <= maxPos);
                inclEnd[i] =  j;
                while (j>0 && (recs[j-1].pos() >= minPos)) {
                    --j;
                }
                start[i] = j;
            }
        }
    }

    private boolean checkData(DoseRec[] recs) {
        boolean isSorted = true;
        for (int j=1; j<recs.length; ++j) {
            Marker m1 = recs[j-1].marker();
            Marker m2 = recs[j].marker();
            if (m1.chromIndex() != m2.chromIndex()) {
                throw new IllegalArgumentException(m2.chrom());
            }
            if (m1.pos() > m2.pos()) {
                isSorted = false;
            }
        }
        return isSorted;
    }

    /**
     * Returns the minimum permitted distance between two consecutive allele
     * dose records in an introgressed segment.
     * @return the minimum permitted distance between two consecutive allele
     * dose records in an introgressed segment
     */
    public static int minDist() {
        return MIN_DIST;
    }

    /**
     * Returns the maximum permitted distance between two consecutive allele
     * dose records in an introgressed segment.
     * @return the maximum permitted distance between two consecutive allele
     * dose records in an introgressed segment
     */
    public static int maxDist() {
        return MAX_DIST;
    }

    /**
     * Returns the number of records.
     * @return the number of records
     */
    public int nRecs() {
        return recs.length;
    }

    /**
     * Returns the specified {@code DoseRec}.
     * @param index a variant index
     * @return the specified {@code DoseRec}
     * @throws IllegalArgumentException if
     * {@code index < 0 || index >= this.nRecs()}
     */
    public DoseRec rec(int index) {
        return recs[index];
    }

    /**
     * Returns the specified pair score.
     * @param i1 the smaller variant index
     * @param i2 the larger variant index
     * @param mutPerCm estimated mutation rate per cM per meioses for the
     * interval between the two allele dose records
     * @return the specified pair score
     * @throws IllegalArgumentException if {@code i1 >= i2}
     * @throws IllegalArgumentException if
     * {@code Double.isFinite(mutPerCm) == false || mutPerCm <= 0.0}
     * @throws IndexOutOfBoundsException if
     * {@code i1 < 0 || i2 >= this.nRecs()}
     */
    public double score(int i1, int i2, double mutPerCm) {
        if (i1 < 0) {
            throw new IndexOutOfBoundsException(String.valueOf(i1));
        }
        if (i1 >= i2) {
            throw new IllegalArgumentException(String.valueOf(i2));
        }
        if (Double.isFinite(mutPerCm)==false || mutPerCm<=0) {
            throw new IllegalArgumentException(String.valueOf(mutPerCm));
        }
        if (i1 < start[i2] || i1 > inclEnd[i2]) {
            return Double.NEGATIVE_INFINITY;
        }
        else {
            int maxDistance = recs[i1].targCnt() + recs[i2].targCnt();
            double d = DoseRec.distance(recs[i1], recs[i2]);
            if (d==maxDistance) {
                // no sample carries both alleles
                return Double.NEGATIVE_INFINITY;
            }
            else {
                double n = Math.min(recs[i1].targCnt(), recs[i2].targCnt());
                double firstTerm = 6000.0*(1-Math.exp(-1.0/(mutPerCm*100)))/(1-Math.exp(-1.0));
                if (recs[i2].outgroupCnt()>0) {
                    firstTerm *= 0.80;
                }
                return firstTerm - 25000*d/n;
            }
        }
    }

    /**
     * Returns the smallest index {@code i} of an allele dose record which
     * satisfies
     * {@code ( (this.rec(i).pos() >= (this.rec(index).pos() - SPairScorer.maxDist()))
     * && (this.rec(i).pos() <= (this.rec(index).pos() - SPairScorer.minDist())) )}.
     * Returns {@code -1} if no such records exist.
     * @param index a variant index
     * @return the specified index
     * @throws IllegalArgumentException if
     * {@code index < 0 || index >= this.nRecs()}
     */
    public int start(int index) {
        return start[index];
    }

    /**
     * Returns the largest index {@code i} of an allele dose record which
     * satisfies
     * {@code ( (this.rec(i).pos() >= (this.rec(index).pos() - SPairScorer.maxDist()))
     * && (this.rec(i).pos() <= (this.rec(index).pos() - SPairScorer.minDist())) )}.
     * Returns {@code -1} if no such records exist.
     * @param index a variant index
     * @return the specified index
     * @throws IllegalArgumentException if
     * {@code index < 0 || index >= this.nRecs()}
     */
    public int inclEnd(int index) {
        return inclEnd[index];
    }
}
