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

/**
 * <p>Class {@code OurRec} represents a variant in an introgressed segment.
 * Class {@code OurRec} has a natural ordering that is not consistent with
 * equals.
 * </p>
 * <p>Instances of {@code SOutRec} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SOutRec implements Comparable<SOutRec> {

    private final DoseRec rec;
    private final int segment;
    private final double score;

    /**
     * Constructs a new {@code SOutRec} object from the specified data.
     * @param rec an allele dose record
     * @param segment the index of the introgressed segment containing the
     * allele dose record
     * @param score the introgressed segment score
     * @throws IllegalArgumentException if
     * {@code Double.isFinite(score) == false}
     * @throws NullPointerException if {@code rec == null}
     */
    public SOutRec(DoseRec rec, int segment, double score) {
        if (rec==null) {
            throw new NullPointerException(DoseRec.class.toString());
        }
        if (Double.isFinite(score)==false) {
            throw new IllegalArgumentException(String.valueOf(score));
        }
        this.rec = rec;
        this.segment = segment;
        this.score = score;
    }

    /**
     * Returns the allele dose record.
     * @return the allele dose record
     */
    public DoseRec rec() {
        return rec;
    }

    /**
     * Returns the index of the introgressed segment containing the allele dose
     * record.
     * @return the index of the introgressed segment containing the allele dose
     * record
     */
    public int segment() {
        return segment;
    }

    /**
     * Returns the introgressed segment score.
     * @return the introgressed segment score
     */
    public double score() {
        return score;
    }

    /**
     * Compares the specified {@code SOutRec} with {@code this} for order, and
     * returns a negative integer, zero, or a positive integer respectively
     * if this object is less than, equal to, or greater than the specified
     * object.  Two {@code SOutRec} objects are first compared on
     * {@code this.rec()}, then on {@code this.segment()}, and finally on
     * {@code this.score()}
     * @param other the {@code SOutRec} to be compared with {@code this}
     * @return a negative integer, zero, or a positive integer respectively
     * if this object is less than, equal to, or greater than the specified
     * object
     * @throws NullPointerException if {@code other == null}
     */
    @Override
    public int compareTo(SOutRec other) {
        int x = this.rec().marker().compareTo(other.rec().marker());
        if (x != 0) {
            return x;
        }
        if (x != 0) {
            if (this.segment() != other.segment()) {
                return this.segment() < other.segment() ? -1 : 1;
            }
            else {
                return Double.compare(this.score(), other.score());
            }
        }
        return 0;
    }

    /**
     * Returns a string representation of {@code this}.  The format
     * is unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(rec.marker());
        sb.append(Const.tab);
        sb.append(segment);
        sb.append(Const.tab);
        sb.append(rec.allele());
        sb.append(Const.tab);
        sb.append((int) Math.round(score));
        return sb.toString();
    }

    /**
     * Returns the header line for SPrime output.  The format of the header
     * line is unspecified and subject to change.
     * @return the header line for SPrime output.
     */
    public static String header() {
        StringBuilder sb = new StringBuilder();
        sb.append("CHROM");
        sb.append(Const.tab);
        sb.append("POS");
        sb.append(Const.tab);
        sb.append("ID");
        sb.append(Const.tab);
        sb.append("REF");
        sb.append(Const.tab);
        sb.append("ALT");
        sb.append(Const.tab);
        sb.append("SEGMENT");
        sb.append(Const.tab);
        sb.append("ALLELE");
        sb.append(Const.tab);
        sb.append("SCORE");
        return sb.toString();
    }
}
