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

import blbutil.IntSet;
import vcf.Marker;
import vcf.GTRec;

/**
 * <p>Class {@code DoseRec} represents an allele dose record.  The
 * {@code DoseRec} class has a natural ordering that is inconsistent with equals.
 * </p>
 * Instances of class {@code DoseRec} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DoseRec implements Comparable<DoseRec> {

    private final Marker marker;
    private final int allele;
    private final IntSet hets;
    private final IntSet homs;
    private final int targCnt;
    private final int outGroupCnt;
    private final float outGroupFreq;

    /**
     * Constructs a new {@code DoseRec} instance from the specified data.
     * @param rec the genotype data
     * @param allele the allele
     * @param inOutgroup an array whose {@code true} elements identify the indices
     * of outgroup samples
     * @throws IllegalArgumentException if
     * {@code allele < 0 || allele >= rec.marker().nAlleles()}
     * @throws IllegalArgumentException if
     * {@code inOutgroup.length != rec.nSamples()}
     * @throws NullPointerException if {@code rec == null || inOutgroup == null}
     */
    public DoseRec(GTRec rec, int allele, boolean[] inOutgroup) {
        if (allele < 0 || allele >= rec.marker().nAlleles()) {
            throw new IllegalArgumentException(String.valueOf(allele));
        }
        if (inOutgroup.length!=rec.nSamples()) {
            throw new IllegalArgumentException(String.valueOf(inOutgroup.length));
        }
        this.marker = rec.marker();
        this.allele = allele;
        this.hets = new IntSet(10);
        this.homs = new IntSet(4);
        int outGrpAlleleCnt = 0;
        int outGrpNonMissingCnt = 0;
        for (int s=0; s<inOutgroup.length; ++s) {
            int a1 = rec.allele1(s);
            int a2 = rec.allele2(s);
            int dose = dose(a1, a2, allele);
            int nMissingAlleles = dose(a1, a2, -1);

            if (inOutgroup[s]) {
                outGrpAlleleCnt += dose;
                outGrpNonMissingCnt += (2 - nMissingAlleles);
            }
            else {
                switch (dose) {
                    case 0: break;
                    case 1: hets.add(s); break;
                    case 2: homs.add(s); break;
                    default: throw new IllegalStateException(String.valueOf(dose));
                }
            }
        }
        this.targCnt = hets.size() + 2*homs.size();
        this.outGroupCnt = outGrpAlleleCnt;
        this.outGroupFreq = (float) outGrpAlleleCnt / outGrpNonMissingCnt;
    }

    private static int dose(int a1, int a2, int allele) {
        int dose = 0;
        if (a1==allele) {
            ++dose;
        }
        if (a2==allele) {
            ++dose;
        }
        return dose;
    }

    /**
     * Returns the sum over all target samples of the absolute value of
     * the allele dose difference between the two specified {@code DoseRec}
     * objects.
     * @param a a variant dose record
     * @param b a variant dose record
     * @return the sum over all target samples of the absolute value of
     * the allele dose difference
     * @throws NullPointerException if {@code a == null || b == null}
     */
    public static int distance(DoseRec a, DoseRec b) {
        int cnt = 0;
        for (int j=0, n=a.hets.size(); j<n; ++j) {
            int sample = a.hets.elementWithIndex(j);
            if (b.hets.contains(sample)==false) {
                ++cnt;
            }
        }
        for (int j=0, n=a.homs.size(); j<n; ++j) {
            int sample = a.homs.elementWithIndex(j);
            if (b.hets.contains(sample)) {
                ++cnt;
            }
            else if (b.homs.contains(sample)==false) {
                cnt+=2;
            }
        }
        for (int j=0, n=b.hets.size(); j<n; ++j) {
            int sample = b.hets.elementWithIndex(j);
            if (a.hets.contains(sample)==false && a.homs.contains(sample)==false) {
                ++cnt;
            }
        }
        for (int j=0, n=b.homs.size(); j<n; ++j) {
            int sample = b.homs.elementWithIndex(j);
            if (a.hets.contains(sample)==false && a.homs.contains(sample)==false) {
                cnt+=2;
            }
        }
        return cnt;
    }

    /**
     * Returns the chromosome position.
     * @return the chromosome position
     */
    public int pos() {
        return marker.pos();
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns the allele whose dose is stored.
     * @return the allele whose dose is stored
     */
    public int allele() {
        return allele;
    }

    /**
     * Returns the number of copies of {@code this.allele()} in the target
     * samples.
     * @return the number of copies of {@code this.allele()} in the target
     * samples
     */
    public int targCnt() {
        return targCnt;
    }

    /**
     * Returns the number of copies of {@code this.allele()} in the outgroup
     * samples.
     * @return the number of copies of {@code this.allele()} in the outgroup
     * samples
     */
    public int outgroupCnt() {
        return outGroupCnt;
    }

    /**
     * Returns the allele frequency of {@code this.allele()} in the outgroup
     * samples.
     * @return the allele frequency of {@code this.allele()} in the outgroup
     * samples
     */
    public float outgroupFreq() {
        return outGroupFreq;
    }

    /**
     * Compares the specified {@code DoseRec} with this for order, and
     * returns a negative integer, zero, or a positive integer respectively
     * if this object is less than, equal to, or greater than the specified
     * object. Two {@code DoseRec} objects are first compared on
     * {@code this.marker()} and then on {@code this.allele()}.
     * @param rec the record to be compared withe {@code this}
     * @return a negative integer, zero, or a positive integer respectively
     * if this object is less than, equal to, or greater than the specified
     * object
     * @throws NullPointerException if {@code rec == null}
     */
    @Override
    public int compareTo(DoseRec rec) {
        int x = this.marker().compareTo(rec.marker());
        if (x==0) {
            if (this.allele()!=rec.allele()) {
                return this.allele() < rec.allele() ? -1 : 1;
            }
        }
        return 0;
    }
}
