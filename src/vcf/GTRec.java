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
package vcf;

import beagleutil.Samples;
import blbutil.Const;

/**
 * <p>Interface {@code GTRec} represents represents genotype data for one
 * marker.
 * </p>
 * <p>All instances of {@code GTRec} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GTRec extends DuplicatesGTRec {

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Returns {@code true} if the value returned by {@code this.gl()} is
     * determined by a called or missing genotype, and returns {@code false}
     * otherwise.
     * @return {@code true} if the value returned by {@code this.gl()} is
     * determined by a called or missing genotype
     *
     * @implSpec The default implementation returns {@code true}
     */
    default boolean isGTData() {
        return true;
    }

    /**
     * Returns the probability of the observed data for the specified sample
     * if the specified pair of ordered alleles is the true ordered genotype.
     * @param sample the sample index
     * @param allele1 the first allele index
     * @param allele2 the second allele index
     * @return the probability of the observed data for the specified sample
     * if the specified pair of ordered alleles is the true ordered genotype.
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || samples >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker().nAlleles()}
     *
     * @implSpec The default implementation returns {@code 1.0f} if the
     * corresponding genotype determined by the {@code isPhased()},
     * {@code allele1()}, and {@code allele2()} methods is consistent
     * with the specified ordered genotype, and returns {@code 0.0f} otherwise.
     */
    default float gl(int sample, int allele1, int allele2) {
        int nAlleles = this.marker().nAlleles();
        if (allele1 < 0 || allele1 >= nAlleles)  {
            String s = "invalid alleles: (" + allele1 + "): " + marker();
            throw new IllegalArgumentException(s);
        }
        if (allele2 < 0 || allele2 >= nAlleles) {
            String s = "invalid alleles: (" + allele2 + "): " + marker();
            throw new IllegalArgumentException(s);
        }
        int a1 = this.allele1(sample);
        int a2 = this.allele2(sample);
        boolean consistent = (a1==-1 || a1==allele1) && (a2==-1 || a2==allele2);
        if (consistent==false && this.isPhased(sample)==false) {
            consistent = (a1==-1 || a1==allele2) && (a2==-1 || a2==allele1);
        }
        return consistent ? 1.0f : 0.0f;
    }

    /**
     * Returns a VCF record corresponding to {@code this}.  The returned
     * VCF record will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return a VCF record corresponding to {@code this}
     */
    default String toVcfRec() {
        StringBuilder sb = new StringBuilder(100);
        sb.append(marker());
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // QUAL
        sb.append(Const.tab);
        sb.append("PASS");                   // FILTER
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // INFO
        sb.append(Const.tab);
        sb.append("GT");                     // FORMAT
        for (int j=0, n=nSamples(); j<n; ++j) {
            int a1 = allele1(j);
            int a2 = allele2(j);
            sb.append(Const.tab);
            if (a1==-1) {
                sb.append(Const.MISSING_DATA_CHAR);
            }
            else {
                sb.append(a1);
            }
            sb.append(GTRec.this.isPhased(j) ? Const.phasedSep : Const.unphasedSep);
            if (a2==-1) {
                sb.append(Const.MISSING_DATA_CHAR);
            }
            else {
                sb.append(a2);
            }
        }
        return sb.toString();
    }
}