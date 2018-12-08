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
import java.util.BitSet;

/**
 * <p>Class {@code BitSetGT} represents genotype emission
 * probabilities for a list of samples at a single marker.
 * The genotype emission probabilities are determined by the called
 * genotypes for the samples.
 * </p>
 * <p>Class {@code BitSetGT} stores alleles using
 * {@code java.util.BitSet} objects.
 * </p>
 * <p>Instances of class {@code BitSetGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitSetGT implements GTRec {

    /**
     * The VCF FORMAT code for genotype data: "GT".
     */
    public static final String GT_FORMAT = "GT";

    private final int bitsPerAllele;
    private final Marker marker;
    private final Samples samples;
    private final boolean isRefData;

    private final BitSet allele1;
    private final BitSet allele2;
    private final BitSet isMissing1;
    private final BitSet isMissing2;
    private final BitSet isPhased;

    /**
     * Constructs a new {@code BitSetGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param vcfHeader meta-information lines and header line for the
     * specified VCF record.
     * @param vcfRecord a VCF record corresponding to the specified
     * {@code vcfHeader} object
     *
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws IllegalArgumentException if {@code rec.nSamples() == 0}
     * @throws IllegalArgumentException if the header line
     * or VCF record does not have a "GT" format field
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRecord == null}
     */
    public BitSetGT(VcfHeader vcfHeader, String vcfRecord) {
        VcfRecGTParser gtp = new VcfRecGTParser(vcfHeader, vcfRecord);
        this.marker = gtp.marker();
        this.samples = vcfHeader.samples();
        this.bitsPerAllele = bitsPerAllele(marker);
        this.allele1 = new BitSet(vcfHeader.nSamples()*bitsPerAllele);
        this.allele2 = new BitSet(vcfHeader.nSamples()*bitsPerAllele);
        this.isMissing1 = new BitSet(vcfHeader.nSamples());
        this.isMissing2 = new BitSet(vcfHeader.nSamples());
        this.isPhased = new BitSet(vcfHeader.nSamples());
        gtp.storeAlleles(allele1, allele2, isMissing1, isMissing2, isPhased);
        this.isRefData = isRef(vcfHeader.nSamples(), isPhased, isMissing1, isMissing2);
    }

    private static boolean isRef(int nSamples, BitSet isPhased,
            BitSet isMissing1, BitSet isMissing2) {
        int nMissing = isMissing1.cardinality() + isMissing2.cardinality();
        int nUnphased = nSamples - isPhased.cardinality();
        return nMissing==0 && nUnphased==0;
    }

//    /**
//     * Constructs a new {@code LowMemGT} instance representing
//     * the specified VCF record's GT format field data.
//     *
//     * @param rec a VCF file record.
//     * @param fam parent-offspring relationships.
//     * @param usePhase {@code true} if phase information in the specified
//     * VCF file record will be used, and {@code false} if phase
//     * information in the specified VCF file record will be ignored.
//     *
//     * @throws IllegalArgumentException if
//     * {@code rec.nSamples()==0|| rec.samples().equals(fam.samples())==false}.
//     * @throws IllegalArgumentException if the VCF record does not have a
//     * GT format field.
//     * @throws NullPointerException if {@code rec==null || fam==null}.
//     */
//    public BitSetGT(VcfRecord rec, NuclearFamilies fam, boolean usePhase) {
//        this(rec);
//        if (rec.samples().equals(fam.samples())==false) {
//            throw new IllegalArgumentException("inconsistent samples");
//        }
//        setBits(rec, usePhase, bitsPerAllele, allele1, allele2, isMissing1,
//                isMissing2, isPhased);
//        removeMendelianInconsistencies(rec, fam, isPhased, isMissing1,
//                isMissing2);
//    }

    private BitSetGT(VcfRecord rec) {
        int nSamples = rec.nSamples();
        if (nSamples==0) {
            String s = "missing sample data: " + rec;
            throw new IllegalArgumentException(s);
        }
        if (rec.hasFormat(GT_FORMAT)==false) {
            String s = "missing GT FORMAT: " + rec;
            throw new IllegalArgumentException(s);
        }
        this.bitsPerAllele = bitsPerAllele(rec.marker());
        this.samples = rec.samples();
        this.marker = rec.marker();
        this.isRefData = isRef(rec);

        this.allele1 = new BitSet(nSamples * bitsPerAllele);
        this.allele2 = new BitSet(nSamples * bitsPerAllele);
        this.isMissing1 = new BitSet(nSamples);
        this.isMissing2 = new BitSet(nSamples);
        this.isPhased = new BitSet(nSamples);
    }

    private static boolean isRef(VcfRecord rec) {
        for (int j=0, n=rec.nSamples(); j<n; ++j) {
            if (rec.isPhased(j)==false || rec.allele1(j)<0 || rec.allele2(j)<0) {
                return false;
            }
        }
        return true;
    }

    private static void setBits(VcfRecord rec, boolean usePhase,
            int bitsPerAllele, BitSet allele1, BitSet allele2,
            BitSet isMissing1, BitSet isMissing2, BitSet isPhased) {
        int index1 = 0;
        int index2 = 0;
        for (int j=0, n=rec.nSamples(); j<n; ++j) {
            if (usePhase && rec.isPhased(j)) {
                isPhased.set(j);
            }
            int a1 = rec.allele1(j);
            int a2 = rec.allele2(j);
            if (a1 < 0) {
                isMissing1.set(j);
                index1 += bitsPerAllele;
            }
            else {
                int mask = 1;
                for (int k=0; k<bitsPerAllele; ++k) {
                    if ((a1 & mask)==mask) {
                        allele1.set(index1);
                    }
                    ++index1;
                    mask <<= 1;
                }
            }

            if (a2 < 0) {
                isMissing2.set(j);
                index2 += bitsPerAllele;
            }
            else {
                int mask = 1;
                for (int k=0; k<bitsPerAllele; ++k) {
                    if ((a2 & mask)==mask) {
                        allele2.set(index2);
                    }
                    ++index2;
                    mask <<= 1;
                }
            }
        }
    }

    private static int bitsPerAllele(Marker marker) {
        int nAllelesM1 = marker.nAlleles() - 1;
        int nStorageBits = Integer.SIZE - Integer.numberOfLeadingZeros(nAllelesM1);
        return nStorageBits;
    }

//    /*
//     * Sets phase to unknown for all parent-offspring relationships, and sets
//     * all genotypes in a duo or trio genotypes to missing if a Mendelian
//     * inconsistency is found.
//     */
//    private static void removeMendelianInconsistencies(VcfRecord rec,
//            NuclearFamilies fam, BitSet isPhased, BitSet isMissing1,
//            BitSet isMissing2) {
//        for (int j=0, n=fam.nDuos(); j<n; ++j) {
//            int p = fam.duoParent(j);
//            int o = fam.duoOffspring(j);
//            isPhased.clear(p);
//            isPhased.clear(o);
//            if (duoIsConsistent(rec, p, o) == false) {
//                logDuoInconsistency(rec, p, o);
//                isMissing1.set(p);
//                isMissing2.set(p);
//                isMissing1.set(o);
//                isMissing2.set(o);
//            }
//        }
//        for (int j=0, n=fam.nTrios(); j<n; ++j) {
//            int f = fam.trioFather(j);
//            int m = fam.trioMother(j);
//            int o = fam.trioOffspring(j);
//            isPhased.clear(f);
//            isPhased.clear(m);
//            isPhased.clear(o);
//            if (trioIsConsistent(rec, f, m, o) == false) {
//                logTrioInconsistency(rec, f, m, o);
//                isMissing1.set(f);
//                isMissing2.set(f);
//                isMissing1.set(m);
//                isMissing2.set(m);
//                isMissing1.set(o);
//                isMissing2.set(o);
//            }
//        }
//    }
//
//    private static boolean duoIsConsistent(VcfRecord rec, int parent,
//            int offspring) {
//        int p1 = rec.gt(parent, 0);
//        int p2 = rec.gt(parent, 1);
//        int o1 = rec.gt(offspring, 0);
//        int o2 = rec.gt(offspring, 1);
//        boolean alleleMissing = (p1<0 || p2<0 || o1<0 || o2<0);
//        return (alleleMissing || p1==o1 || p1==o2 || p2==o1 || p2==o2);
//    }
//
//    private static boolean trioIsConsistent(VcfRecord rec, int father,
//            int mother, int offspring) {
//        int f1 = rec.gt(father, 0);
//        int f2 = rec.gt(father, 1);
//        int m1 = rec.gt(mother, 0);
//        int m2 = rec.gt(mother, 1);
//        int o1 = rec.gt(offspring, 0);
//        int o2 = rec.gt(offspring, 1);
//        boolean fo1 = (o1<0 || f1<0 || f2<0 || o1==f1 || o1==f2);
//        boolean mo2 = (o2<0 || m1<0 || m2<0 || o2==m1 || o2==m2);
//        if (fo1 && mo2) {
//            return true;
//        }
//        else {
//            boolean fo2 = (o2<0 || f1<0 || f2<0 || o2==f1 || o2==f2);
//            boolean mo1 = (o1<0 || m1<0 || m2<0 || o1==m1 || o1==m2);
//            return (fo2 && mo1);
//        }
//    }
//
//    private static void logDuoInconsistency(VcfRecord rec, int parent,
//            int offspring) {
//        StringBuilder sb = new StringBuilder(80);
//        sb.append("WARNING: Inconsistent duo genotype set to missing");
//        sb.append(Const.tab);
//        sb.append(rec.marker());
//        sb.append(Const.colon);
//        sb.append(rec.samples().id(parent));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(offspring));
//        main.Logger.getInstance().println(sb.toString());
//    }
//
//    private static void logTrioInconsistency(VcfRecord rec, int father,
//            int mother, int offspring) {
//        StringBuilder sb = new StringBuilder(80);
//        sb.append("WARNING: Inconsistent trio genotype set to missing");
//        sb.append(Const.tab);
//        sb.append(rec.marker());
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(father));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(mother));
//        sb.append(Const.tab);
//        sb.append(rec.samples().id(offspring));
//        main.Logger.getInstance().println(sb.toString());
//    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nHaps() {
        return 2*samples.nSamples();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public boolean isPhased() {
        return isRefData;
    }

    @Override
    public boolean isGTData() {
        return true;
    }

    @Override
    public boolean isPhased(int sample) {
        return isPhased.get(sample);
    }

    @Override
    public int allele1(int sample) {
        return isMissing1.get(sample) ? -1 : allele(allele1, sample);
    }

    @Override
    public int allele2(int sample) {
        return isMissing2.get(sample) ? -1 : allele(allele2, sample);
    }

    @Override
    public float gl(int sample, int a1, int a2) {
        if (a1 < 0 || a1 >= marker.nAlleles())  {
            String s = "invalid alleles: (" + a1 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        if (a2 < 0 || a2 >= marker.nAlleles()) {
            String s = "invalid alleles: (" + a2 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        int obsA1 = allele1(sample);
        int obsA2 = allele2(sample);
        boolean consistent = (obsA1==-1 || obsA1==a1) && (obsA2==-1 || obsA2==a2);
        if (consistent==false && isPhased(sample)==false) {
            consistent = (obsA1==-1 || obsA1==a2) && (obsA2==-1 || obsA2==a1);
        }
        return consistent ? 1.0f : 0.0f;
    }

    private int allele(BitSet bits, int sample) {
        int start = bitsPerAllele*sample;
        int end = start + bitsPerAllele;
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (bits.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    @Override
    public int nAlleles() {
        return this.marker().nAlleles();
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return toVcfRec();
    }
}
