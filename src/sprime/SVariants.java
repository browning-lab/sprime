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

import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import vcf.GeneticMap;
import blbutil.IntList;
import blbutil.SampleFileIt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import vcf.FilterUtil;
import vcf.Marker;
import vcf.PlinkGenMap;
import vcf.GTRec;
import vcf.VcfIt;

/**
 * <p>Class {@code SVariants} stores variant positions.
 * </p>
 * <p>Instances of class {@code SVariants} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SVariants {

    private final GeneticMap map;
    private final int[][] pos;
    private final double mu;    // global mutation rate / base pair / meiosis
    private final double globalDensity;

    /**
     * Constructs a new {@code SVariants} instance from the specified analysis
     * parameters.
     *
     * @param par the analysis parameters
     * @throws NullPointerException if {@code par == null}
     */
    public SVariants(SPar par) {
        List<IntList> list = new ArrayList<>(30);
        try (SampleFileIt<GTRec> it = vcfIt(par)) {
            while (it.hasNext()) {
                Marker m = it.next().marker();
                int chromIndex = m.chromIndex();
                while (chromIndex >= list.size()) {
                    list.add(new IntList(1000));
                }
                for (int j=1; j<m.nAlleles(); ++j) {
                    list.get(chromIndex).add(m.pos());
                }
            }
        }
        this.map = PlinkGenMap.fromPlinkMapFile(par.map());
        this.pos = convertToArray(list);
        this.mu = par.mu();
        this.globalDensity = globalDensity(pos);
    }

    private static SampleFileIt<GTRec> vcfIt(SPar par) {
        // want to read all chromosomes -- don't apply chrom filter
        Filter<String> sampleFilter = FilterUtil.sampleFilter(
                par.excludesamples());
        Filter<Marker> markerFilter = FilterUtil.markerFilter(
                par.excludemarkers());
        FileIt<String> it = InputIt.fromGzipFile(par.gt());
        SampleFileIt<GTRec> vcfIt
                = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toBitSetGT);
        return vcfIt;
    }

    private static int[][] convertToArray(List<IntList> list) {
        int[][] posMap = new int[list.size()][];
        for (int j=0; j<posMap.length; ++j) {
            posMap[j] = list.get(j).toArray();
        }
        return posMap;
    }

    /**
     * Returns the number of chromosomes used to determine the global
     * variant density, which is the number of ALT variants per base pair.
     * @return the number of chromosomes used to determine the global
     * variant density
     */
    public int nChrom() {
        return pos.length;
    }

    /**
     * Returns the estimated mutation rate per cM per meioses for the
     * specified closed chromosome interval.
     * @param chrom a chromosome index
     * @param startPos the starting base coordinate (inclusive)
     * @param inclEndPos the ending base coordinate (inclusive)
     * @return the estimated mutation rate per cM per meioses for the
     * specified closed chromosome interval
     * @throws IllegalArgumentException if {@code this.nVariants(chrom) == 0}
     * @throws IllegalArgumentException if {@code inclEndPos < startPos}
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= this.nChrom()}
     */
    public double mutPerCmPerGen(int chrom, int startPos, int inclEndPos) {
        if (inclEndPos < startPos) {
            throw new IllegalArgumentException(String.valueOf(inclEndPos));
        }
        if (pos[chrom].length==0) {
            throw new IllegalArgumentException(String.valueOf(chrom));
        }
        double cmPerBp = cmPerBp(chrom, startPos, inclEndPos);
        double bpPerCm = 1.0/cmPerBp;
        double localDensity = localDensity(chrom, startPos, inclEndPos);
        return (localDensity / globalDensity) * mu * bpPerCm;
    }

    private double cmPerBp(int chrom, int startPos, int inclEndPos) {
        final int STEP = 5000;
        final double MAX_CM = 0.01;
        int lastIndex = pos[chrom].length - 1;
        double minCmPerBp = Double.MAX_VALUE;
        double cm = 0.0;
        for (int n=0; (n<=20 || cm==0); ++n) {
            int p1 = Math.max(pos[chrom][0], startPos - n*STEP);
            int p2 = Math.min(pos[chrom][lastIndex], inclEndPos + n*STEP);
            double cm1 = map.genPos(chrom, p1);
            double cm2 = map.genPos(chrom, p2);
            cm = cm2 - cm1;
            if (cm > 0) {
                double cmPerBp = cm / (p2 - p1 + 1);
                if (cmPerBp < minCmPerBp) {
                    minCmPerBp = cmPerBp;
                }
            }
            if (cm>=MAX_CM && minCmPerBp<Double.MAX_VALUE) {
                break;
            }
        }
        if (minCmPerBp==Double.MAX_VALUE) {
            String s = "ERROR: local cmPerBp estimated to be 0.0";
            throw new IllegalArgumentException(s);
        }
        return minCmPerBp;
    }

    private double localDensity(int chrom, int startPos, int inclEndPos) {
        final int STEP = 5000;
        final int MIN_NVAR = 6;
        final int MAX_NVAR = 10;
        final int MAX_ITS = 20;
        int lastIndex = pos[chrom].length - 1;
        double maxDensity = -Double.MAX_VALUE;
        for (int n=0; n<=MAX_ITS; ++n) {
            int p1 = Math.max(pos[chrom][0], startPos - n*STEP);
            int p2 = Math.min(pos[chrom][lastIndex], inclEndPos + n*STEP);
            int nVar = nVariants(chrom, p1,  p2);
            if (nVar >= MIN_NVAR || n==MAX_ITS) {
                double density = (double) nVar / (p2 - p1 + 1);
                if (density >= maxDensity) {
                    maxDensity = density;
                }
            }
            if (nVar >= MAX_NVAR) {
                break;
            }
        }
        if (maxDensity == -Double.MAX_VALUE) {
            String s = "ERROR: Too few variants to estimate local variant density";
            throw new IllegalArgumentException(s);
        }
        return maxDensity;
    }

    private int nVariants(int chr, int startPos, int inclEndPos) {
        int i1 = Arrays.binarySearch(pos[chr], startPos);
        if (i1>=0) {
            while (i1>0 && pos[chr][i1]==pos[chr][i1-1]) {
                --i1;
            }
        }
        else {
            i1 = -i1 - 1;   // first variant in closed interval
        }
        int i2 = Arrays.binarySearch(pos[chr], i1, pos[chr].length, inclEndPos);
        if (i2>=0) {
            while ((i2+1)<pos[chr].length && pos[chr][i2]==pos[chr][i2+1]) {
                ++i2;
            }
        }
        else {
            i2 = -i2 - 2;   // last variant in closed interval
        }
        return i2 - i1 + 1; // number of varints in closed interval
    }

    private static double globalDensity(int[][] pos) {
        long bp = 0;
        long nVar = 0;
        for (int chr=0; chr<pos.length; ++chr) {
            int pos1 = pos[chr][0];
            int pos2 = pos[chr][pos[chr].length - 1];
            nVar += pos[chr].length;
            bp += pos2 - pos1 + 1;
        }
        return (double) nVar / bp;
    }
}
