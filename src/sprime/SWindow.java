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

import beagleutil.ChromInterval;
import beagleutil.Samples;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.stream.IntStream;
import vcf.FilterUtil;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.GTRec;
import vcf.VcfIt;

/**
 * <p>Class {@code SWindow} iterates through the chromosomes of data in an
 * input VCF file.
 * </p>
 * <p>Instances of class {@code SWindow} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SWindow implements FileIt<DoseRec[]> {

    private final File file;
    private final SampleFileIt<GTRec> it;
    private final boolean[] inOutgroup;
    private final int[] refIndices;
    private final int[] targIndices;
    private final int maxCnt;

    private GTRec next;

    /**
     * Constructs a new {@code SWindow} instance from the specified data.
     * @param par the SPrime analysis parameters
     * @throws NullPointerException if {@code par == null}
     */
    public SWindow(SPar par) {
        this.file = par.gt();
        this.it = vcfIt(par);
        this.inOutgroup = inOutgroup(it.samples(), par.outgroup());
        this.refIndices = indices(inOutgroup, true);
        this.targIndices = indices(inOutgroup, false);
        this.maxCnt = (int) Math.floor(par.maxfreq()*refIndices.length);
        this.next = it.hasNext() ? it.next() : null;
    }

    private static SampleFileIt<GTRec> vcfIt(SPar par) {
        Filter<String> sampleFilter = FilterUtil.sampleFilter(
                par.excludesamples());
        Filter<Marker> markerFilter = FilterUtil.markerFilter(
                par.excludemarkers());
        ChromInterval chromInterval = ChromInterval.parse(par.chrom());
        FileIt<String> it = InputIt.fromGzipFile(par.gt());
        SampleFileIt<GTRec> vcfIt
                = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toBitSetGT);
        if (chromInterval!=null) {
            vcfIt = new IntervalVcfIt<>(vcfIt, chromInterval);
        }
        return vcfIt;
    }

    private static boolean[] inOutgroup(Samples samples, File outgroupFile) {
        Set<String> outgroup = Utilities.idSet(outgroupFile);
        boolean[] ba = new boolean[samples.nSamples()];
        for (int j=0; j<ba.length; ++j) {
            ba[j] = outgroup.contains(samples.id(j));
        }
        return ba;
    }

    private static int[] indices(boolean[] ba, boolean value) {
        return IntStream.range(0, ba.length)
                .filter(j -> ba[j]==value)
                .toArray();
    }

    /**
     * Returns the number of outgroup samples
     * @return the number of outgroup samples
     */
    public int nOutgroupSamples() {
        return refIndices.length;
    }

    /**
     * Returns the number of target samples
     * @return the number of target samples
     */
    public int nTargetSamples() {
        return targIndices.length;
    }

    @Override
    public File file() {
        return file;
    }

    @Override
    public void close() {
        next = null;
        it.close();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return next!=null;
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public DoseRec[] next() {
        assert next!=null;
        List<DoseRec> list = new ArrayList<>();
        int chromIndex = next.marker().chromIndex();
        boolean sameChrom = true;
        while (next!=null && sameChrom) {
            checkNoMissingAlleles(next);
            processRec(next, list);
            if (it.hasNext()) {
                next = it.next();
                sameChrom = next.marker().chromIndex()==chromIndex;
            }
            else {
                next = null;
            }
        }
        return list.toArray(new DoseRec[0]);
    }

    private static void checkNoMissingAlleles(GTRec rec) {
        for (int h=0; h<rec.nHaps(); ++h) {
            if (rec.allele(h) == -1) {
                Marker m = rec.marker();
                System.out.println(Const.nl + "ERROR: VCF record has missing alleles:"
                        + " CHROM=" + m.chrom() + " POS=" + m.pos());
                System.exit(1);
            }
        }
    }

    private void processRec(GTRec rec, List<DoseRec> list) {
        int[] alCnts = outgroupAlleleCnts(rec);
        for (int al=0; al<alCnts.length; ++al) {
            if (alCnts[al] <= maxCnt) {
                list.add(new DoseRec(rec, al, inOutgroup));
            }
        }
    }

    private int[] outgroupAlleleCnts(GTRec ve) {
        int[] alCnts = new int[ve.marker().nAlleles()];
        for (int j=0; j<refIndices.length; ++j) {
            int a1 = ve.allele1(refIndices[j]);
            int a2 = ve.allele2(refIndices[j]);
            ++alCnts[a1];
            ++alCnts[a2];
        }
        return alCnts;
    }
}
