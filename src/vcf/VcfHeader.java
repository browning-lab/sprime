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
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.StringUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code VcfHeader} represents the Variant Call Format (VCF)
 * meta-information lines and the Variant Call Format header line
 * that precede the first Variant Call Format record.
 * </p>
 * <p>Instances of class {@code VcfHeader} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfHeader  {

    private static final String SHORT_HEADER_PREFIX= "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO";

    /**
     * A string equal to the first nine tab-delimited fields of a VCF header
     * line that contains sample data.
     */
    public static final String HEADER_PREFIX =
            SHORT_HEADER_PREFIX + Const.tab + "FORMAT";

    private static final int sampleOffset = 9;

    private final File file;   // null if source is standard input
    private final VcfMetaInfo[] metaInfoLines;
    private final String headerLine;
    private final int nHeaderFields;
    private final int[] includedIndices;
    private final Samples samples;

    /**
     * Constructs a new {@code VcfHeader} object from the VCF
     * meta-information lines and the VCF header line returned by the
     * specified {@code FileIterator<String>}.  This constructor will advance
     * the {@code FileIterator<String>} to the point before the first VCF record
     * in the file.  The {@code VcfHeader} object will have no excluded samples.
     * @param it an iterator that returns lines of VCF file
     *
     * @throws IllegalArgumentException if any of the meta-information lines
     * returned by the specified {@code FileIterator<String>} does not conform
     * to the VCF specification
     * @throws IllegalArgumentException if the header lines returned by the
     * specified {@code FileIterator<String>} does not conform to the VCF
     * specification
     * @throws IllegalArgumentException if no header line is returned by the
     * specified {@code FileIterator<String>}
     *
     * @throws NullPointerException if {@code it == null}
     */
    public VcfHeader(FileIt<String> it) {
        this(it, Filter.acceptAllFilter());
    }
    /**
     * Constructs a new {@code VcfHeader} object from the VCF
     * meta-information lines and the VCF header line returned by the
     * specified {@code FileIterator<String>}.  This constructor will advance
     * the {@code FileIterator<String>} to the point before the first VCF record in the file.
     * @param it an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     *
     * @throws IllegalArgumentException if any of the meta-information lines
     * returned by the specified {@code FileIterator<String>} does not conform
     * to the VCF specification
     * @throws IllegalArgumentException if the header lines returned by the
     * specified {@code FileIterator<String>} does not conform to the VCF
     * specification
     * @throws IllegalArgumentException if no header line is returned by the
     * specified {@code FileIterator<String>}
     *
     * @throws NullPointerException if {@code it == null}
     */
    public VcfHeader(FileIt<String> it, Filter<String> sampleFilter) {
        if (sampleFilter==null) {
            sampleFilter = Filter.acceptAllFilter();
        }
        List<VcfMetaInfo> metaInfo = new ArrayList<>(20);
        String candidateHeader = null;
        while (it.hasNext() && candidateHeader==null) {
            String line = it.next().trim();
            if (line.startsWith(VcfMetaInfo.PREFIX)) {
                metaInfo.add(new VcfMetaInfo(line));
            }
            else {
                candidateHeader = line;
            }
        }
        checkHeaderLine(candidateHeader, it.file());
        String[] headerFields = StringUtil.getFields(candidateHeader, Const.tab);

        this.file = it.file();
        this.metaInfoLines = metaInfo.toArray(new VcfMetaInfo[0]);
        this.headerLine = candidateHeader;
        this.nHeaderFields = headerFields.length;
        this.includedIndices = includedIndices(headerFields, sampleFilter);
        this.samples = samples(headerFields, includedIndices);
    }

    private static void checkHeaderLine(String line, File file) {
        if (line == null || line.startsWith("#")==false) {
            String s = "Missing line (#CHROM ...) after meta-information lines"
                    + Const.nl + "File source: " + (file==null ? "stdin" : file)
                    + Const.nl + line;
            throw new IllegalArgumentException(s);
        }
        if (line.startsWith(HEADER_PREFIX) == false) {
            if (line.equals(SHORT_HEADER_PREFIX)==false) {
                String s = "Missing header line (file source: "
                        + (file==null ? "stdin" : file) + ")"
                        + Const.nl + "The first line after the initial meta-information lines"
                        + Const.nl + "does not begin with: "
                        + Const.nl + HEADER_PREFIX
                        + Const.nl + line
                        + Const.nl + "The data fields in the header line must be tab-separated.";
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static int[] includedIndices(String[] headerFields,
            Filter<String> sampleFilter) {
        int nUnfilteredSamples = Math.max(headerFields.length - sampleOffset, 0);
        int[] includedIndices = new int[nUnfilteredSamples];
        int index = 0;
        for (int j=0; j<nUnfilteredSamples; ++j) {
            if (sampleFilter.accept(headerFields[sampleOffset + j])) {
                includedIndices[index++] = j;
            }
        }
        if (index < includedIndices.length) {
            includedIndices = Arrays.copyOf(includedIndices, index);
        }
        return includedIndices;
    }

    private Samples samples(String[] headerFields, int[] includedIndices) {
        String[] ids = new String[includedIndices.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = headerFields[sampleOffset + includedIndices[j]];
        }
        return Samples.fromIds(ids);
    }

    /**
     * Returns the file from which data are read, or returns
     * {@code null} if the source is standard input.
     * @return the file from which data are read, or
     * {@code null} if the source is standard input
     */
    public File file() {
        return file;
    }

    /**
     * Returns the number of VCF meta-information lines. VCF meta-information
     * lines are lines that precede the VCF header line. A VCF meta-information
     * line must begin with "##".
     *
     * @return the number of VCF meta-information lines
     */
     public int nMetaInfoLines() {
         return metaInfoLines.length;
     }

    /**
      * Returns the specified VCF meta-information line.

      * @param index a VCF meta-information line index
      * @return the specified VCF meta-information line
      *
      * @throws IndexOutOfBoundsException if
      * {@code index < 0 || index >= this.nMetaInfoLines()}
      */
     public VcfMetaInfo metaInfoLine(int index) {
         return metaInfoLines[index];
     }

     /**
      * Returns the VCF header line.  The VCF header line begins with "#CHROM".
      * @return the VCF header line
      */
     public String headerLine() {
         return headerLine;
     }

     /**
      * Returns the number of fields in the VCF header line before sample
      * exclusions.
      * @return the number of fields in the VCF header line before sample
      * exclusions
      */
     public int nHeaderFields() {
         return nHeaderFields;
     }

     /**
      * Returns the number of samples before sample exclusions.
      * @return the number of samples before sample exclusions
      */
     public int nUnfilteredSamples() {
         return Math.max(0, nHeaderFields - sampleOffset);
     }

     /**
      * Returns the index of the specified sample in the the list original
      * list of samples before sample exclusions.
      * @param sample a sample index
      * @return the index of the specified sample in the the list original
      * list of samples before sample exclusions
      * @throws IndexOutOfBoundsException if
      * {@code sample < 0 || sample >= this.nSamples()}
      */
     public int unfilteredSampleIndex(int sample) {
         return includedIndices[sample];
     }

     /**
      * Returns the number of samples after sample exclusions.
      * @return the number of samples after sample exclusions
      */
     public int nSamples() {
         return samples.nSamples();
     }

    /**
     * Return the list of samples after sample exclusions.
     * @return the list of samples after sample exclusions
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns {@code this.sample().ids()}.
     * @return {@code this.sample().ids()}
     */
    public String[] sampleIds() {
        return samples.ids();
    }

    /**
     * Returns the VCF meta-information lines and the VCF header line used to
     * construct {@code this}.
     * @return the VCF meta-information lines and the VCF header line used to
     * construct {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(400);
        for (int j=0; j<metaInfoLines.length; ++j) {
            sb.append(metaInfoLines[j]);
            sb.append(Const.nl);
        }
        sb.append(headerLine);
        sb.append(Const.nl);
        return sb.toString();
    }
}
