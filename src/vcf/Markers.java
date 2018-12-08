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

import beagleutil.ChromIds;
import blbutil.Const;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>Class {@code Markers} represent a list of markers in chromosome order.
 * </p>
 * <p>Instances of class {@code Markers} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Markers {

    private final Set<Marker> markerSet;

    private final Marker[] fwdMarkerArray;
    private final int[] fwdSumAlleles;
    private final int[] fwdSumGenotypes;
    private final int[] fwdSumHaplotypeBits;
    private final int fwdHashCode;

    private final Marker[] revMarkerArray;
    private final int[] revSumAlleles;
    private final int[] revSumGenotypes;
    private final int[] revSumHaplotypeBits;
    private final int revHashCode;
    private final Markers revMarkers;

    /**
     * Returns a new {@code Markers} instance that is constructed from
     * the specified data.
     * @param markers a list of markers in chromosome order
     * @return a new {@code Markers} instance corresponding to the
     * specified list of markers
     *
     * @throws IllegalArgumentException if markers on a chromosome are not
     * in chromosome order
     * @throws IllegalArgumentException if there are duplicate markers
     * @throws IllegalArgumentException if the markers on a chromosome
     * do not form a contiguous set of entries within the array
     *
     * @throws NullPointerException if
     * {@code markers == null} or if {@code markers[j] == null}
     * for any {@code j} satisfying {@code (0 <= j && j < markers.length)}
     */
    public static Markers create(Marker[] markers) {
        Markers fwd = new Markers(markers);
        return new Markers(fwd.reverse());
    }

    /**
     * Construct a new {@code Markers} instance that represents the
     * specified list of markers.
     * @param markers a list of markers in chromosome order
     *
     * @throws IllegalArgumentException if markers on a chromosome are not
     * in chromosome order
     * @throws IllegalArgumentException if there are duplicate markers
     * @throws IllegalArgumentException if the markers on a chromosome
     * do not form a contiguous set of entries within the array
     *
     * @throws NullPointerException if
     * {@code markers == null} or if {@code markers[j] == null}
     * for any {@code j} satisfying {@code (0 <= j && j < markers.length)}
     */
    private Markers(Marker[] markers) {
        checkMarkerPosOrder(markers);
        this.fwdMarkerArray = markers.clone();
        this.revMarkerArray = reverse(this.fwdMarkerArray);
        this.markerSet = markerSet(fwdMarkerArray);

        this.fwdSumAlleles = cumSumAlleles(fwdMarkerArray);
        this.fwdSumGenotypes = cumSumGenotypes(fwdMarkerArray);
        this.fwdSumHaplotypeBits = cumSumHaplotypeBits(fwdMarkerArray);
        this.fwdHashCode = Arrays.deepHashCode(fwdMarkerArray);

        this.revSumAlleles = cumSumAlleles(revMarkerArray);
        this.revSumGenotypes = cumSumGenotypes(revMarkerArray);
        this.revSumHaplotypeBits = cumSumHaplotypeBits(revMarkerArray);
        this.revHashCode = Arrays.deepHashCode(revMarkerArray);
        this.revMarkers = null;
    }

    /**
     * Constructs a new {@code Markers} instance whose {@code reverse()}
     * method returns the specified {@code Markers}
     * @param revMarkers a list of markers
     */
    private Markers(Markers revMarkers) {
        this.markerSet = revMarkers.markerSet;
        this.fwdMarkerArray = revMarkers.revMarkerArray;
        this.revMarkerArray = revMarkers.fwdMarkerArray;

        this.fwdSumAlleles = revMarkers.revSumAlleles;
        this.fwdSumGenotypes = revMarkers.revSumGenotypes;
        this.fwdSumHaplotypeBits = revMarkers.revSumHaplotypeBits;
        this.fwdHashCode = revMarkers.revHashCode;

        this.revSumAlleles = revMarkers.fwdSumAlleles;
        this.revSumGenotypes = revMarkers.fwdSumGenotypes;
        this.revSumHaplotypeBits = revMarkers.fwdSumHaplotypeBits;
        this.revHashCode = revMarkers.fwdHashCode;
        this.revMarkers = revMarkers;
    }

    private static void checkMarkerPosOrder(Marker[] markers) {
        if (markers.length < 2) {
            return;
        }
        Set<Integer> chromIndices = new HashSet<>();
        chromIndices.add(markers[0].chromIndex());
        chromIndices.add(markers[1].chromIndex());
        for (int j=2; j<markers.length; ++j) {
            int chr0 = markers[j-2].chromIndex();
            int chr1 = markers[j-1].chromIndex();
            int chr2 = markers[j].chromIndex();
            if (chr0 == chr1 && chr1==chr2) {
                int pos0 = markers[j-2].pos();
                int pos1 = markers[j-1].pos();
                int pos2 = markers[j].pos();
                if ( (pos1<pos0 && pos1<pos2) || (pos1>pos0 && pos1>pos2) ) {
                    String s = "markers not in chromosomal order: "
                            + Const.nl + markers[j-2]
                            + Const.nl + markers[j-1]
                            + Const.nl + markers[j];
                    throw new IllegalArgumentException(s);
                }
            }
            else if (chr1!=chr2) {
                if (chromIndices.contains(chr2)) {
                    String s = "markers on chromosome are not contiguous: "
                            + ChromIds.instance().id(chr2);
                    throw new IllegalArgumentException(s);
                }
                chromIndices.add(chr2);
            }
        }
    }

    private static Marker[] reverse(Marker[] markers) {
        int lastIndex = markers.length - 1;
        Marker[] rev = new Marker[markers.length];
        for (int j=0; j<markers.length; ++j) {
            rev[j] = markers[lastIndex - j];
        }
        return rev;
    }

    private static Set<Marker> markerSet(Marker[] markers) {
        Set<Marker> markerSet = new HashSet<>(markers.length);
        for (Marker m : markers) {
            if (markerSet.add(m)==false) {
                throw new IllegalArgumentException("Duplicate marker: " + m);
            }
        }
        return markerSet;
    }

    private static int[] cumSumAlleles(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + markers[j-1].nAlleles();
        }
        return ia;
    }

    private static int[] cumSumGenotypes(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + markers[j-1].nGenotypes();
        }
        return ia;
    }

    private static int[] cumSumHaplotypeBits(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            int nAllelesM1 = markers[j-1].nAlleles() - 1;
            int nStorageBits = Integer.SIZE
                    - Integer.numberOfLeadingZeros(nAllelesM1);
            ia[j] = ia[j-1] + nStorageBits;
        }
        return ia;
    }

    /**
     * Returns a hash code value for the object.
     * The returned hash code equals
     * {@code Arrays.deepHashCode(this.markers())}.
     * @return a hash code value for the object
     */
    @Override
    public int hashCode() {
        return fwdHashCode;
    }

    /**
     * Returns {@code true} if the specified object is a {@code Markers}
     * instance which represents the same list of markers as {@code this},
     * and returns {@code false} otherwise. Two lists of markers are
     * the same if the lists have the same size and if markers with the
     * same index in the two lists are equal.
     *
     * @param obj the object to be tested for equality with {@code this}
     *
     * @return {@code true} if the specified object is a {@code Markers}
     * instance which represents the same list of markers as {@code this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Markers other = (Markers) obj;
        return Arrays.deepEquals(this.fwdMarkerArray, other.fwdMarkerArray);
    }

    /**
     * Returns a {@code Markers} instance that is obtained by reversing the
     * order of markers in {@code this}.
     * @return a new {@code Markers} instance that is obtained by reversing
     * the order of markers in {@code this}
     */
    public Markers reverse() {
        return revMarkers==null ? new Markers(this) : revMarkers;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return fwdMarkerArray.length;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public Marker marker(int marker) {
        return fwdMarkerArray[marker];
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Marker[] markers() {
        return fwdMarkerArray.clone();
    }

    /**
     * Returns {@code true} if the specified marker is not {@code null}
     * and is an element in the list of markers represented by {@code this},
     * and returns {@code false} otherwise.
     *
     * @param marker a marker
     *
     * @return {@code true} if the specified marker is not {@code null} and
     * is an element in the list of markers represented by {@code this}
     */
    public boolean contains(Marker marker) {
        return markerSet.contains(marker);
    }

    /**
     * Returns a {@code Markers} instance that represents
     * the specified range of marker indices.
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @return a {@code Markers} instance that represents
     * the specified range of marker indices
     *
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > this.nMarkers()}
     * @throws IllegalArgumentException if {@code start >= end}.
     */
    public Markers restrict(int start, int end) {
        if (end > fwdMarkerArray.length) {
            throw new IndexOutOfBoundsException("end > this.nMarkers(): " + end);
        }
        return new Markers(Arrays.copyOfRange(fwdMarkerArray, start, end));
    }

    /**
     * Returns the sum of the number of alleles for
     * the markers with index less than the specified index.
     * @param marker a marker index
     * @return the sum of the number of alleles for
     * the markers with index less than the specified index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker > this.nMarkers()}
     */
    public int sumAlleles(int marker) {
        return fwdSumAlleles[marker];
    }

    /**
     * Returns {@code this.sumAlleles(this.nMarkers())}.
     * @return {@code this.sumAlleles(this.nMarkers())}
     */
    public int sumAlleles() {
        return fwdSumAlleles[fwdMarkerArray.length];
    }

    /**
     * Returns the sum of the number of possible genotypes for the markers
     * with index less than the specified index.
     * @param marker a marker index
     * @return the sum of the number of possible genotypes for the markers
     * with index less than the specified index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker > this.nMarkers()}
     */
    public int sumGenotypes(int marker) {
        return fwdSumGenotypes[marker];
    }

    /**
     * Returns {@code this.sumGenotypes(this.nMarkers())}.
     * @return {@code this.sumGenotypes(this.nMarkers())}
     */
    public int sumGenotypes() {
        return fwdSumGenotypes[fwdMarkerArray.length];
    }

    /**
     * Returns the number of bits requires to store a haplotype for the
     * markers with index less than the specified index.
     * @param marker a marker index
     * @return the number of bits requires to store a haplotype for the
     * markers with index less than the specified index
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker > this.nMarkers()}
     */
    public int sumHaplotypeBits(int marker) {
        return fwdSumHaplotypeBits[marker];
    }

    /**
     * Returns {@code this.sumHaplotypeBits(this.nMarkers())}.
     * @return {@code this.sumHaplotypeBits(this.nMarkers())}
     */
    public int sumHaplotypeBits() {
        return fwdSumHaplotypeBits[fwdMarkerArray.length];
    }

    /**
     * Returns a string representation of {@code this}.
     * The exact details of the representation are unspecified and
     * subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(fwdMarkerArray);
    }
}
