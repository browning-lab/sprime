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
package beagleutil;

import java.util.Arrays;

/**
 * <p>Class {@code Samples} stores a list of samples.
 * </p>
 * Instances of class {@code Samples} are immutable.
 *
 * @author Brian L. Browning
 */
public final class Samples {

    private static final SampleIds sampleIds = SampleIds.instance();
    private final int[] idIndexToIndex;
    private final int[] indexToIdIndex;

    /**
     * Constructs a new instance of {@code Samples} corresponding to
     * the specified list of sample identifier indices.
     * @param idIndices an array of sample identifier indices
     *
     * @throws IllegalArgumentException if the specified array
     * has two or more elements that are equal
     * @throws IndexOutOfBoundsException if any element of the specified
     * array is negative or greater than or equal to
     * {@code beagleutil.SampleIds.size()}
     * @throws NullPointerException if {@code idIndices == null}
     */
    public Samples(int[] idIndices) {
        int[] copy = idIndices.clone();
        this.idIndexToIndex = idIndexToIndex(copy);
        this.indexToIdIndex = copy;
    }

    private static int[] idIndexToIndex(int[] idIndices) {
        int[] idIndexToIndex = new int[sampleIds.size()];
        Arrays.fill(idIndexToIndex, -1);
        for (int j=0; j<idIndices.length; ++j) {
            int idIndex = idIndices[j];
            if (idIndexToIndex[idIndex] != -1) {
                String s = "duplicate sample: " + sampleIds.id(idIndex) +
                        " (ID index: " + idIndex + ")";
                throw new IllegalArgumentException(s);
            }
            else {
                idIndexToIndex[idIndex] = j;
            }
        }
        return idIndexToIndex;
    }

    /**
     * Constructs and returns a {@code Samples} instance
     * corresponding to the specified list of sample identifiers.
     * @param ids an array of sample identifiers.
     * @return a {@code Samples} instance corresponding to the specified
     * list of sample identifiers
     *
     * @throws IllegalArgumentException if the specified array
     * has two or more elements that are equal as strings
     * @throws NullPointerException if {@code ids == null}
     */
    public static Samples fromIds(String[] ids) {
        return new Samples(sampleIds.getIndices(ids));
    }

    /**
     * <p>Returns a hash code value for the object.
     * </p>
     * <p>The hash code is defined by the following calculation:
     * </p>
     * <pre>
        int hash = 1;
        for (int j, n=this.nSamples(); j&lt;n; ++j) {
            hash = 31 * hash + this.idIndex(j);
        }
     * </pre>
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(this.indexToIdIndex);
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}, and returns {@code false}
     * otherwise.
     * @param obj the object to be tested for equality with {@code this}
     * @return {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null || this.getClass() != obj.getClass()) {
            return false;
        }
        final Samples other = (Samples) obj;
        return Arrays.equals(this.indexToIdIndex, other.indexToIdIndex);
    }

    /**
     * Returns the sample identifier index corresponding to the sample
     * with the specified index in this list of samples.
     * @param index a sample index
     * @return the sample identifier index corresponding to the sample
     * with the specified index in this list of samples
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nSamples()}
     */
    public int idIndex(int index) {
        return indexToIdIndex[index];
    }

    /**
     * Returns the index of the sample that corresponds to the
     * specified sample identifier index, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @param idIndex a sample identifier index
     * @return the index of the sample that corresponds to the
     * specified sample identifier index, or returns {@code -1}
     * if there is no corresponding sample in this list of samples
     * @throws IndexOutOfBoundsException if {@code index < 0}
     */
    public int index(int idIndex) {
        if (idIndex >= idIndexToIndex.length) {
            return -1;
        }
        return idIndexToIndex[idIndex];
    }

    /**
     * Returns the index of the sample that corresponds to the
     * specified sample identifier, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @param id a sample identifier
     * @return the index of the sample that corresponds to the
     * specified sample identifier, or returns {@code -1}
     * if there is no corresponding sample in this list of samples
     * @throws NullPointerException if {@code id == null}
     */
    public int index(String id) {
        int idIndex = SampleIds.instance().getIndexIfIndexed(id);
        if (idIndex != -1) {
            return index(idIndex);
        }
        else {
            return -1;
        }
    }

    /**
     * Returns the number of samples in this list.
     * @return the number of samples in this list
     */
    public int nSamples() {
        return indexToIdIndex.length;
    }

    /**
     * Returns the identifier for the sample with the specified
     * index in this list of samples.
     * @param index a sample index
     * @return the identifier for the sample with the specified
     * index in this list of samples
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nSamples()}
     */
    public String id(int index) {
        return sampleIds.id(indexToIdIndex[index]);
    }

    /**
     * Returns this list of samples as an array of sample identifiers.
     * The returned array has length {@code this.nSamples()}, and it
     * satisfies {@code this.ids()[j].equals(this.id(j))} for
     * {@code 0 <= j && j < this.nSamples()}
     * @return this list of samples as an array of sample identifiers
     */
    public String[] ids() {
        return sampleIds.ids(indexToIdIndex);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.ids())}.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(ids());
    }
}
