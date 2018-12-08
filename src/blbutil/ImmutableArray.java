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
package blbutil;

import java.util.Arrays;

/**
 * <p>Class {@code ImmutableArray} represents an immutable array of
 * immutable objects.
 * </p>
 *
 * @param <E> the type of element in this {@code ImmutableArray}
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImmutableArray<E> {

    private final E[] oa;

    /**
     * Constructs a new {@code ImmutableArray} instance by cloning the
     * specified array.  The behavior of this class is undefined
     * if the elements of the specified array are not immutable.
     * @param oa the array of elements
     * @throws NullPointerException if {@code oa == null}
     */
    public ImmutableArray(E[] oa) {
        this.oa = oa.clone();
    }

    /**
     * Returns the number of elements in this {@code ImmutableArray}.
     * @return the number of elements in this {@code ImmutableArray}
     */
    public int size() {
        return oa.length;
    }

    /**
     * Returns the specified array element.
     * @param index an array index
     * @return the specified array element
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    public E get(int index) {
        return oa[index];
    }

    /**
     * Returns a string representation of this {@code ImmutableArray} by applying
     * {@code java.utils.Arrays.toString()} to an equivalent {@code Object[]}
     * instance.
     *
     * @return a string representation of this {@code ImmutableArray}
     */
    @Override
    public String toString() {
        return Arrays.toString(oa);
    }
}
