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
 * <p>Class {@code IntSet} represents an indexed set of integers.
 * </p>
 * <p>Class {@code IntSet} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IntSet {

    private static final int NIL = -1;
    private static final float LOAD_FACTOR = 0.75f;

    private int size;
    private int nBuckets;

    private int[] next;
    private int[] indices; // stores list index of element
    private int[] list;
    private int firstFreeIndex;

    /**
     * Creates a new {@code IntSet} instance.
     *
     * @param capacity the initial capacity of the set
     * @throws IllegalArgumentException if
     * {@code capacity < 0 || (capacity > (1 << 30))}
     */
    public IntSet(int capacity) {
        if (capacity < 0 || capacity > (1<<30)) {
            throw new IllegalArgumentException(String.valueOf(capacity));
        }
        int numBuckets = (int) Math.ceil((capacity+1)/LOAD_FACTOR);

        assert numBuckets>0;
        allocateArrays(capacity, numBuckets);
        initializeFields(numBuckets);
    }

    private void allocateArrays(int capacity, int buckets) {
        this.next = new int[buckets + capacity];
        this.indices = new int[buckets + capacity];
        this.list = new int[capacity];
    }

    private void initializeFields(int numBuckets) {
        size = 0;
        nBuckets = numBuckets;
        firstFreeIndex = nBuckets;
        Arrays.fill(next, 0, nBuckets, NIL);
        for (int j=nBuckets; j<next.length; ++j) {
            next[j] = j+1;
        }
    }

    /*
     * Increases the capacity of the internal hash table.
     */
    private void rehash(int newCapacity) {
        if (newCapacity > size) {
            int oldSize = size;
            int[] oldList = list;
            int newNumBuckets = (int) Math.ceil(newCapacity/LOAD_FACTOR);
            allocateArrays(newCapacity, newNumBuckets);
            initializeFields(newNumBuckets);
            for (int j=0; j<oldSize; ++j) {
                add(oldList[j]);
            }
        }
    }

    /**
     * Returns {@code true} if the set contains the specified element,
     * and returns {@code false} otherwise.
     * @param element an nonnegative integer
     * @return {@code true} if the set contains the specified element
     */
    public boolean contains(int element) {
        int index = next[bucket(element)];
        while (index!=NIL && list[indices[index]]<element) {
            index = next[index];
        }
        return (index!=NIL && list[indices[index]]==element);
    }

    /**
     * Adds the specified element to this set.  The indexing of set elements
     * immediately before and after this command is invoked may differ if
     * the set is changed by the operation.
     * @param element an integer to add to this set
     * @return {@code true} if the set was changed by the operation, and
     * {@code false} otherwise
     */
    public boolean add(int element) {
        int prevIndex = prevIndex(element);
        int nextIndex = next[prevIndex];
        if (nextIndex==NIL || list[indices[nextIndex]]!=element) {
            int index = firstFreeIndex;
            firstFreeIndex = next[firstFreeIndex];
            next[prevIndex] = index;
            indices[index] = size;
            next[index] = nextIndex;
            list[size++] = element;
            if (size == list.length) {
                int newCapacity = 3*list.length/2 + 1;
                rehash(newCapacity);
            }
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Removes the specified element from this set. The indexing of set elements
     * immediately before and after this command is invoked may differ if
     * the set is changed by the operation.
     *
     * @param element an integer to remove this set
     * @return {@code true} if the set was changed by the operation, and
     * {@code false} otherwise
     */
    public boolean remove(int element) {
        int prevIndex = prevIndex(element);
        int index = next[prevIndex];
        if (index==NIL || list[indices[index]]!=element) {
            return false;
        }
        else {
            int oldListIndex = indices[index];
            next[prevIndex] = next[index];
            next[index] = firstFreeIndex;
            firstFreeIndex = index;

            --size;
            if (oldListIndex!=size) {
                index = index(list[size]);
                indices[index] = oldListIndex;
                list[oldListIndex] = list[size];   // overwrite removed element
            }
            return true;
        }
    }

    private int bucket(int element) {
        return Math.abs((71*element) % nBuckets);
    }

    private int prevIndex(int element) {
        int prevIndex = bucket(element);
        int index = next[prevIndex];
        while (index!=NIL && list[indices[index]]<element) {
            prevIndex = index;
            index = next[index];
        }
        return prevIndex;
    }

    private int index(int element) {
        int index = next[bucket(element)];
        while (index!=NIL && list[indices[index]]<element) {
            index = next[index];
        }
        return (index!=NIL && list[indices[index]]==element) ? index : NIL;
    }

    /**
     * Returns the specified element.
     * @param index an index of an element in this set
     * @return the specified element
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int elementWithIndex(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return list[index];
    }

    /**
     * Removes all elements from this set.
     */
    public void clear() {
        initializeFields(nBuckets);
    }

    /**
     * Returns the number of elements in this set.
     *
     * @return the number of elements in this set
     */
    public int size() {
        return size;
    }

    /**
     * Returns the capacity of this set.  The capacity of this set
     * is the maximum number of elements that may be stored without
     * allocating more memory.
     *
     * @return the capacity of this set
     */
    public int capacity() {
        return list.length;
    }

    /**
     * Sets the capacity of this list to the specified value. The capacity
     * of this set is the maximum number of elements that may be stored
     * without allocating more memory.
     * @param capacity the desired capacity
     * @throws IllegalArgumentException if {@code capacity < this.size()}
     */
    public void setCapacity(int capacity) {
        if (capacity < size) {
            throw new IllegalArgumentException(String.valueOf(capacity));
        }
        if (capacity != list.length) {
            rehash(capacity);
        }
    }

    /**
     * Returns an array containing the elements in this set. The returned
     * array will satisfy:
     * {@code this.toArray()[j]==this.elementWithIndex(j)} for each
     * {@code j} satisfying {@code 0 < j && j < this.size()}
     * @return an array containing the elements in this set
     */
    public int[] toArray() {
        return Arrays.copyOf(list, size);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.toArray())}.
     *
     * @return {@code java.util.Arrays.toString(this.toArray())}
     */
    @Override
    public String toString() {
        return Arrays.toString(toArray());
    }
}
