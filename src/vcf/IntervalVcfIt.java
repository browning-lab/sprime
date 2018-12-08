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

import blbutil.SampleFileIt;
import beagleutil.ChromInterval;
import beagleutil.Samples;
import blbutil.Const;
import java.io.File;
import java.util.NoSuchElementException;

/**
 * <p>Class {@code IntervalVcfIterator} is a sample file iterator whose
 * {@code next()} method returns a marker container.
 * </p>
 *
 * @param <E> the type parameter
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IntervalVcfIt<E extends MarkerContainer>
        implements SampleFileIt<E> {

    private final SampleFileIt<E> it;
    private final ChromInterval interval;
    private E next;

    /**
     * Constructs a new {@code IntervalVcfIterator} instance.
     * @param it an iterator whose {@code next()} method returns a marker
     * container
     * @param interval a chromosome interval
     * @throws NullPointerException if {@code it == null || interval == null}
     */
    public IntervalVcfIt(SampleFileIt<E> it, ChromInterval interval) {
        E firstRecord = readFirstRecord(it, interval);
        if (firstRecord==null) {
            String s = "No VCF records found in specified interval."
                    + Const.nl + "Check chromosome identifier ["
                    + interval.chrom() + "] and interval ["
                    + interval.start() + "-" + interval.end() + "]";
            throw new IllegalArgumentException(s);
        }
        this.it = it;
        this.interval = interval;
        this.next = firstRecord;
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements.
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        E current = next;
        this.next = readNextRecord(it, interval);
        return current;
    }

    private E readFirstRecord(SampleFileIt<E> it, ChromInterval interval) {
        E nextRecord = null;
        while (nextRecord==null && it.hasNext()) {
            E candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private E readNextRecord(SampleFileIt<E> it, ChromInterval interval) {
        E nextRecord = null;
        if (it.hasNext()) {
            E candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private static boolean inInterval(ChromInterval interval, Marker marker) {
        return (marker.chromIndex() == interval.chromIndex()
                && interval.start() <= marker.pos()
                && marker.pos() <= interval.end());
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public void close() {
        it.close();
    }
}
