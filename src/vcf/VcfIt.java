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
import beagleutil.Samples;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * <p>Class {@code VcfIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record.
 * </p>
 * <p>Instances of class {@code VcfIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 * @param <E> the type parameter
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfIt<E extends MarkerContainer> implements SampleFileIt<E> {

    private static final float DEFAULT_MAX_LR = Float.MAX_VALUE;

    private final VcfHeader vcfHeader;
    private final FileIt<String> it;
    private final Function<String, E> mapper;
    private final Filter<Marker> markerFilter;
    private final Thread fileReaderThread;
    private volatile boolean stopFileReadingThread = false;

    private final BlockingQueue<String[]> stringBuffers;
    private final Deque<E> emBuffer;

    /**
     * The default number of VCF records stored in a buffer, which is 1000.
     */
    public static final int MAX_EM_BUFFER_SIZE = 1000;

    /**
     * A function mapping a string VCF record with GT format fields
     * to a {@code GTRec} object.
     */
    public static final BiFunction<VcfHeader, String, GTRec> toBitSetGT
            = (VcfHeader h, String s) -> new BitSetGT(h, s);

    /**
     * A function mapping a string VCF record with GL format fields
     * to a {@code VcfRecord} object.
     */
    public static final BiFunction<VcfHeader, String, GTRec> toGLRec
            = (VcfHeader h, String s) -> VcfRecord.fromGL(h, s, DEFAULT_MAX_LR);

    /**
     * A function mapping a string VCF record with GT or GL format fields
     * to a {@code VcfRecord} object.
     */
    public static final BiFunction<VcfHeader, String, GTRec> toGTGLRec
            = (VcfHeader h, String s) -> VcfRecord.fromGTGL(h, s, DEFAULT_MAX_LR);

    /**
     * A function mapping a string VCF record with GT or GL format fields
     * to a {@code VcfRecord} object.
     */
    public static final BiFunction<VcfHeader, String, VcfRecord> toVcfRecord
            = (VcfHeader h, String s) -> VcfRecord.fromGTGL(h, s, DEFAULT_MAX_LR);

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, Filter.acceptAllFilter(), recMapper);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, sampleFilter, Filter.acceptAllFilter(),
                recMapper);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, sampleFilter, markerFilter, recMapper,
                MAX_EM_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @param bufferSize the buffer size
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper, int bufferSize) {
        VcfIt<R> vcfIt = new VcfIt<>(strIt, sampleFilter, markerFilter,
                recMapper, bufferSize);
        vcfIt.start();
        return vcfIt;
    }

    private VcfIt(FileIt<String> it, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, E> recMapper, int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        if (markerFilter==null) {
            markerFilter = Filter.acceptAllFilter();
        }
        this.vcfHeader = new VcfHeader(it, sampleFilter);
        this.it = it;
        this.mapper = (String s) -> recMapper.apply(vcfHeader, s);
        this.markerFilter = markerFilter;
        this.stringBuffers = new ArrayBlockingQueue<>(1);
        this.emBuffer = new ArrayDeque<>(bufferSize);
        this.fileReaderThread = fileReadingThread();
    }

    private void start() {
        this.fileReaderThread.setDaemon(true);
        this.fileReaderThread.start();
        fillEmissionBuffer();
        if (emBuffer.isEmpty()) {
            noRecordFoundError(it);
        }
    }

    private void noRecordFoundError(FileIt<String> it) {
        if (it.hasNext()==false) {
            StringBuilder sb = new StringBuilder(100);
            sb.append("No VCF records found (data source: ");
            sb.append(it.file()==null ? "stdin" : it.file());
            sb.append(")");
            sb.append(Const.nl);
            sb.append("Check that the chromosome identifiers are the same in each input VCF");
            sb.append(Const.nl);
            sb.append("file and in the \'chrom=\' command line argument (if \'chrom=\' is used).");
            throw new IllegalArgumentException(sb.toString());
        }
    }

    private Thread fileReadingThread() {
        Runnable runnable = () -> {
            String line = readLine(it);
            int bufferSize = stringBufferSize(line);
            while (line != null && stopFileReadingThread == false) {
                String chromPlusTab = chromFieldPlusTab(line);
                String[] sa = new String[bufferSize];
                int size = 0;
                while (line != null && size < bufferSize
                        && line.startsWith(chromPlusTab)) {
                    sa[size++] = line;
                    line = readLine(it);
                }
                if (size < bufferSize) {
                    sa = Arrays.copyOf(sa, size);
                }
                putInBlockingQueue(stringBuffers, sa);
            }
            if (stopFileReadingThread == false) {
                putInBlockingQueue(stringBuffers, new String[0]);    // sentinel
            }
        };
        return new Thread(runnable);
    }

    private static int stringBufferSize(String line) {
        if (line == null) {
            return 0;
        }
        long nBytesPerLine = 2*line.length();
        Runtime rt = Runtime.getRuntime();
        long maxMem = rt.maxMemory();
        if (maxMem == Long.MAX_VALUE) {
            maxMem = 500 * (1 << 30);
        }
        long bufferSize = 1 + (maxMem / (20*nBytesPerLine));
        if (bufferSize > MAX_EM_BUFFER_SIZE) {
            bufferSize = MAX_EM_BUFFER_SIZE;
        }
        return (int) bufferSize;
    }

    private static <E> void putInBlockingQueue(BlockingQueue<E> q, E e) {
        try {
            q.put(e);
        } catch (InterruptedException ex) {
            Utilities.exit("ERROR: InterruptedException", ex);
        }
    }

    private static <E> E takeFromBlockingQueue(BlockingQueue<E> q) {
        try {
            return q.take();
        } catch (InterruptedException ex) {
            Utilities.exit("ERROR: InterruptedException", ex);
        }
        assert false;
        return null;
    }

    private static String chromFieldPlusTab(String vcfRecord) {
        int tabIndex = vcfRecord.indexOf(Const.tab);
        if (tabIndex == -1) {
            String s = Const.nl + "ERROR: Missing tab delimiter in VCV Record:"
                    + Const.nl + vcfRecord
                    + Const.nl + "Exiting Program";
            Utilities.exit(s);
        }
        return vcfRecord.substring(0, tabIndex + 1);
    }

    private void fillEmissionBuffer() {
        assert emBuffer.isEmpty();
        int lastLength = -1;
        while (lastLength != 0 && emBuffer.size() < MAX_EM_BUFFER_SIZE) {
            String[] stringBuffer = takeFromBlockingQueue(stringBuffers);
            lastLength = stringBuffer.length;
            if (stringBuffer.length>0) {
                List<E> list = Arrays.stream(stringBuffer)
                        .parallel()
                        .map(mapper)
                        .filter(e -> markerFilter.accept(e.marker()))
                        .collect(Collectors.toList());
                emBuffer.addAll(list);
            }
            else {
                // put sentinel element back
                putInBlockingQueue(stringBuffers, stringBuffer);
            }
        }
    }

    private static String readLine(FileIt<String> it) {
        if (it.hasNext()==false) {
            return null;
        }
        String line = it.next();
        while (line.trim().isEmpty() && it.hasNext()) {
            line = it.next();
        }
        return line;
    }

    @Override
    public void close() {
        stopFileReadingThread = true;
        stringBuffers.poll();  // unblock file reading thread
        try {
            fileReaderThread.join();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        it.close();
        emBuffer.clear();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !emBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        E first = emBuffer.removeFirst();
        if (emBuffer.isEmpty()) {
            fillEmissionBuffer();
        }
        return first;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        String s = "remove() is not supported by VcfIterator";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(it.file()==null ? "stdin" : it.file().toString());
        return sb.toString();
    }
}
