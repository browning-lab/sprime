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
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.Utilities;
import java.io.File;
import java.io.PrintWriter;
import java.util.Locale;

/**
 * <p>Class {@code SMain} contains the {@code static void main(String[] args)}
 * entry method for the SPrime program. See {@code SPar.usage()} for usage
 * instructions.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SMain {

    private static final String PROGRAM = "sprime.jar";
    private static final String VERSION = "__REV__";

    /**
     * The java command to run this version of the sprime program.
     */
    static final String COMMAND = "java -jar sprime.jar";

    /**
     * The copyright string.
     */
    private static final String COPYRIGHT = "Copyright (C) 2017 Brian L. Browning";

    /**
     * A brief help message.
     */
    private static final String SHORT_HELP = PROGRAM + "(version: " + VERSION + ")"
            + Const.nl + SMain.COPYRIGHT
            + Const.nl + "Enter \"" + COMMAND + "\" to print a list of"
                    + "command line arguments";

    private SMain() {
        // private constructor to prevent instantiation
    }

    /**
     * Entry point to the SPrime program.  See {@code SPar.usage()} for
     * usage instructions. The SPrime program will exit with an error message
     * if the command line arguments are incorrectly specified.
     *
     * @param args command line arguments
     *
     * @throws NullPointerException if {@code args == null} or if there
     * exists {@code j} such that {@code (0 <= j && j < args.length)} and
     * {@code (args[j] == null)}
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        if (args.length==0) {
            System.out.print(PROGRAM);
            System.out.print(" (version: ");
            System.out.print(VERSION);
            System.out.println(")");
            System.out.println();
            System.out.println(SPar.usage());
            System.exit(0);
        }
        SPar par = new SPar(args);
        checkOutputPrefix(par);
        checkChromParameter(par);

        try (PrintWriter log = FileUtil.printWriter(new File(par.out() + ".log"))) {
            long t0 = System.nanoTime();
            printStartInfo(log, par);

            SAnalyzer sAnalyzer = new SAnalyzer(par);
            checkNumChroms(log, sAnalyzer);
            sAnalyzer.run();

            printEndInfo(log, sAnalyzer, System.nanoTime() - t0);
        }
    }

    private static void checkOutputPrefix(SPar par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(SPar.usage() + s);
        }

        File out = new File(par.out() + ".seg");
        if (out.equals(par.gt())
                || out.equals(par.map())
                || out.equals(par.outgroup())
                || out.equals(par.excludemarkers())
                || out.equals(par.excludesamples())) {
            String s = "ERROR: Output file equals input file: " + out;
            Utilities.exit(SPar.usage() + s);
        }
    }

    private static void checkChromParameter(SPar par) {
        if (par.chrom()!=null && ChromInterval.parse(par.chrom())==null) {
            String s = SHORT_HELP + Const.nl
                    + Const.nl + "ERROR: invalid \"chrom\" parameter: \""
                    + par.chrom() + "\""
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
    }

    private static void checkNumChroms(PrintWriter log, SAnalyzer sAnalyzer) {
        int nChrom = sAnalyzer.nChrom();
        if (nChrom==1) {
            Utilities.duoPrint(log,
                      Const.nl + "WARNING: The input VCF file contains only one chromosome."
                    + Const.nl + "All autosomes must be included in the input VCF file in"
                    + Const.nl + "order for SPrime to estimate a genomewide variant density."
                    + Const.nl);
        }
    }

    private static void printStartInfo(PrintWriter log, SPar par) {
        Utilities.duoPrint(log, SMain.SHORT_HELP + Const.nl);
        Utilities.duoPrintln(log, Const.nl + "Start time: " + Utilities.timeStamp());
        Utilities.duoPrint(log, Utilities.commandLine(PROGRAM, par.args()));
    }

    private static void printEndInfo(PrintWriter log, SAnalyzer sAnalyzer,
            long nanos) {
        Utilities.duoPrint(log, Const.nl + "Number of outgroup samples:      " +
                sAnalyzer.nOutgroupSamples());
        Utilities.duoPrint(log, Const.nl + "Number of target samples:        " +
                sAnalyzer.nTargetSamples());
        Utilities.duoPrint(log, Const.nl + "Variant analyzed:                " +
                sAnalyzer.variantCnt());
        Utilities.duoPrint(log, Const.nl + "Segments detected:               " +
                sAnalyzer.segmentCnt());
        Utilities.duoPrint(log, Const.nl + "Run time:                        ");
        Utilities.duoPrintln(log, Utilities.elapsedNanos(nanos));
        Utilities.duoPrintln(log, Const.nl + "End time: "
                + Utilities.timeStamp());
        Utilities.duoPrintln(log, Const.nl + SMain.PROGRAM + " finished");
    }
}
