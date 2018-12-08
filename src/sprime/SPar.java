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

import blbutil.Const;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code SPar} represents the parameters for an SPrime analysis.
 * </p>
 * <p>Instances of class {@code SPar} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SPar {

    private final String[] args;

    // data input/output parameters
    private final File gt;
    private final File outgroup;
    private final File map;
    private final String out;

    private final File excludesamples;
    private final File excludemarkers;
    private final String chrom;

    // algorithm parameters
    private final double maxfreq;
    private final double minscore;
    private final double mu;

    private static final double DEF_MAX_FREQ = 0.01;
    private static final double DEF_MIN_SCORE = 100_000;
    private static final double DEF_MU = 1.2e-8;

    /**
     * Constructs a new {@code SPar} instance from the specified
     * command line arguments.
     * @param args the SPrime command line arguments
     * @throws IllegalArgumentException if a command line argument
     * is incorrectly specified
     * @throws NumberFormatException if a numeric value for a parameter
     * is incorrectly specified
     * @throws NullPointerException if {@code args == null} or if there
     * exists {@code j} such that {@code (0 <= j && j < args.length)} and
     * {@code (args[j] == null)}
     */
    public SPar(String[] args) {
        double DMIN = Double.MIN_VALUE;
        double DMAX = Double.MAX_VALUE;

        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        // data input/output parameters
        gt = Validate.getFile(
                Validate.stringArg("gt", argsMap, true, null, null));
        outgroup = Validate.getFile(
                Validate.stringArg("outgroup", argsMap, true, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, true, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        excludesamples = Validate.getFile(
                Validate.stringArg("excludesamples", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));

        chrom = Validate.stringArg("chrom", argsMap, false, null, null);

        // algorithm parameters

        maxfreq = Validate.doubleArg("maxfreq", argsMap, false, DEF_MAX_FREQ, 0.0f, 1.0f);
        minscore = Validate.doubleArg("minscore", argsMap, false, DEF_MIN_SCORE, DMIN, DMAX);
        mu = Validate.doubleArg("mu", argsMap, false, DEF_MU, DMIN, DMAX);
        Validate.confirmEmptyMap(argsMap);
    }

    /**
     * Returns the SPrime command line arguments.
     * @return the SPrime command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns SPrime usage instructions.
     * @return SPrime usage instructions
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Syntax: " + SMain.COMMAND + " [arguments in format: parameter=value]" + nl + nl

                + "  gt=[VCF file with no missing genotypes]            (required)" + nl
                + "  outgroup=[file with 1 sample ID per line)]         (required)" + nl
                + "  map=[PLINK map file with cM units]                 (required)" + nl
                + "  out=[output file prefix]                           (required)" + nl +nl

                + "  excludesamples=[file with 1 sample ID per line]    (optional)" + nl
                + "  excludemarkers=[file with 1 marker ID per line]    (optional)" + nl
                + "  chrom=[ [chrom] or [chrom]:[start]-[end] ]         (optional)" + nl + nl

                + "  maxfreq=[max outgroup variant frequency)]          (default="
                        + DEF_MAX_FREQ + ")" + nl
                + "  minscore=[min score of an introgressed segment]    (default="
                        + DEF_MIN_SCORE + ")" + nl
                + "  mu=[mutation rate (mutations/bp/meiosis)]          (default="
                        + DEF_MU + ")" + nl;
    }

    // data input/output parameters

    /**
     * Returns the gt parameter.
     * @return the gt parameter
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the outgroup parameter.
     * @return the outgroup parameter
     */
    public File outgroup() {
        return outgroup;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
     */
    public File map() {
        return map;
    }

    /**
     * Returns the chrom parameter or {@code null}
     * if no chrom parameter was specified.
     *
     * @return the chrom parameter or {@code null}
     * if no chrom parameter was specified
     */
    public String chrom() {
        return chrom;
    }

    // general parameters

    /**
     * Returns the maxfreq parameter.
     * @return the maxfreq parameter
     */
    public double maxfreq() {
        return maxfreq;
    }

    /**
     * Returns the maxscore parameter.
     * @return the minscore parameter
     */
    public double minscore() {
        return minscore;
    }

    /**
     * Returns the mu parameter.
     * @return the mu parameter
     */
    public double mu() {
        return mu;
    }
}
