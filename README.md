The **SPrime** program identifies variants that are introgressed
from archaic populations.

If you use **SPrime** in a published analysis, please report the program
version printed in the first line of the output **log** file, and please
cite the article that describes the SPrime method:

> S R Browning, B L Browning, Y Zhou, S Tucci, J M Akey (2018).
> Analysis of human sequence data reveals two pulses of archaic
> Denisovan admixture. Cell 173(1):53-61.
> doi: [10.1016/j.cell.2018.02.031](http://dx.doi.org/10.1016/j.cell.2018.02.031)

Sharon Browning and Brian Browning

Last updated Dec 4, 2018

## Installation

You can download the latest executable file,
[sprime.jar](https://faculty.washington.edu/browning/sprime.jar) with the
command:

    wget https://faculty.washington.edu/browning/sprime.jar

or you can download the source files and create the executable file
with the commands:

    git clone https://github.com/browning-lab/sprime.git
    javac -cp sprime/src/ sprime/src/sprime/SMain.java
    jar cfe sprime.jar sprime/SMain -C sprime/src/ ./
    jar -i sprime.jar

## Running SPrime

The SPrime program requires Java version 1.8 (or a later version). Use of an
earlier Java version will produce an "Unsupported Class Version" error.

The command:

    java -jar sprime.jar

prints a summary of the command line arguments.

To run SPrime, enter the following command:

    java -Xmx[GB]g -jar sprime.jar [arguments]

where **[GB]** is the maximum number of gigabytes of memory to use, and
**[arguments]** is a space-separated list of parameter values, each expressed as
**parameter=value**.

The shell script
[run.sprime.test](https://raw.githubusercontent.com/browning-lab/sprime/master/test/run.sprime.test)
will run a test sprime analysis.

### Required Parameters

SPrime has four required parameters:

* **gt=[file]** where **[file]** is a
[Variant Call Format](https://faculty.washington.edu/browning/intro-to-vcf.html)
(VCF) file with genotypes for **all** autosomes.  All VCF records must
include a GT FORMAT subfield, and all genotypes must have non-missing
alleles.  A VCF record may have multiple ALT alleles.
Files with name ending in ".gz" are assumed to be gzip-compressed.
* **outgroup=[file]** where **[file]** is a text file with one sample
identifier per line. The outgroup file identifies the samples from an
outgroup population that is not expected to contain introgressed variants.
Outgroup samples must have genotype data in the VCF file specified with the
**gt** parameter.
* **map=[file]** where **[file]** is a
[PLINK format genetic map](http://zzz.bwh.harvard.edu/plink/data.shtml#map)
with cM units. SPrime will use linear interpolation to estimate genetic
positions between map positions. Chromosome identifiers in the genetic map and
input VCF file must match. HapMap genetic maps in cM units are available for
[GRCh36](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh36.map.zip),
[GRCh37](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip), and
[GRCh38](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip).

* **out=[string]** where **[string]** is a filename prefix for the output files.

### Optional Parameters
* **chrom=[chrom]:[start]‑[end]** limits the output to a chromosome interval
where **[chrom]** is the CHROM identifier in the input VCF file, and 
**[start]** and **[end]** are the start and end positions. An entire chromosome, 
the beginning of a chromosome, or the end of a chromosome may be specified 
with **chrom=[chrom]**, **chrom=[chrom]:‑[end]**, and 
**chrom=[chrom]:[start]‑** respectively. If the **chrom** parameter is used, 
the input VCF file must still contain **all** autosomes so that the
global variant density can be accurately estimated.
* **maxfreq=[nonnegative number ≤ 1.0]** specifies the maximum frequency of an
introgressed variant in the outgroup (default: **maxfreq=0.01**),  Variants with
outgroup frequency less than or equal to **maxfreq** are ignored.
* **minscore=[nonnegative number]** specifies the minimum score of an introgressed
segment (default: **minscore=1.0e5**).
* **mu=[positive number]** specifies the genomewide mutation rate per base pair
per meiosis (default: **mu=1.2e-8**).
* **excludesamples=[file]** specifies a file containing samples to be excluded
from the analysis (one sample identifier per line).
* **excludemarkers=[file]** specifies a file containing markers to be excluded
from the analysis (one marker per line). Each line of the file can be a
marker identifier from a VCF record’s ID field or a genomic coordinate expressed
as CHROM:POS.

## Output files
SPrime produces two output files:

The **log** file (.log) contains a summary of the analysis.

The **score** file (.score) lists all introgressed variants with segment
score greater than or equal to the **minscore** parameter. The score file has
eight tab-delimited columns. The first five columns are the variant's CHROM,
POS, ID, REF, and ALT fields from the input VCF file. The last three columns are
the segment index (SEGMENT), the allele (ALLELE), and the segment score (SCORE).

## License
The SPrime program is licensed under the Apache License, Version 2.0 (the License).
You may obtain a copy of the License from http://www.apache.org/licenses/LICENSE-2.0

Source files in the net/sf/samtools/ directory are from the Broad Institute
and are licensed under the [MIT License](https://opensource.org/licenses/MIT).
These Broad Institute files are used to perform BGZIP decompression.

