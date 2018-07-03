package cox.external;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import cox.utils.ProcessUtils;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Bowtie2 mapper
 * 
 * @author <a href="http://www.unifi.it/dblage/CMpro-v-p-65.html">Giovanni
 *         Bacci</a>
 *
 */
public class BowTie2 implements Mapper, ExternalApplication {

	/**
	 * The logger
	 */
	private Logger log = Logger.getLogger(this.getClass().getName());

	/**
	 * Verbose property
	 */
	private boolean verbose = false;

	/**
	 * The reference database (once built)
	 */
	private String reference = null;
	
	/**
	 * If <code>true</code> the reference database will be rebuilt even if it exists
	 */
	private boolean forceRebuild = false;

	/**
	 * The name of the bowtie2 log file
	 */
	private String logfile = "bowtie2.log";

	private String outputDir = ".";

	@Override
	public boolean isVerbose() {
		return verbose;
	}

	@Override
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	@Override
	public void setOutputDir(String outputDir) {
		this.outputDir = outputDir;

	}

	@Override
	public String getOutputDir() {
		return outputDir;
	}

	@Override
	public String mapPairedEnd(String forward, String reverse, String... otherOpts) {
		// Guard block
		if (reference == null) {
			log.severe("Reference database has not been instantiated");
			return null;
		}

		String opt = Arrays.stream(otherOpts).collect(Collectors.joining(" "));
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "paired_", ".sam").toString();
			String cmd = String.format("bowtie2 %s -1 %s -2 %s -x %s -S %s", opt, forward, reverse, reference, output);
			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logfile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s building bowtie2 db");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot map reads with bowtie2, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Reads mapping has been interrupted");
		}
		return null;
	}

	@Override
	public String mapSingleEnd(String single, String... otherOpts) {
		// Guard block
		if (reference == null) {
			log.severe("Reference database has not been instantiated");
			return null;
		}

		String opt = Arrays.stream(otherOpts).collect(Collectors.joining(" "));
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "single_", ".sam").toString();
			String cmd = String.format("bowtie2 %s -U %s -x %s -S %s", opt, single, reference, output);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logfile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s building bowtie2 db");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot map reads with bowtie2, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Reads mapping has been interrupted");
		}
		return null;
	}

	@Override
	public boolean buildReference(String reference) {

		// Check if reference database exists. If forceRebuild is set to true
		// the reference will be rebuilt in any case
		if (Files.exists(Paths.get(reference + ".1.bt2")) && !forceRebuild) {
			this.reference = reference;
			return true;
		}
		String cmd = String.format("bowtie2-build %1$s %1$s", reference);

		// Building process
		ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

		try {
			int res = ProcessUtils.runProcess(pb, true, verbose, logfile);
			if (res != 0) {
				throw new IOException("Error/s building bowtie2 db, see log file for details");
			} else {
				this.reference = reference;
				return true;
			}
		} catch (IOException e) {
			log.severe("Cannot build Bowtie2 db, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Bowtie2 build has been interrupted");
		}
		return false;
	}

	/**
	 * Filter a sam or bam file based on mapping quality.
	 * 
	 * @param sam_or_bam
	 *            the input file
	 * @param quality_cutoff
	 *            the quality cutoff
	 * @return the name of the filtered file
	 */
	public String filterLowQualityRecords(String sam_or_bam, int quality_cutoff) {
		// Making reader and writer
		SamReader reader = SamReaderFactory.makeDefault().open(new File(sam_or_bam));
		SAMFileWriter writer = null;

		String extension = (reader.type() == SamReader.Type.BAM_TYPE) ? ".bam" : ".sam";

		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "quality_", extension).toString();
			writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(output));

			reader.iterator().stream().filter(r -> r.getMappingQuality() >= quality_cutoff)
					.forEach(writer::addAlignment);
			return output;
		} catch (IOException e) {
			log.severe("Cannot filter record based on mapping quality, " + e.getMessage());
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (writer != null)
				writer.close();

		}

		return null;
	}

	/**
	 * Produces a sam/bam file reporting only redundant alignment. An alignment
	 * is considered redundant, according to the bowtie2 manual, if:<br>
	 * 1 - It has the AS flag (basically the read aligns somewhere)<br>
	 * 2 - It has the XS flag (basically the read has (at least) another valid
	 * alignment)<br>
	 * 3 - The AS flag is equal to the XS flag and they are both equal to the
	 * best AS flag of the batch (basically: AS == XS == best(AS && XS))<br>
	 * <br>
	 * If the read has the YS flag (basically if it is paired), then it has to
	 * satisfy another condition:<br>
	 * 1 - The YS flag has to be equal to the best YS flag of the batch (the
	 * alignment of the pair has to be equally good)<br>
	 * <br>
	 * This method will collect each alignment based on its ID and then it will
	 * parse all batches using the {@link BowTie2#parseBatch(List)} method
	 * 
	 * @param sam_or_bam
	 *            the sam/bam file
	 * @return a sam/bam file reporting only redundant alignments
	 */
	public String getRedundantAlignments(String sam_or_bam) {

		// Making reader and writer
		SamReader reader = SamReaderFactory.makeDefault().open(new File(sam_or_bam));
		SAMFileWriter writer = null;

		String extension = (reader.type() == SamReader.Type.BAM_TYPE) ? ".bam" : ".sam";

		Iterable<SAMRecord> filtered = reader.iterator().stream().filter(r -> r.getAttribute("AS") != null)::iterator;
		String output = null;

		try {
			// Create output file and setting writer
			output = Files.createTempFile(Paths.get(getOutputDir()), "redundant_", extension).toString();
			writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(output));

			// This code reads the sam/bam file line by line selecting batches
			// of
			// reads having the same id. Each batch is passed to the parseBatch
			// method.
			SAMRecord previous = null;
			List<SAMRecord> batch = new ArrayList<>();

			for (SAMRecord current : filtered) {

				if (previous != null) {
					if (!current.getReadName().equals(previous.getReadName())) {
						this.parseBatch(batch).stream().forEach(writer::addAlignment);
						batch.clear();
					}
				}

				previous = current;
				batch.add(previous);
			}
			this.parseBatch(batch).stream().forEach(writer::addAlignment);			
		} catch (IOException e) {
			log.severe("Cannot write redundant alignment file , " + e.getMessage());
		} finally {
			if (reader != null)
				try {
					reader.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (writer != null)
				writer.close();

		}
		return output;
	}

	/**
	 * Each batch of reads is divided into paired and unpaired. The paired reads
	 * are in turn divided into forward and reverse reads. Then each group is
	 * parsed using the {@link BowTie2#getMultimappedReads(List)} method.
	 * 
	 * @param batch
	 *            the batch of reads to parse
	 * @return a list of redundant reads
	 */
	private List<SAMRecord> parseBatch(List<SAMRecord> batch) {
		// Dividing paired reads from single reads
		Map<Boolean, List<SAMRecord>> whole = batch.stream()
				.collect(Collectors.groupingBy(SAMRecord::getReadPairedFlag));

		List<SAMRecord> paired = whole.get(true);
		List<SAMRecord> single = whole.get(false);

		ArrayList<SAMRecord> redundant = new ArrayList<>();

		// Parsing single reads with the getMultimappedReads method
		if (single != null && single.size() > 1) {
			redundant.addAll(this.getMultimappedReads(single));
		}

		// Dividing paired reads into forward and reverse and parsing them with
		// the getMultimappedReads method
		if (paired != null)
			paired.stream().collect(Collectors.groupingBy(SAMRecord::getFirstOfPairFlag)).values().stream()
					.filter(l -> l.size() > 1).map(this::getMultimappedReads).forEach(redundant::addAll);

		// Redundant reads are added to the redundant list
		return redundant;
	}

	/**
	 * Each batch of reads passed by the {@link BowTie2#parseBatch(List)} method
	 * will be parsed collecting redundant reads
	 * 
	 * @param list
	 *            a list of reads all with the same ID
	 * @return a list of redundant reads
	 */
	private List<SAMRecord> getMultimappedReads(List<SAMRecord> list) {
		// Soritng list (best record at the top)
		list = sortList(list);

		// Gets the first record (the best one if the input has been sorted
		// using the sortList method)
		SAMRecord first = list.get(0);

		// If the reads is not mapped in any other region of the reference
		// genome (XS flag is null) or the other alignment is not good as the
		// primary alignment (XS < AS) then an empty list will be returned
		if (first.getAttribute("XS") == null || !first.getAttribute("XS").equals(first.getAttribute("AS")))
			return new ArrayList<>();

		// Getting the best value of AS and YS (if the read is paired)
		final Object as = first.getAttribute("AS");
		final Object ys = first.getAttribute("YS");

		// Filtering out reads with an AS lower than the best AS obtained before
		list = list.stream().filter(r -> r.getAttribute("AS").equals(as)).collect(Collectors.toList());
		if (ys != null)
			// Filtering out alignments with an YS values lower than the best YS
			// obtained before
			list = list.stream().filter(r -> r.getAttribute("YS").equals(ys)).collect(Collectors.toList());

		// If the list of redundant alignments contains only one element an
		// empty list will be returned
		if (list.size() > 1)
			return list;

		return new ArrayList<>();
	}

	/**
	 * Sorts a list of records based on primary and secondary alignmens (primary
	 * first), quality mapping, the AS flag, the XS flag, and finally the YS
	 * flag of Bowtie2 aligners.
	 * 
	 * @param toSort
	 *            list of records to sort
	 * @return a sorted list in descending order
	 */
	private List<SAMRecord> sortList(List<SAMRecord> toSort) {
		// Comparator
		Comparator<SAMRecord> comparator = new Comparator<SAMRecord>() {

			@Override
			public int compare(SAMRecord o1, SAMRecord o2) {

				// Testing for primary alignment
				if (o1.isSecondaryOrSupplementary() && !o2.isSecondaryOrSupplementary())
					return -1;
				if (!o1.isSecondaryOrSupplementary() && o2.isSecondaryOrSupplementary())
					return 1;

				// Testing for mapping quality
				int res = o1.getMappingQuality() - o2.getMappingQuality();

				if (res == 0 && o1.getAttribute("AS") != null && o2.getAttribute("AS") != null) {
					Integer as1 = (Integer) o1.getAttribute("AS");
					Integer as2 = (Integer) o2.getAttribute("AS");

					// Testing AS flag
					res = as1.compareTo(as2);

					if (res == 0 && o1.getAttribute("XS") != null && o2.getAttribute("XS") != null) {
						Integer xs1 = (Integer) o1.getAttribute("XS");
						Integer xs2 = (Integer) o2.getAttribute("XS");

						// Testing XS flag
						res = xs1.compareTo(xs2);

						if (res == 0 && o1.getAttribute("YS") != null && o2.getAttribute("YS") != null) {
							Integer ys1 = (Integer) o1.getAttribute("YS");
							Integer ys2 = (Integer) o2.getAttribute("YS");

							// Testing YS flag
							res = ys1.compareTo(ys2);
						}
					}

				}

				// Return
				return res;
			}

		};

		Collections.sort(toSort, comparator);
		Collections.reverse(toSort);
		return toSort;
	}

	/**
	 * @return the forceRebuild
	 */
	public boolean isForceRebuild() {
		return forceRebuild;
	}

	/**
	 * @param forceRebuild the forceRebuild to set
	 */
	public void setForceRebuild(boolean forceRebuild) {
		this.forceRebuild = forceRebuild;
	}

}
