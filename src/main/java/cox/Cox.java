package cox;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.DNASequenceCreator;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.util.InputStreamProvider;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import cox.external.BedTools;
import cox.external.BowTie2;
import cox.external.ExternalApplication;
import cox.external.PicardWrapper;
import cox.external.SamTools;
import cox.utils.NameFormatter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class Cox {

	// Main logger
	private static Logger log = Logger.getLogger("cox");

	private static final String USAGE = "" + "Java jar usage:" + System.lineSeparator()
			+ "java -jar cox.jar ((-1 <forward> -2 <reverse) |& -s <single>) -r <reference> -o <out_dir> [OPTIONS]"
			+ System.lineSeparator() + "Script usage:" + System.lineSeparator()
			+ "./cox ((-1 <forward> -2 <reverse) |& -s <single>) -r <reference> -o <out_dir> [OPTIONS]"
			+ System.lineSeparator();

	/**
	 * Help option
	 */
	@Option(name = "-h", aliases = "--help", usage = "print help message", help = true)
	private boolean help = false;

	/**
	 * Forward read file
	 */
	@Option(name = "-1", aliases = "--fwd", usage = "forward read file", depends = { "-2" })
	private String forward = null;

	/**
	 * Reverse read file
	 */
	@Option(name = "-2", aliases = "--rew", usage = "reverse read file", depends = { "-1" })
	private String reverse = null;

	/**
	 * Single read file
	 */
	@Option(name = "-s", aliases = "--single", usage = "single read file")
	private String single = null;

	/**
	 * References
	 */
	@Option(name = "-r", aliases = "--ref", usage = "reference sequence file", required = true)
	private String reference;

	/**
	 * Output directory
	 */
	@Option(name = "-o", aliases = "--output", usage = "the output directory", required = true)
	private String out;

	/**
	 * Verbose
	 */
	@Option(name = "-v", aliases = "--verbose", usage = "the output of all commands will be print to the consolle")
	private boolean verbose = false;

	/**
	 * Flag for coverage map computation
	 */
	@Option(name = "-c", aliases = "--cov-map", usage = "a map reporting the coverage value for each position in the reference genome/s will be generated")
	private boolean cov_map = false;

	/**
	 * Flag for mean coverage computation
	 */
	@Option(name = "-m", aliases = "--mean-cov", usage = "the average covergae for each reference genome/s will be reported")
	private boolean avg_cov = false;

	/**
	 * Flag for keeping temporary files
	 */
	@Option(name = "-k", aliases = "--keep-tmp", usage = "temp files will be saved in a temporary directory inside the output folder")
	private boolean keep = false;

	/**
	 * Filter alignment quality
	 */
	@Option(name = "-f", aliases = "--filter-qual", usage = "alignment with a mapping quality lower than this value will not be reported")
	private int filter = -1;

	/**
	 * Estimate redundancy
	 */
	@Option(name = "-d", aliases = "--redundancy", usage = "redundant alignments will be saved and analyzed following specified criteria (see -c and -m options)")
	private boolean red = false;

	/**
	 * Calculating GC count
	 */
	@Option(name = "-g", aliases = "--gc-count", usage = "GC count will be reported for each reference sequence in a tab separated file")
	private boolean gc_count = false;

	// External Applications
	private BowTie2 bw2 = new BowTie2();
	private PicardWrapper picard = new PicardWrapper();
	private SamTools samtools = new SamTools();
	private BedTools bedtools = new BedTools();

	private ExternalApplication[] applications = { bw2, picard, samtools, bedtools };

	/**
	 * List of temporary files
	 */
	private List<String> temp_files = new ArrayList<>();

	/**
	 * List of file to save
	 */
	private Map<String, String> toSave = new HashMap<>();

	/**
	 * Starts the analysis with the specified parameters
	 */
	public void startAnalysis() {

		// Set up applications and check
		this.setUpApplications();

		// Building and indexing reference database
		String faidx = this.buildDBandIndexing();

		// Computing GC count if needed
		this.getGCcontent();

		// Mapping
		String uber_map = this.map();

		// Filtering out low quality alignments if needed
		uber_map = this.filtering(uber_map);

		// Mark duplicates and remove them
		uber_map = this.markDuplicates(uber_map);

		// Estimating redundancy if needed
		this.redundancyEstimation(uber_map, faidx);

		// Removing supplementary alignments
		String uber_map_indexed = this.indexingAndSorting(faidx, uber_map, true);
		if (uber_map_indexed != null) {
			toSave.put(uber_map_indexed, "clean_alignment" + NameFormatter.getExtension(uber_map_indexed));
			uber_map = uber_map_indexed;
		}

		// Coverage estimation
		if (avg_cov) {
			String avg = meanCoverageEstimation(uber_map);
			if (avg != null) {
				toSave.put(avg, "mean.coverage");
			}
		}

		if (cov_map) {
			String map = perBaseCoverageEstimation(uber_map);
			if (map != null) {
				toSave.put(map, "map.coverage");
			}
		}

		// Clear temp and save
		clearTempAndCollect(keep);
	}

	/**
	 * Computes gc count for each reference sequence if needed
	 */
	private void getGCcontent() {
		if (gc_count) {
			log.info("Calculating GC content for reference sequences");
			InputStreamProvider isp = new InputStreamProvider();
			InputStream is = null;
			Path output = null;
			try {
				output = Files.createTempFile(Paths.get(out), "gc_", ".metrics");
				is = isp.getInputStream(reference);
				FastaReader<DNASequence, NucleotideCompound> reader = new FastaReader<>(is,
						new GenericFastaHeaderParser<>(), new DNASequenceCreator(DNACompoundSet.getDNACompoundSet()));

				LinkedHashMap<String, DNASequence> seq = null;
				while ((seq = reader.process(10)) != null) {
					List<String> results = seq.entrySet().stream().map(e -> {
						DNASequence s = e.getValue();
						double gc = (double) s.getGCCount() / (double) s.getLength();
						String res = String.format("%s\t%.2f", e.getKey(), gc);
						return res;
					}).collect(Collectors.toList());
					Files.write(output, results, StandardOpenOption.APPEND);
				}
				toSave.put(output.toString(), "gc_reference.tsv");
			} catch (IOException e) {
				if (output != null)
					try {
						Files.deleteIfExists(output);
					} catch (IOException e1) {
						// NOTHING TO DO
					}
				log.info("Error/s computing GC count:" + e.getMessage());
			}
		}
	}

	/**
	 * Sets up all external applications
	 */
	private void setUpApplications() {

		// Chack if output folder exists
		Path output = Paths.get(out);
		if (Files.notExists(output)) {
			log.info("Output folder does not exist, trying to create it");
			try {
				Files.createDirectory(output);
			} catch (IOException e) {
				log.log(Level.SEVERE, "Cannot create output directory {0}: {1}", new String[] { out, e.getMessage() });
			}
		}

		log.log(Level.INFO, "Writing output to: {0}", out);
		log.log(Level.INFO, "Setting verbose to: {0}", verbose);

		Arrays.stream(applications).forEach(app -> {
			app.setOutputDir(out);
			app.setVerbose(verbose);
		});

		// Checking inptu files
		if (forward == null && reverse == null && single == null) {
			log.severe("You have to specify at least a single read file");
			System.exit(-1);
		}
	}

	/**
	 * Mean coverage estimation
	 * 
	 * @param alignment
	 *            an alignment sorted and indexed based on reference
	 * @see Cox#indexingAndSorting(String, String, boolean)
	 * @return the name of the mean coverage file or <code>null</code> if
	 *         something did not work
	 */
	private String meanCoverageEstimation(String alignment) {
		String mean_coverage = null;
		log.info("Computing mean coverage per reference");
		String cov = bedtools.genomeCoverage(alignment);

		if (cov == null) {
			log.warning("Cannot compute mean coverage per reference, see BedTools log for details");
		} else {
			temp_files.add(cov);
			mean_coverage = bedtools.getMeanCoveragePerContigs(cov);
			if (mean_coverage == null) {
				log.warning("Cannot compute mean coverage per reference, see BedTools log for details");
			}
		}
		return mean_coverage;
	}

	/**
	 * Per base coverage estimation
	 * 
	 * @param alignment
	 *            an alignment sorted and indexed based on reference
	 * @see Cox#indexingAndSorting(String, String, boolean)
	 * @return the name of the coverage map file or <code>null</code> if
	 *         something did not work
	 */
	private String perBaseCoverageEstimation(String alignment) {
		String coverage_map = null;
		log.info("Computing coverage map for each reference sequence");

		coverage_map = bedtools.genomeCoveragePerBase(alignment);

		if (coverage_map == null) {
			log.warning("Cannot compute coverage map, see BedTools log for details");
		}
		return coverage_map;
	}

	/**
	 * Indexing and sorting alignment.
	 * 
	 * @param faidx
	 *            the faidx of the reference
	 * @param alignment
	 *            the alignment
	 * @param removeSupplementaryAlignments
	 *            if <code>true</code> supplementary alignments will be removed
	 * @return a clean alignment
	 */
	private String indexingAndSorting(String faidx, String alignment, boolean removeSupplementaryAlignments) {
		// Removing supplementary alignments if needed
		if (removeSupplementaryAlignments) {
			log.info("Removing supplementary alignments and indexing based on reference");
		} else {
			log.info("Indexing based on reference");
		}

		String bam = (removeSupplementaryAlignments) ? samtools.view(faidx, alignment, "-F 0x0100")
				: samtools.view(faidx, alignment);

		if (bam == null && removeSupplementaryAlignments) {
			log.warning("Cannot remove supplementray alignments, coverage could be overestimated");
			bam = alignment;
		} else if (bam == null && !removeSupplementaryAlignments) {
			log.warning("Cannot determine sort order based on reference");
			bam = alignment;
		} else {
			temp_files.add(bam);
		}

		// Sorting and indexing
		log.info("Sorting and indexing");

		String bam_smds = samtools.sort(bam);

		// If sorting didn't work logging and exit
		if (bam_smds == null) {
			log.severe("Cannot sort alignments, see Samtools log for details");
			clearTempAndCollect(keep);
			System.exit(-1);
		}

		String bam_smds_index = samtools.index(bam_smds);

		// If indexing didn't work, continue without indexing (slower analyses)
		if (bam_smds_index == null) {
			log.warning("Problem/s indexing alignments");
		} else {
			temp_files.add(bam_smds_index);
		}

		// Return sorted file
		return bam_smds;
	}

	/**
	 * Estimates redundancy of an alignment file. If the redundancy flag has not
	 * been specified this method will not do anything
	 * 
	 * @param uber_map
	 *            the alignment file
	 * @param faidx
	 *            index file
	 */
	private void redundancyEstimation(String uber_map, String faidx) {
		if (red) {
			// Estimating redundancy
			log.info("Estimating redundancy");

			// Get redundant alignments
			String redundant = bw2.getRedundantAlignments(uber_map);
			if (redundant == null) { // If redundancy esitmation did not
										// succeede
										// return
				log.warning("Cannot estimate redundancy, redundancy file will not be reported");
				return;
			}

			// Open redundant file
			SamReader reader = SamReaderFactory.makeDefault().open(Paths.get(redundant));

			// Getting the number of redundant alignments per read
			Collection<Long> results = reader.iterator().stream()
					.collect(Collectors.groupingBy(SAMRecord::getReadName, Collectors.counting())).values();

			// If redundant alignment file is empty, return
			if (results.size() == 0) {
				log.info("No redundant alignment was found");
				temp_files.add(redundant);
				return;
			}

			// Open input file
			SamReader uber_reader = SamReaderFactory.makeDefault().open(Paths.get(uber_map));

			// Getting the number of aligned sequences
			int uber_size = uber_reader.iterator().stream().filter(s -> !s.getReadUnmappedFlag())
					.collect(Collectors.groupingBy(SAMRecord::getReadName)).keySet().size();

			// Statistics of redundant alignment
			DoubleSummaryStatistics statsWithin = results.stream().mapToDouble(Long::doubleValue).summaryStatistics();

			// Percentage of redundant alignments
			double redundant_perc = (double) statsWithin.getCount() / (double) uber_size;

			// Redundancy reate per alignment
			double redundant_rate = (double) ((uber_size - statsWithin.getCount()) + statsWithin.getSum())
					/ (double) uber_size;

			// Formatting output
			List<String> stats = new ArrayList<>();
			stats.add(String.format("Overall:"));
			stats.add(String.format("Redundant alignments: %d (%.2f%%)", statsWithin.getCount(), redundant_perc));
			stats.add(String.format("Redundant alignment rate: %.2f%n", redundant_rate));
			stats.add(String.format("Redundant alignments:"));
			stats.add(String.format("Average redundancy: %.2f", statsWithin.getAverage()));
			stats.add(String.format("Minimum number of redundant alignments: %.0f", statsWithin.getMin()));
			stats.add(String.format("Maximum number of redundant alignments: %.0f", statsWithin.getMax()));

			try {
				// Writing stats
				log.info("Writing redundancy stats");
				Path output = Files.createTempFile(Paths.get(out), "redundancy_", ".metrics");
				Charset cs = Charset.defaultCharset();

				Files.write(output, stats, cs, StandardOpenOption.APPEND);

				if (Files.exists(output)) // If output file exists mark it for
											// saving
					toSave.put(output.toString(), "redundancy.stats");
			} catch (IOException e) {
				log.severe("Cannot write redundancy stats: " + e.getMessage());
			}

			// Indexing and sorting redundant alignments. This step needs to be
			// done
			// for coverage estimation
			log.info("Indexing and sorting redundant alignment file");
			// log.log(Level.INFO, "FAIDX = {0}, REDUNDANT = {1}", new String[]
			// {
			// faidx, redundant });
			String indexed = this.indexingAndSorting(faidx, redundant, false);

			if (indexed == null) { // If indexing didn't work return
				log.warning("Cannot indexing redundant alignment file");
				return;
			}

			// Add redundant file to temp list, and mark the indexed file for
			// saving
			temp_files.add(redundant);
			toSave.put(indexed, "redundancy_alignments" + NameFormatter.getExtension(indexed));

			// Coverage estimation
			if (avg_cov) {
				String avg = meanCoverageEstimation(indexed);
				if (avg != null) {
					toSave.put(avg, "redundancy_mean.coverage");
				}
			}

			if (cov_map) {
				String map = perBaseCoverageEstimation(indexed);
				if (map != null) {
					toSave.put(map, "redundancy_map.coverage");
				}
			}
		}
	}

	/**
	 * Removes duplicated alignments
	 * 
	 * @param duplicated
	 *            the file to check for duplicates
	 * @return a file with no duplicated alignments
	 */
	private String markDuplicates(String duplicated) {
		// Mark duplicates and remove them
		log.info("Removing duplicated reads from alignment");

		// Removing duplicated sequences assume sorting by queryname
		String[] md = picard.markDuplicates(duplicated, "queryname");

		String alignment = duplicated;
		if (md == null) {
			log.warning("Duplicates were not removed due to error/s, see Picard log for details");
		} else {
			alignment = md[0];
			toSave.put(md[1], "remove_duplicates.metrics");
			temp_files.add(alignment);
		}

		return alignment;
	}

	/**
	 * Removes alignments that have a lower mapping quality than the one
	 * specified in the option. If no quality cutoff has been specified, this
	 * method will not do anything.
	 * 
	 * @param alignment
	 *            the alignment file
	 * @return the name of the filtered file
	 */
	private String filtering(String alignment) {
		// Filtering out low quality alignments if needed
		if (filter >= 0) {

			log.log(Level.INFO, "Filtering out low quality alignments (quality cutoff = {0})", filter);

			String filtered = bw2.filterLowQualityRecords(alignment, filter);
			if (filtered == null) { // If filtering didn't succeede returning
									// the input file
				log.warning("Cannot filter records based on mapping quality");
			} else {
				temp_files.add(filtered);
				alignment = filtered;
			}

		}

		return alignment;
	}

	/**
	 * Maps given read files onto the reference genome provided
	 * 
	 * @return the name of the file with the alignment in sam format
	 */
	private String map() {
		// Mapping
		String uber_map = null;

		log.info("Start mapping");

		// If forward, reverse and single files have been specified, mapping all
		// files
		if (forward != null && reverse != null && single != null) {

			log.log(Level.INFO, "Forward reads: {0}", Paths.get(forward).toAbsolutePath().toString());
			log.log(Level.INFO, "Reverse reads: {0}", Paths.get(reverse).toAbsolutePath().toString());
			log.log(Level.INFO, "Single reads: {0}", Paths.get(single).toAbsolutePath().toString());

			String single_map = bw2.mapSingleEnd(single, "-a");
			String paired_map = bw2.mapPairedEnd(forward, reverse, "-a");

			// Handling errors
			if (single_map != null && paired_map == null) { // Single mapped,
															// paired not
				uber_map = single_map;
				log.warning("Problems mapping paired reads, using only single alignments");
			} else if (single_map == null && paired_map != null) { // Paired
																	// mapped,
																	// single
																	// not
				uber_map = paired_map;
				log.warning("Problems mapping single reads, using only paired alignments");
			} else if (single_map != null && paired_map != null) { // Both
																	// mapped

				log.info("Merging paired and single files");

				// Merging
				uber_map = picard.mergeMaps(paired_map, single_map);
				if (uber_map == null) { // Not merged, using paired data and
										// adding single map to temp files
					log.warning("Problems merging paired and single reads, using only paired alignments");
					uber_map = paired_map;
					temp_files.add(single_map);
				} else { // Merged, adding paired and single file to temp files
					temp_files.add(single_map);
					temp_files.add(paired_map);
				}
			}
			// If only single reads have beeen specified, mapping only single
			// reads
		} else if (forward == null && reverse == null && single != null) {

			log.log(Level.INFO, "Single reads: {0}", Paths.get(single).toAbsolutePath().toString());

			uber_map = bw2.mapSingleEnd(single, "-a");

			// If only paired files have been specified, mapping only paired
			// reads
		} else if (forward != null && reverse != null && single == null) {

			log.log(Level.INFO, "Forward reads: {0}", Paths.get(forward).toAbsolutePath().toString());
			log.log(Level.INFO, "Reverse reads: {0}", Paths.get(reverse).toAbsolutePath().toString());

			uber_map = bw2.mapPairedEnd(forward, reverse, "-a");
		}

		// If mapping didn't work log message and leave
		if (uber_map == null) {
			log.severe("Cannot map reads to reference genome/s, see Bowtie2 log for details");
			clearTempAndCollect(keep);
			System.exit(-1);
		}

		// Adding mapping to temp files and returning the mapped file
		temp_files.add(uber_map);
		return uber_map;
	}

	/**
	 * Builds a Bowtie2 database and index it.
	 * 
	 * @return the index file
	 */
	private String buildDBandIndexing() {
		// Building reference database
		log.log(Level.INFO, "Building Bowtie2 database: {0}", Paths.get(reference).toAbsolutePath().toString());

		// If building didn't work, logging, removing temp files, and leave
		if (!bw2.buildReference(reference)) {
			log.severe("Cannot build Bowtie2 database");
			this.clearTempAndCollect(keep);
			System.exit(-1);
		}

		// Indexing DB
		log.info("Indexing reference");

		String faidx = samtools.faidx(reference);

		// If idexing didn't work, logging the event and continue
		if (faidx == null) {
			log.warning("Cannot index reference sequence/s");
		} else {
			temp_files.add(faidx);
		}

		// Return index file
		return faidx;
	}

	/**
	 * Clear temp file or move them in a temporary folder
	 * 
	 * @param files
	 *            the list of temp file/s to remove
	 * @param keep
	 *            if <code>true</code> temp files will be stored in a temporary
	 *            folder
	 */
	private void clearTempAndCollect(boolean keep) {
		if (!keep) {
			log.info("Removing temporary files");
			temp_files.stream().filter(s -> s != null).map(Paths::get).forEach(t -> {
				try {
					Files.deleteIfExists(t);
				} catch (IOException e) {
					log.severe("Cannot delete temporary file/s: " + e.getMessage());
				}
			});
		} else {
			try {
				Path tmp_dir = Files.createTempDirectory(Paths.get(out), "tmp");
				log.log(Level.INFO, "Moving temporary files into {0}", tmp_dir);
				temp_files.stream().filter(s -> s != null).map(Paths::get).forEach(p -> {
					Path target = tmp_dir.resolve(p.getFileName());
					try {
						Files.move(p, target, StandardCopyOption.REPLACE_EXISTING);
					} catch (IOException e) {
						log.severe("Cannot move temporary file/s: " + e.getMessage());
					}
				});
			} catch (IOException e) {
				log.severe("Cannot create temporary directory: " + e.getMessage());
			}

		}

		toSave.forEach((k, v) -> {
			Path old = Paths.get(k);
			Path nnew = Paths.get(v);

			try {
				Files.move(old, old.resolveSibling(nnew));
			} catch (IOException e) {
				log.log(Level.WARNING, "Cannot rename {0} in {1}, using old name", new String[] { k, v });
			}
		});
	}

	@SuppressWarnings("unused")
	public static void main(String[] args) {

		// Setting log
		log.setUseParentHandlers(false);
		ConsoleHandler handler = new ConsoleHandler();
		handler.setFormatter(new ConsoleLogFormatter());
		log.addHandler(handler);

		// Creating main class
		Cox cox = new Cox();

		// Parsing args
		CmdLineParser parser = new CmdLineParser(cox);
		try {
			parser.parseArgument(args);

			// Print help and exit
			if (cox.help) {
				System.out.println(USAGE);
				parser.printUsage(System.out);
				System.exit(0);
			}

			// Start analysis
			cox.startAnalysis();
		} catch (CmdLineException e) {
			System.err.println(USAGE);
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.exit(-1);
		}

		// Original Pipeline
		// TODO remove this part after debugging
		if (false) {
			String out = Paths.get(args[0]).getParent().toString();
			BowTie2 bw2 = new BowTie2();
			bw2.setVerbose(true);
			bw2.setOutputDir(out);

			PicardWrapper picard = new PicardWrapper();
			picard.savePicardPathFromProperties(args[0]);
			picard.setVerbose(true);
			picard.setOutputDir(out);

			SamTools samtools = new SamTools();
			samtools.setVerbose(true);
			samtools.setOutputDir(out);

			BedTools bedtools = new BedTools();
			bedtools.setVerbose(true);
			bedtools.setOutputDir(out);

			// Building reference DB
			bw2.buildReference(args[0]);

			// Mapping paired and unpaired
			String map1 = bw2.mapPairedEnd(args[1], args[2], "-a");
			String map2 = bw2.mapSingleEnd(args[3], "-a");
			// Merging mapped files
			String map = picard.mergeMaps(map1, map2);

			// Indexing DB
			String faidx = samtools.faidx(args[0]);

			// Removing supplementary alignments
			String bam = samtools.view(faidx, map, "-F 0x0100");

			// Sorting and indexing
			String bam_s = samtools.sort(bam);
			samtools.index(bam_s);

			// Getting stats
			samtools.flagStats(bam_s);

			// Removing duplicates, sorting and indexing
			String[] mdresults = picard.markDuplicates(bam_s, "queryname");
			String bam_smd = mdresults[0];
			String bam_smds = samtools.sort(bam_smd);
			samtools.index(bam_smds);

			// Getting stats again
			samtools.flagStats(bam_smds);

			// Coverage estimation
			String cov = bedtools.genomeCoverage(bam_smds);
			// bedtools.genomeCoveragePerBase(bam_smds);
			bedtools.getMeanCoveragePerContigs(cov);

			System.out.println("results in: " + bam_smds);

			bw2.getRedundantAlignments(args[0]);
		}
	}

}
