package cox.external;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;
import java.util.logging.Logger;

import cox.utils.NameFormatter;
import cox.utils.ProcessUtils;

/**
 * Wrapper for PicardTools
 * 
 * @author <a href="http://www.unifi.it/dblage/CMpro-v-p-65.html">Giovanni
 *         Bacci</a>
 *
 */
public class PicardWrapper implements ExternalApplication {

	/**
	 * Path to Picard JAR file
	 */
	private final String pathToPicard;

	/**
	 * Verbose property
	 */
	private boolean verbose = false;

	/**
	 * The name of the log file
	 */
	private String logfile = "picard.log";

	/**
	 * The logger
	 */
	private Logger log = Logger.getLogger(this.getClass().getName());

	/**
	 * The name of the property file
	 */
	private final String PROPERTY_FILE = "picard.properties";

	/**
	 * The name of the picard.jar path property
	 */
	private final String PICARD_HOME = "picard_home";

	/**
	 * Output directory
	 */
	private String outputDir = ".";

	/**
	 * Creates a {@link PicardWrapper} with the specified path to picard.jar
	 * 
	 * @param pathToPicard
	 *            the path to picard.jar
	 */
	public PicardWrapper(String pathToPicard) {
		this.pathToPicard = pathToPicard;
	}

	/**
	 * Creates a {@link PicardWrapper} using the path to picard.jar specified in
	 * the property file
	 */
	public PicardWrapper() {
		this.pathToPicard = loadPicardPathFromProperties();
	}

	/**
	 * Load the picard_home property from the file specified
	 * 
	 * @return the path to picard.jar
	 */
	private String loadPicardPathFromProperties() {
		Properties props = new Properties();
		InputStream is = null;

		String res = "picard.jar";

		try {
			// is = PicardWrapper.class.getResourceAsStream(PROPERTY_FILE);
			is = PicardWrapper.class.getClassLoader().getResourceAsStream(PROPERTY_FILE);
			props.load(is);
			res = props.getProperty(PICARD_HOME);
		} catch (IOException e) {
			log.severe("Cannot load Picard settings from "
					+ PicardWrapper.class.getClassLoader().getResource(PROPERTY_FILE).getFile());
		} finally {
			if (is != null) {
				try {
					is.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			}
		}
		return res;
	}

	/**
	 * Saves the path to picard.jar in the property file
	 * 
	 * @param pathToPicard
	 *            the path to picard.jar
	 */
	public void savePicardPathFromProperties(String pathToPicard) {
		Properties props = new Properties();
		props.setProperty(PICARD_HOME, pathToPicard);
		// File o = new
		// File(PicardWrapper.class.getResource(PROPERTY_FILE).getFile());
		File o = new File(PicardWrapper.class.getClassLoader().getResource(PROPERTY_FILE).getFile());
		OutputStream out;
		try {
			out = new FileOutputStream(o);
			props.store(out, "");
		} catch (FileNotFoundException e) {
			this.log.severe("Cannot write Picard properties to " + PROPERTY_FILE + ", " + e.getMessage());
		} catch (IOException e) {
			this.log.severe("I/O error writing Picard properties to " + PROPERTY_FILE + ", " + e.getMessage());
		}
	}

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

	public String mergeMaps(String map1, String map2) {
		String out = null;
		try {
			out = Files.createTempFile(Paths.get(getOutputDir()), "merged_", ".sam").toString();
			String cmd = String.format(
					"java -jar %s MergeSamFiles I=%s I=%s O=%s SORT_ORDER=queryname USE_THREADING=true", pathToPicard,
					map1, map2, out);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logfile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(out));
				throw new IOException("Error/s merging sam files");
			} else {
				return out;
			}
		} catch (IOException e) {
			log.severe("Cannot merge sam files, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Merging has been interrupted");
		}
		return out;
	}

	/**
	 * Search for and removes duplicate sequences from bam/sam alignment file.
	 * 
	 * @param bam
	 *            the alignment
	 * @param sort_order
	 *            one of [unsorted, queryname, coordinate, duplicate]
	 * @return an array containing the output file (index 0) and the metrics
	 *         file (index 1)
	 */
	public String[] markDuplicates(String bam_or_sam, String sort_order) {
		String[] results = new String[2];
		String extension = NameFormatter.getExtension(bam_or_sam);
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "markduplicates_", extension)
					.toString();
			String metrics = output.replaceFirst("\\." + extension + "$", ".metrics");

			String cmd = String.format("java -jar %s MarkDuplicates INPUT=%s OUTPUT=%s METRICS_FILE=%s "
					+ "ASSUME_SORT_ORDER=%s VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 "
					+ "REMOVE_DUPLICATES=TRUE", pathToPicard, bam_or_sam, output, metrics, sort_order);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logfile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s removing duplicate sequences");
			} else {
				results[0] = output;
				results[1] = metrics;
				return results;
			}
		} catch (IOException e) {
			log.severe("Cannot mark duplicate sequences, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Mark duplicates has been interrupted");
		}
		return null;
	}

}
