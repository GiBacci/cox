package cox.external;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import cox.utils.ProcessUtils;

/**
 * Samtools
 * 
 * @author <a href="http://www.unifi.it/dblage/CMpro-v-p-65.html">Giovanni
 *         Bacci</a>
 *
 */
public class SamTools implements ExternalApplication {

	private String outputDir = ".";

	private boolean verbose = false;

	private String logFile = "samtools.log";

	private final Logger log = Logger.getLogger(this.getClass().getName());

	@Override
	public void setOutputDir(String outputDir) {
		this.outputDir = outputDir;
	}

	@Override
	public String getOutputDir() {
		return outputDir;
	}

	@Override
	public boolean isVerbose() {
		return verbose;
	}

	@Override
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * computes the faidx of a fasta file
	 * 
	 * @param reference
	 *            a fasta file
	 * @return the name of the faidx file
	 */
	public String faidx(String reference) {
		String cmd = String.format("samtools faidx %s", reference);
		ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));
		String output = reference + ".fai";
		try {
			int res = ProcessUtils.runProcess(pb, true, verbose, logFile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error generating fai index");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("I/O error generating fai index, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Samtools faidx has been interrupted, " + e.getMessage());
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * converts a sam file into bam format
	 * 
	 * @param faiIndex
	 *            the faidx index
	 * @param sam
	 *            the alignment file
	 * @return a bam file
	 */
	public String view(String faiIndex, String sam, String... optionalArgs) {
		String opt = Arrays.stream(optionalArgs).collect(Collectors.joining(" "));
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "view_", ".bam").toString();
			String cmd = (optionalArgs.length > 0)
					? String.format("samtools view %s -bt %s -o %s %s", opt, faiIndex, output, sam)
					: String.format("samtools view -bt %s -o %s %s", faiIndex, output, sam);
			// String cmd = String.format("samtools view %s -bt %s -o %s %s",
			// opt, faiIndex, output, sam);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logFile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s viewing sam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot view alignment with samtools, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Viewing alignment has been interrupted, " + e.getMessage());
		}
		return null;
	}

	/**
	 * Sorts a bam file
	 * 
	 * @param bam
	 *            the alignment file
	 * @return the sorted bam file
	 */
	public String sort(String bam) {
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "sorted_", ".bam").toString();
			String cmd = String.format("samtools sort -o %s %s", output, bam);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logFile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s sorting bam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot sort alignment with samtools, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Sorting alignment has been interrupted, " + e.getMessage());
		}
		return null;
	}

	/**
	 * Creates a bai index for a bam file
	 * 
	 * @param bam
	 *            the bam file
	 * @return the fai index
	 */
	public String index(String bam) {
		try {
			String cmd = String.format("samtools index %s", bam);
			String output = bam + ".bai";

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcess(pb, true, verbose, logFile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s indexing bam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot index alignment with samtools, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Indexing alignment has been interrupted, " + e.getMessage());
		}
		return null;
	}

	/**
	 * Gets stats from a bam file
	 * 
	 * @param bam
	 *            the bam file
	 * @return the stat file
	 */
	public String flagStats(String bam) {
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "stats_", ".txt").toString();
			String cmd = String.format("samtools flagstat %s", bam);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcessToFile(pb, verbose, output, logFile);
			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s getting stats bam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot get stats from alignment with samtools, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Flagstat has been interrupted, " + e.getMessage());
		}
		return null;
	}

}
