package cox.external;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Logger;

import cox.utils.ProcessUtils;

public class BedTools implements ExternalApplication {

	private String outputDir = ".";

	private boolean verbose = false;

	private String logFile = "bedtools.log";

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
	 * Computes genome coverage using standard format
	 * 
	 * @param bam
	 *            the alignment file
	 * @return the coverage map as reported by genomeCoverageBed command
	 */
	public String genomeCoverage(String bam) {
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "map_", ".coverage").toString();
			String cmd = String.format("genomeCoverageBed -ibam %s", bam);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcessToFile(pb, verbose, output, logFile);

			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s viewing sam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot compute coverage, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Covergae computation has been interrupted, " + e.getMessage());
		}
		return null;
	}

	/**
	 * Computes genome coverage using the per-base format
	 * 
	 * @param bam
	 *            the alignment file
	 * @return the coverage map as reported by genomeCoverageBed command with
	 *         the -d flag
	 */
	public String genomeCoveragePerBase(String bam) {
		try {
			String output = Files.createTempFile(Paths.get(getOutputDir()), "perbase_", ".coverage").toString();
			String cmd = String.format("genomeCoverageBed -d -ibam %s", bam);

			// Building process
			ProcessBuilder pb = new ProcessBuilder(cmd.split("\\s"));

			int res = ProcessUtils.runProcessToFile(pb, verbose, output, logFile);

			if (res != 0) {
				Files.deleteIfExists(Paths.get(output));
				throw new IOException("Error/s viewing sam alignment");
			} else {
				return output;
			}
		} catch (IOException e) {
			log.severe("Cannot compute coverage, " + e.getMessage());
		} catch (InterruptedException e) {
			log.warning("Covergae computation has been interrupted, " + e.getMessage());
		}
		return null;
	}

	/**
	 * Computes the mean coverage per contigs based on the coverage file
	 * generated by genomeCoverageBed with standard output format
	 * 
	 * @param coverageFile
	 *            the coverage file
	 * @return the mean coverage per sample
	 */
	public String getMeanCoveragePerContigs(String coverageFile) {
		String output = null;
		try (BufferedReader reader = Files.newBufferedReader(Paths.get(coverageFile))) {
			output = Files.createTempFile(Paths.get(getOutputDir()), "mean_", ".coverage").toString();

			BufferedWriter writer = new BufferedWriter(new FileWriter(output));

			String line = null;
			String pc = null;
			double cov = 0.0;

			while ((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if (pc == null)
					pc = fields[0];
				if (pc.equals(fields[0])) {
					cov += Double.valueOf(fields[1]) * Double.valueOf(fields[4]);
				} else {
					writer.write(String.format("%s\t%.3f", pc, cov));
					writer.newLine();
					writer.flush();
					cov = Double.valueOf(fields[1]) * Double.valueOf(fields[4]);
					pc = fields[0];
				}

			}
			// Do not write the last coverage as it is the "genome coverage"
			// introduced by bedtools
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return output;
	}

}
