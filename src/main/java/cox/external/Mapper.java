package cox.external;

/**
 * Interface for mapper objects (Bowtie2, BWA ...)
 * 
 * @author <a href="http://www.unifi.it/dblage/CMpro-v-p-65.html">Giovanni
 *         Bacci</a>
 *
 */
public interface Mapper {

	/**
	 * Maps paired reads
	 * 
	 * @param forward
	 *            file with mate1 sequences
	 * @param reverse
	 *            file with mate2 sequences
	 * @param otherOpts
	 *            other options.
	 * @return the name of the output file produced or <code>null</code> if
	 *         something went wrong
	 */
	public String mapPairedEnd(String forward, String reverse, String... otherOpts);

	/**
	 * Maps single end reads
	 * 
	 * @param single
	 *            the reads file
	 * @param otherOpts
	 *            other options. 
	 * @return the name of the output file produced or <code>null</code> if
	 *         something went wrong
	 */
	public String mapSingleEnd(String single, String... otherOpts);

	/**
	 * Builds the reference database. Once built the database will be used for
	 * mapping
	 * 
	 * @param reference
	 *            the path to the sequence file that must be used as reference
	 * @return <code>true</code> if the database has been built
	 */
	public boolean buildReference(String reference);
}
