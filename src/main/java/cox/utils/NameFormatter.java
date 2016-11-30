package cox.utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Util class for formatting file names. The base name depth refers to the
 * number of strings separated by punctuation characters in the file name. <br>
 * 
 * A name like this one:<br>
 * CF_ABB01SS_t2M16_CTCTCTAC-TATCCTCT_quality_control_screenSensLoc_hg19_screenSensLoc_phiX174_trim_adapters.fastq.bz2<br>
 * 
 * Has a basename depth of 12
 * 
 * @author <a href="http://www.unifi.it/dblage/CMpro-v-p-65.html">Giovanni
 *         Bacci</a>
 *
 */
public class NameFormatter {

	/**
	 * Gets the base name of a file based on the number of special characters in
	 * its name. If the depth os too big thsi method will return the original
	 * name fo the file<br>
	 * e.g.: foo_bar_bar.foo with depth 2 becomes -> foo_bar <br>
	 * foo_bar_bar.foo with depth 0 becomes -> foo <br>
	 * foo_bar_bar.foo with depth 4 (too big) -> foo_bar_bar <br>
	 * 
	 * @param name
	 *            the name
	 * @param depth
	 *            the number of special character
	 * @return the basename of the file
	 */
	public static String getBaseName(String name, int depth) {
		Pattern p = Pattern.compile("\\p{Punct}+");

		String baseName = getBaseName(name);

		Matcher m = p.matcher(baseName);

		int index = 0;

		boolean found = false;
		while (m.find()) {
			if (index == depth) {
				index = m.start();
				found = true;
				break;
			}
			index++;
		}

		if (!found)
			return baseName;

		return baseName.substring(0, index);
	}

	/**
	 * Returns the name of the file without the extension. If the file has not
	 * an extension this method will return the name without modifying it.
	 * 
	 * @param name
	 *            the name of the file
	 * @return the name of the file without the extension
	 */
	public static String getBaseName(String name) {
		int ext = name.lastIndexOf(".");
		if (ext >= 0)
			return name.substring(0, ext);
		return name;
	}

	/**
	 * Returns the extension of the file. If the file has no extension this
	 * method will return <code>null</code>.
	 * 
	 * @param name
	 *            the name of the file
	 * @return the extension
	 */
	public static String getExtension(String name) {
		int ext = name.lastIndexOf(".");
		if (ext >= 0)
			return name.substring(ext);
		return null;
	}

	/**
	 * Returns <code>true</code> if the file has an extension
	 * 
	 * @param name
	 *            the file name
	 * @return <code>true</code> if the file has an extension
	 */
	public static boolean hasExtension(String name) {
		return name.lastIndexOf('.') >= 0;
	}

	/**
	 * Finds the greatest common prefix of two strings
	 * 
	 * @param a
	 *            the first string
	 * @param b
	 *            the second string
	 * @return the greatest common prefix
	 */
	public static String greatestCommonPrefix(String a, String b) {
		int minLength = Math.min(a.length(), b.length());
		for (int i = 0; i < minLength; i++) {
			if (a.charAt(i) != b.charAt(i)) {
				return a.substring(0, i);
			}
		}
		return a.substring(0, minLength);
	}
}
