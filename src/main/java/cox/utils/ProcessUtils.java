package cox.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class ProcessUtils {

	public static int runProcess(ProcessBuilder pb, boolean redirectError, boolean verbose, String log)
			throws IOException, InterruptedException {
		// Process and output reader
		Process process = null;
		BufferedReader br = null;
		BufferedWriter writer = null;
		int res = -1;
		pb.redirectErrorStream(redirectError);

		// Building process
		try {
			process = pb.start(); // starting process
			InputStreamReader isr = new InputStreamReader(process.getInputStream());
			br = new BufferedReader(isr);
			FileWriter logWriter = new FileWriter(log, true);
			writer = new BufferedWriter(logWriter);
			String line = null;
			while ((line = br.readLine()) != null) {
				writer.write(line);
				writer.newLine();
				writer.flush();
				if (verbose)
					System.out.println(line);
			}
			
			res = process.waitFor();
		} finally {
			// cleaning up everything
			if (process != null)
				process.destroy();
			if (br != null)
				try {
					br.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (writer != null)
				try {
					writer.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
		}
		return res;
	}

	public static int runProcessToFile(ProcessBuilder pb, boolean verbose, String file, String log)
			throws IOException, InterruptedException {
		// Process and output reader
		Process process = null;

		BufferedReader br = null;
		BufferedReader ebr = null;

		BufferedWriter writer = null;
		BufferedWriter output = null;

		int res = -1;
		pb.redirectErrorStream(false);

		// Building process
		try {
			process = pb.start(); // starting process
			InputStreamReader isr = new InputStreamReader(process.getInputStream());
			br = new BufferedReader(isr);

			InputStreamReader esr = new InputStreamReader(process.getErrorStream());
			ebr = new BufferedReader(esr);

			FileWriter logWriter = new FileWriter(log, true);
			FileWriter outWriter = new FileWriter(file);

			writer = new BufferedWriter(logWriter);
			output = new BufferedWriter(outWriter);

			String line = null;

			while ((line = br.readLine()) != null) {
				output.write(line);
				output.newLine();
				output.flush();
				if (verbose)
					System.out.println(line);
			}
			
			while ((line = ebr.readLine()) != null) {
				writer.write(line);
				writer.newLine();
				writer.flush();
				if (verbose)
					System.out.println(line);
			}

			res = process.waitFor();
		} finally {
			// cleaning up everything
			if (process != null)
				process.destroy();
			if (br != null)
				try {
					br.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (writer != null)
				try {
					writer.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (ebr != null)
				try {
					ebr.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
			if (output != null)
				try {
					output.close();
				} catch (IOException e) {
					// NOTHING TO DO
				}
		}
		return res;
	}
}
