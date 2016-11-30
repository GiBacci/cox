package cox.external;

public interface ExternalApplication {
	
	public void setOutputDir(String outputDir);
	
	public String getOutputDir();
	
	public boolean isVerbose();

	public void setVerbose(boolean verbose);
}
