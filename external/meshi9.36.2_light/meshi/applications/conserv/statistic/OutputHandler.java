package meshi.applications.conserv.statistic;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

public class OutputHandler {
	
	private static final String SEPERATOR = ", ";
	
	BufferedWriter	bufferedWriter;
	List<String> 	options;
	
	public OutputHandler(String fileName, List<String> headers, String xValue) throws IOException {
		this.options = headers;
		bufferedWriter = new BufferedWriter(new FileWriter(new File(fileName)));
		String headersString = "File".concat(SEPERATOR.concat(xValue));
		for (String header : headers) {
			headersString = headersString.concat(SEPERATOR.concat(header));
		}
		headersString = headersString.concat("\n");
		bufferedWriter.write(headersString);
	}
	
	public void writeLine(double rms, Map<String, Double> results, String row_name) throws IOException {
		String result = row_name.concat(SEPERATOR.concat(String.valueOf(rms)));
		for (String i : options) {
			Double element = results.get(i);
			if (element != null) {
				result = result.concat(SEPERATOR.concat(String.valueOf(element)));
			} else {
				result = result.concat(SEPERATOR.concat(String.valueOf("null")));
			}
		}
		result = result.concat("\n");
		bufferedWriter.write(result);
	}
	
	public void close() throws IOException {
		bufferedWriter.close();
	}
}
