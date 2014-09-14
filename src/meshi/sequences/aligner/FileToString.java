package meshi.sequences.aligner;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 30/07/2005
 * Time: 14:29:19
 * To change this template use File | Settings | File Templates.
 */
public class FileToString {
	String str = new String();

	public FileToString() {
	}

	public FileToString(String file) {

		try {
			FileReader fr = new FileReader(file);
			BufferedReader is = new BufferedReader(fr);
			String s1 = is.readLine();
			while (s1 != null) {
				str = str.concat("\n" + s1);
				s1 = is.readLine();
			}
			is.close();

		}

		catch (IOException e) {
			System.out.print("Error: while reading" + e);
			System.exit(2);
		}
	}

	public String getText() {
		return str;
	}
}