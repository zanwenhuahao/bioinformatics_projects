import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;


public class KmerFreqTreat {

	private static String AAvalid = "ACDEFGHIKLMNPQRSTVWY";
	private static int totalfivemer = 181535219;
	
	public static void main(String[] args) throws IOException {
		PrintWriter writer = new PrintWriter("kmerFreqTreated.txt");
		writer.println("Total "+ totalfivemer);
		
		File kmerfreqfile = new File("kmerFreq.txt");
		BufferedReader br = new BufferedReader(new FileReader(kmerfreqfile));
		String line;
		
		while((line = br.readLine())!=null)
		{
			;
			String[]linesep = line.trim().split("\\s+");
			
			boolean contains = false;
			for(int i=0; i<linesep[0].length() && !contains; i++)
			{
				char curchar = linesep[0].charAt(i);
				if(curchar=='B'||curchar=='J'||curchar=='O'||curchar=='U'||curchar=='X'||curchar=='Z')
				{
					contains = true;
				}
			}
			
			int fivemerfreq = Integer.parseInt(linesep[1]);
			
			if(!contains && fivemerfreq>=110)
			{
				writer.println(line);
			}
		}
		writer.close();
		 
		
	}

}
