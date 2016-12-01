import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Enumeration;
import java.util.Hashtable;



public class pepetidesIdentification {
	//The amino acids frequency in the entire uniprot_Swiss_protein file
	/*I found that taking 100,000 randomly generated 20-50 character long lines of true amino acids
	Using the sampler.java (i.e. remove the decoy ones) gives a good approximation up to 4 decimal points
	Of the frequency of the 20 amino acids.*/ 
	private static Hashtable<String, Double> aatable;
	
	//The 3-mer probability table
	private static Hashtable<String, Double> kmertable;
	private static long totalkmer = 0;
	private static double totalloglikelihoodscore = 0;
	
	static void Makeaatable() throws NumberFormatException, IOException
	{
		File AAfreqfile = new File("AAFreq.txt");

		BufferedReader br =  new BufferedReader(new FileReader(AAfreqfile));
		String line;
		aatable = new Hashtable<String, Double>();
		
		while((line = br.readLine()) !=null)
		{
			if(line.trim().length()>0)
			{
				String[]linesep = line.trim().split("\\s+");
				//System.out.println("DEBUG: "+ linesep[0] + Double.parseDouble(linesep[1]));
				aatable.put(linesep[0], Double.parseDouble(linesep[1]));
			}
		}
	}
	
	static void Makekmertable() throws NumberFormatException, IOException
	{
		File kmerfreqfile = new File("kmerFreqTreated.txt");
		BufferedReader br = new BufferedReader(new FileReader(kmerfreqfile));
		String line;
		
		kmertable = new Hashtable<String, Double>();
		
		while((line = br.readLine()) !=null)
		{
			if(line.trim().length()>0)
			{
				String[]linesep = line.trim().split("\\s+");
				if(linesep[0].equals("Total"))
				{
					totalkmer = Long.parseLong(linesep[1]);
				}
				else
				{
					kmertable.put(linesep[0], (Double)(Double.parseDouble(linesep[1])/totalkmer));
					//System.out.println("DEBUG: " + (Double)(Double.parseDouble(linesep[1])/totalkmer));
				}
			}
			
		}
	}
	
	static double scoreProtein(String PPline)
	{
		int linelen = PPline.length();
		double score = 0;
		for(int i=0; i<=linelen-5; i++)
		{
			String fourMer = PPline.substring(i, i+5);
			//System.out.println("DEBUG: " + fourMer);
			if(kmertable.get(fourMer) == null)
			{
				score += 0;
			}
			else
			{
				double probreal = kmertable.get(fourMer);
				double probgarb = aatable.get(fourMer.substring(0,1))*
						aatable.get(fourMer.substring(1,2))*
						aatable.get(fourMer.substring(2,3))*
						aatable.get(fourMer.substring(3,4))*
						aatable.get(fourMer.substring(4,5));
				double loglikelihood = Math.log(probreal/probgarb);
				score += loglikelihood;
				//System.out.println("DEBUG: "+score);
			}
		}
		return score;
	}
	
	static void DisplayResult(int index, String protein, double score, String type)
	{
		//String twodec = String.format("%.2f",score);
		System.out.println(""+index+"\t"+protein+"\t"+score+"\t"+type);
	}
	
	static void Printkmertable() throws FileNotFoundException, UnsupportedEncodingException
	{
		Enumeration fourMers;
		fourMers = kmertable.keys();
        //PrintWriter writer = new PrintWriter("kmerFreq.txt", "UTF-8");
        String mer;
        
        while(fourMers.hasMoreElements()) {
            mer = (String) fourMers.nextElement();
           //writer.println(mer+ " "+kmertable.get(mer));
            System.out.println(mer+ " "+kmertable.get(mer));
        }
        //writer.println("Total 3-mers: "+totalkmer);
        
        //writer.close();
	}
	
	public static void main(String[] args) throws IOException
	{
		/*if (args.length != 1) {
            System.out.println("Arguments: fastaInputFile");
            return;
		}*/
		
		Makeaatable();
		Makekmertable();
		//DEBUG:
		//Printaatable();
		//Printkmertable();
		
		File infile = new File("A4_training_set.txt");
		BufferedReader br = new BufferedReader(new FileReader(infile));
		PrintWriter writer = new PrintWriter("results.txt");
		
		String line;
		String type="";
		double proteinscore = 0;
		int proteincount = 1;
		while((line = br.readLine())!=null)
		{
			if(line.trim().length()>0) 		
			{
				if(line.substring(0, 1).equals(">")) 
				{
					type = line;
				}
				else
				{
					proteinscore = scoreProtein(line);
					DisplayResult(proteincount, line, proteinscore, type);
					writer.println(""+proteincount+"\t"+line+"\t"+proteinscore+"\t"+type);
					proteincount++;
				}
				//The totalloglikelihoodscore is used to compare the scoring heuristics of the 3-mer with regression testing of the average score
				//I consistently get that with a sample of anywhere between 5,000-100,000 amino acid sequences, the 
				//average loglikelihoodscore for all decoys is consistently around -0.15 while for all real proteins is consistently around +0.5
				//Showing some effectiveness in the scoring heuristics.
				
				totalloglikelihoodscore += proteinscore;
			}	
		}
		writer.close();
		//The average likelihoodratio
		System.out.println(""+totalloglikelihoodscore/20000);
	}
}
