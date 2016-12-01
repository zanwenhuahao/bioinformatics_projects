
import java.io.*;
import java.util.Enumeration;
import java.util.Hashtable;

public class ProteinFrequency {
	
	private static long[] AAarray = new long[26];
	private static long total = 0;
	private static BufferedReader br;
	private static Hashtable<String, Integer> kmertable;
	private static long totalkmer = 0;
	private static String AAvalid = "ACDEFGHIKLMNPQRSTVWY"; 
	
	static void FreqCalc (String AAline)
	{
		long AAlen=AAline.length();
		for(int i=0;i<AAlen;i++)
		{
			int asciinum = (int)AAline.charAt(i); 
			if(asciinum-65>=0 && asciinum-65<=25)
			{
				AAarray[asciinum-65]++;
			}
			if(asciinum-97>=0 && asciinum-97<=25)
			{
				AAarray[asciinum-97]++;
			}
		}
	}
	
	static void kmerHashBuilder()
	{
		kmertable = new Hashtable<String, Integer>();
		for(int i=0; i<20; i++)
		{
			for(int j=0; j<20; j++)
			{
				for(int k=0; k<20; k++)
				{
					for(int l=0; l<20; l++)
					{
						for(int m=0; m<20;m++)
						{
							kmertable.put(""+AAvalid.substring(i,i+1)+AAvalid.substring(j,j+1)+AAvalid.substring(k,k+1)+AAvalid.substring(l,l+1)+AAvalid.substring(m,m+1), 0);
						}
					}
				}
			}
		}
	}
	
	static void KmerFreq(String PPline)
	{
		int linelen = PPline.length();
		for(int i=0; i<=linelen-5; i++)
		{
			String fourMer = PPline.substring(i, i+5);
			//System.out.println("DEBUG: " + threeMer);
			int curoccur = 0;
			if(kmertable.get(fourMer) != null)
			{
				curoccur = (Integer) kmertable.get(fourMer);
			}
			kmertable.put(fourMer, curoccur+1);
			totalkmer++;
		}
	}
	
	static void DisplayFreq()
	{
		for(int i=0; i<AAarray.length; i++)
		{
			if(i!=1 && i!=9 && i!=14 && i!=20 && i!=23 && i!=25)
			{
				total+=AAarray[i];
			}
		}
		
		for(int i=0; i<AAarray.length; i++)
		{
			if(i!=1 && i!=9 && i!=14 && i!=20 && i!=23 && i!=25)
			{
				double eachfreq;
				if (total == 0)
				{	
					eachfreq = (double) AAarray[i];
				}
				else 
				{
					eachfreq = (double) AAarray[i]/total;
				}
				String fourdec = String.format("%.4f", eachfreq);
				System.out.println((char)(i+65) + " " + fourdec);
			}
		}
	}
	
	static void WriteKmerFreq() throws FileNotFoundException, UnsupportedEncodingException
	{
		Enumeration fourMers;
		fourMers = kmertable.keys();
        PrintWriter writer = new PrintWriter("kmerFreq.txt", "UTF-8");
        String mer;
        
        while(fourMers.hasMoreElements()) {
            mer = (String) fourMers.nextElement();
            writer.println(mer+ " "+kmertable.get(mer));
            //System.out.println(mer+ " "+kmertable.get(mer));
        }
        writer.println("Total 4-mers: "+totalkmer);
        
        writer.close();
	}
	
	public static void main(String[] args) throws IOException 
	{
		/*File infile = new File("uniprot_sprot.txt");
		br = new BufferedReader(new FileReader(infile));
		String line;
		while((line = br.readLine())!=null)
		{
			if((line.trim().length()>0) 
					&& !(line.substring(0, 1).equals(">")) 
					&& !(line.substring(0,1).equals(";"))) //If it's not a header line
			{
				FreqCalc(line);
			}
		}	
		DisplayFreq();*/
		
		//File databank = new File("real_protein_training_set.txt");
		File databank = new File("uniprot_sprot.txt");
		br = new BufferedReader(new FileReader(databank));
		String line;
		
		kmerHashBuilder();
		
		while((line = br.readLine())!=null)
		{
			if(line.trim().length()>4 && !(line.substring(0, 1).equals(">")))
			{
				KmerFreq(line);
			}
			
			if((line.trim().length()>0) 
					&& !(line.substring(0, 1).equals(">")) 
					&& !(line.substring(0,1).equals(";"))) //If it's not a header line
			{
				FreqCalc(line);
			}
		}
		
		WriteKmerFreq();
		DisplayFreq();
	}

}

