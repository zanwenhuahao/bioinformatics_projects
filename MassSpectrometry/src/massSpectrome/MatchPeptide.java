package massSpectrome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Enumeration;
import java.util.Hashtable;


class trypepinfo {
	 public String pep;
	 public String protein;
	 public double mass;
	 
	 public trypepinfo(String pep, String protein, double mass)
	 {
		 this.pep = pep;
		 this.protein = protein;
		 this.mass = mass;
	 }
	 
	 public void print_info()
	 {
		 System.out.println(this.protein+" "+this.pep+" "+this.mass);
	 }
}

public class MatchPeptide {
	
	private static int index;
	private static double pepmass;
	private static int charge;
	private static int scans;
	private static double rtinseconds;
	
	private static double actual_mass;
	private static double tallest_intensity;
	private static String matchedpeptide;
	private static double peptopscore;
	
	private static Hashtable<Character, Double> aamasstable;
	private static Hashtable<String, trypepinfo> peptable;
	private static Hashtable<String, String> msms;
	
	private static int total_matched = 0;
	
	static void initialize_aamasstable()
	{
		aamasstable = new Hashtable<Character, Double>();
		aamasstable.put('A', 71.0371138);
		aamasstable.put('R', 156.1011110);
		aamasstable.put('N', 114.0429274);
		aamasstable.put('D', 115.0269430);
		aamasstable.put('C', 160.0306482);
		aamasstable.put('E', 129.0425931);
		aamasstable.put('Q', 128.0585775);
		aamasstable.put('G', 57.0214637);
		aamasstable.put('H', 137.0589119);
		aamasstable.put('I', 113.0840640);
		aamasstable.put('L', 113.0840640);
		aamasstable.put('K', 128.0949630);
		aamasstable.put('M', 131.0404846);
		aamasstable.put('F', 147.0684139);
		aamasstable.put('P', 97.0527639);
		aamasstable.put('S', 87.0320284);
		aamasstable.put('T', 101.0476785);
		aamasstable.put('W', 186.0793130);
		aamasstable.put('Y', 163.0633285);
		aamasstable.put('U', 168.9641990);
		aamasstable.put('V', 99.0684139);
	}
	
	static void initialize_peptable()
	{
		peptable = new Hashtable<String, trypepinfo>();
	}
	
	static double calculate_mass(String pep)
	{
		int pep_len = pep.length();
		double totalmass = 0;
		for(int i=0;i<pep_len;i++)
		{
			totalmass += aamasstable.get(pep.charAt(i));
		}
		totalmass += 18.01;
		return totalmass;
	}
	
	static void parse_to_trypep(String protein_name, String protein)
	{
		int line_len = protein.length();
		int prev_break = 0;
		String temp_peptide = "";
		double temp_pep_mass = 0;
		
		for(int i=0;i<line_len-1;i++)
		{
			if((protein.charAt(i) == 'K' || protein.charAt(i) == 'R') && protein.charAt(i+1) !='P')
			{
				temp_peptide = protein.substring(prev_break, i+1);
				temp_pep_mass = calculate_mass(temp_peptide);
				trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
				peptable.put(temp_peptide, pep);
				prev_break = i+1;
			}
		}
		
		if(prev_break!=line_len)
		{
			temp_peptide = protein.substring(prev_break, line_len);
			temp_pep_mass = calculate_mass(temp_peptide);
			trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
			peptable.put(temp_peptide, pep);
			prev_break = line_len;
		}
		
		prev_break = 0;
		int one_miss = 0;
		int middle_break = 0;
		int counter = 0;
		
		while(counter < line_len-1)
		{
			if((protein.charAt(counter) == 'K' || protein.charAt(counter) == 'R')  && one_miss == 0)
			{
				one_miss = 1;
				middle_break=counter+1;
				counter++;
			}
			else if((protein.charAt(counter) == 'K' || protein.charAt(counter) == 'R') 
					&& protein.charAt(counter+1) !='P' && one_miss == 1)
			{
				temp_peptide = protein.substring(prev_break, counter+1);
				temp_pep_mass = calculate_mass(temp_peptide);
				trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
				peptable.put(temp_peptide, pep);
				counter = middle_break;
				prev_break = counter;
				one_miss = 0;
			}
			else
			{
				counter++;
			}
		}
		
		if(prev_break!=line_len)
		{
			temp_peptide = protein.substring(prev_break, line_len);
			temp_pep_mass = calculate_mass(temp_peptide);
			trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
			peptable.put(temp_peptide, pep);
			prev_break = line_len;
		}
		
		prev_break = 0;
		counter = 0;
		int two_miss = 0;
		int middle_break_one = 0;
		//int middle_break_two = 0;
		
		while(counter < line_len-1)
		{
			if((protein.charAt(counter) == 'K' || protein.charAt(counter) == 'R') 
					&& two_miss == 0)
			{
				two_miss++;
				middle_break_one=counter+1;
				counter++;
			}
			else if((protein.charAt(counter) == 'K' || protein.charAt(counter) == 'R') 
					&& two_miss == 1)
			{
				two_miss++;
				//middle_break_two=counter+1;
				counter++;
			}
			else if((protein.charAt(counter) == 'K' || protein.charAt(counter) == 'R') 
					&& two_miss == 2)
			{
				temp_peptide = protein.substring(prev_break, counter+1);
				temp_pep_mass = calculate_mass(temp_peptide);
				trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
				peptable.put(temp_peptide, pep);
				counter = middle_break_one;
				prev_break = counter;
				two_miss = 0;
			}
			else
			{
				counter++;
			}
		}
		if(prev_break!=line_len)
		{
			temp_peptide = protein.substring(prev_break, line_len);
			temp_pep_mass = calculate_mass(temp_peptide);
			trypepinfo pep = new trypepinfo(temp_peptide, protein_name, temp_pep_mass);
			peptable.put(temp_peptide, pep);
			prev_break = line_len;
		}
		
			
	}
	
	static void initialize_ms_table()
	{
		msms = new Hashtable<String, String>();
	}
	
	static void ms_array_maker(String line, int nth)
	{
		String []line_vector = line.split("\\s+");
		//System.out.println(line_vector[0]+" "+line_vector[1]);
		msms.put(line_vector[0], line_vector[1]);
	}
	
	static void calculate_actual_mass()
	{
		actual_mass = pepmass*charge-1.00728*charge;
	}
	
	static void find_tallest_peak()
	{
		Enumeration spectra_masses;
		spectra_masses = msms.keys();
		double most_intense = 0;
		
		String a_mass;
		double a_intensity;
		while(spectra_masses.hasMoreElements())
		{
			a_mass = (String) spectra_masses.nextElement();
			a_intensity = Double.parseDouble(msms.get(a_mass));
			
			if(a_intensity > most_intense)
			{
				most_intense = a_intensity;
			}
		}
		tallest_intensity = most_intense;
	}
	
	static double score_y_ion(String yion)
	{
		double yion_mass = calculate_mass(yion) - 18.01 + 19.0178;
		//System.out.println("This yion has mass:" + yion_mass);
		double tallest_peak = 0;
		double intensity_tallest_peak = 0;
		double final_y_ion_score = 0;
		
		Enumeration spectra_masses;
		spectra_masses = msms.keys();
		
		String a_mass;
		double mass_num;
		while(spectra_masses.hasMoreElements())
		{
			a_mass = (String) spectra_masses.nextElement();
			mass_num = Double.parseDouble(a_mass);
			//System.out.println("A spectra's mass: " + mass_num);
			if(mass_num>=yion_mass-0.5 && mass_num<=yion_mass+0.5)
			{
				double this_peak = mass_num;
				double intensity_this_peak = Double.parseDouble(msms.get(a_mass));
				
				if(intensity_this_peak > intensity_tallest_peak)
				{
					intensity_tallest_peak = intensity_this_peak;
					tallest_peak = this_peak;
				}
			}
		}
		double ratio_x = intensity_tallest_peak/tallest_intensity;
		final_y_ion_score = Math.max(0,Math.log10(100*ratio_x));
		
		return final_y_ion_score;
	}
	
	static double score_function(String a_pep)
	{
		//System.out.println("Trying pep: "+a_pep);
		double total_score = 0;
		String yion = "";
		int pep_len = a_pep.length();
		
		for(int i=1; i<pep_len;i++)
		{
			yion = a_pep.substring(i,pep_len);
			//System.out.println("The Yion is: "+yion);
			total_score += score_y_ion(yion);
		}
		//System.out.println("Total sum of y-ion score:" + total_score);
		return total_score;
	}
	
	static void find_matching_peptide()
	{
		String best_match_pep = "";
		double best_match_score = 0;
		
		Enumeration all_peptides;
		all_peptides = peptable.keys();
		
		String a_pep;
		while(all_peptides.hasMoreElements())
		{
			a_pep = (String) all_peptides.nextElement();
			
			trypepinfo a_pep_info = peptable.get(a_pep);
			if(a_pep_info.mass>=actual_mass-0.1 && a_pep_info.mass<=actual_mass+0.1)
			{	
				//System.out.println("a_pep within range: "+a_pep);
				double match_score = score_function(a_pep);
				if(match_score>best_match_score)
				{
					best_match_score = match_score;
					best_match_pep = a_pep;
				}
			}
		}
		matchedpeptide = best_match_pep;
		//System.out.println("matched_peptide: " + matchedpeptide);
		peptopscore = best_match_score;
	}
	
	static void print_trypep_table()
	{
		Enumeration trypeps;
		trypeps = peptable.keys();
        
		String the_pep;
		
        while(trypeps.hasMoreElements()) {
            the_pep = (String) trypeps.nextElement();
            trypepinfo pep_to_print = peptable.get(the_pep);
            pep_to_print.print_info();
        }
	}
	
	static void print_msms_info()
	{
		System.out.println("Index: "+index);
		System.out.println("Pepmass: "+pepmass);
		System.out.println("Charge: "+charge);
		System.out.println("Scans: "+scans);
		System.out.println("Rtinseconds: "+rtinseconds);
		System.out.println("Tallest Peak Intensity: "+tallest_intensity);
		System.out.println("Actual mass: "+actual_mass);
		Enumeration table_keys;
		table_keys = msms.keys();
        
		String the_mass;
		
        while(table_keys.hasMoreElements()) {
            the_mass = (String) table_keys.nextElement();
            String intensity = msms.get(the_mass);
            System.out.println(the_mass+" "+intensity);
        }

	}
	
	static void print_match_result()
	{
		String protein = "";
		if (!(matchedpeptide.equals("")))
		{
			trypepinfo the_matched_pep_info = peptable.get(matchedpeptide);
			protein = the_matched_pep_info.protein;	
			total_matched++;
		}
		
		
		System.out.println(index + "\t" + pepmass + "\t" + charge + "\t" +
							matchedpeptide + "\t" + protein + "\t" + peptopscore);
	}
	
	public static void main(String[] args) throws IOException {
			/*if (args.length != 2) {
	        System.out.println("Arguments: fastaFile nSamples [label|nolabel]");
	        System.out.println(" label will label the target and decoy in header lines");
	        System.out.println(" nolabel will not label the target and decoy in header lines");
	        return;
	    	} */
		
		initialize_aamasstable();
		
		//Store the fasta file's protein information in the tryptic peptide hashtable
		
		
	    File fastaFile = new File("test_fasta.txt");
	    BufferedReader br = new BufferedReader(new FileReader(fastaFile));
		PrintWriter writer = new PrintWriter("results.txt");
		
		String line;
		
		boolean read_protein = false;
		
		initialize_peptable();
		
		String protein_name = "";
		String peptide_name = "";
		double peptide_mass = 0;
		while((line = br.readLine())!=null)
		{
			if((line.trim().length() > 0) && (line.trim().substring(0, 1).equals(">"))) 
				//Name of the protein
			{
				protein_name = line.substring(0,10);
				read_protein = true;
			}
			
			if((line.trim().length()>0) 
					&& !(line.substring(0, 1).equals(">")) && read_protein) 
				//The protein line
			{
				//System.out.println(protein_name);
				parse_to_trypep(protein_name, line);
				read_protein = false;
			}
		}
		
		//print_trypep_table();
		//Process the ms/ms file
		
		File mgfFile = new File("test_mgf.txt");
	    br = new BufferedReader(new FileReader(mgfFile));
	    
	    int nth=0;
	    
		while((line = br.readLine())!= null)
		{
			line=line.trim();
			String[]linesep = line.split("\\s+");
			if(line.length()>0 && linesep[0].equals("END"))
			{
				find_tallest_peak();
				//print_msms_info();
				find_matching_peptide();
				
				print_match_result();
				writer.println(index + "\t" + pepmass + "\t" + charge + "\t" +
						matchedpeptide + "\t" + peptopscore);
				nth=0;
			}
			
			else if(line.length()==0 || linesep[0].equals("BEGIN"))
			{
				initialize_ms_table();
				continue;
			}
			
			else
			{
				linesep = line.split("=");
				if(linesep[0].equals("TITLE"))
				{
					index = Integer.parseInt(linesep[2]);
				}
				else if(linesep[0].equals("PEPMASS"))
				{
					pepmass = Double.parseDouble(linesep[1]);
				}
				else if(linesep[0].equals("CHARGE"))
				{
					String[] just_charge = linesep[1].split("\\+");
					charge = Integer.parseInt(just_charge[0]);
				}
				else if(linesep[0].equals("SCANS"))
				{
					scans = Integer.parseInt(linesep[1]);
				}
				else if(linesep[0].equals("RTINSECONDS"))
				{
					rtinseconds = Double.parseDouble(linesep[1]);
					calculate_actual_mass();
				}
				else
				{
						ms_array_maker(line, nth);
						nth++;
				}
			}
		}
		
		System.out.println("Total Matched: "+total_matched);
		writer.close();
		//System.out.println("Program finished.");
	}

}
