import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class GetRandomTrainingSet {
    List<String> proteins;
    int[] accumLen;
    Random rand;
    boolean [] isAminoAcid;

    GetRandomTrainingSet(File file) throws IOException {
        //read sequences
        BufferedReader reader = new BufferedReader(new FileReader(file));
        proteins = new ArrayList<String>();
        String line;
        StringBuilder b = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.startsWith(">")) {
                if (b.length() > 0) {
                    proteins.add(b.toString());
                    b = new StringBuilder();
                }
            } else {
                b.append(line);
            }
        }
        if (b.length() > 0) {
            proteins.add(b.toString());
        }
        reader.close();

        //compute accumlen
        accumLen = new int[proteins.size() + 1];
        accumLen[0] = 0;
        for (int i = 0; i < proteins.size(); i++) {
            accumLen[i + 1] = accumLen[i] + proteins.get(i).length();
        }

        rand = new Random();

        //determine if a letter is an amino acid
        isAminoAcid = new boolean[26];
        String aa = "ARNDCEQGHILKMFPSTWYV";
        for (int i = 0; i < aa.length(); i++) {
            char ch = aa.charAt(i);
            isAminoAcid[ch-'A'] = true;
        }
    }

    static class Sample {
        boolean isTarget;
        String peptide;

        public Sample(boolean isTarget, String peptide) {
            this.isTarget = isTarget;
            this.peptide = peptide;
        }
    }

    static String shuffle(String str) {
        List<Character> list = new ArrayList<Character>();
        for (int i = 0; i < str.length(); i++) {
            list.add(str.charAt(i));
        }
        Collections.shuffle(list);
        StringBuilder b = new StringBuilder();
        for (Character ch : list) {
            b.append(ch);
        }
        return b.toString();
    }

    boolean isPeptide(String peptide){
        for (int i = 0; i < peptide.length(); i++) {
            if (!isAminoAcid[peptide.charAt(i)-'A']){
                return false;
            }
        }
        return true;
    }

    //create a random peptide with peptide length.
    String randomPeptide(int peptideLength) {
        while(true){
            int proteinId;
            int peptideStart;
            do {
                int randPos = rand.nextInt(accumLen[accumLen.length - 1]);
                proteinId = Arrays.binarySearch(accumLen, randPos);
                if (proteinId < 0) {
                    proteinId = -proteinId - 2;
                }
                peptideStart = randPos - accumLen[proteinId];
                assert peptideStart >= 0 && peptideStart < proteins.get(proteinId).length();
            } while (proteins.get(proteinId).length() < peptideStart + peptideLength);
            String peptide = proteins.get(proteinId).substring(peptideStart, peptideStart + peptideLength);
            if (isPeptide(peptide)){
                return peptide;
            } else {
                System.out.println(peptide);
            }
        }
    }

    //produce n pairs of samples
    List<Sample> randomSamples(int nSamples){
        List<Sample> sampleList = new ArrayList<Sample>();
        for (int i = 0; i < nSamples; i++) {
            int peptideLength = rand.nextInt(31)+20;
            String target = randomPeptide(peptideLength);
            //No decoys, only real ones
            //String decoy = randomPeptide(peptideLength);
            //decoy = shuffle(decoy);
            sampleList.add(new Sample(true, target));
            //sampleList.add(new Sample(false, decoy));
        }
        return sampleList;
    }

    public static void main(String[] args) throws IOException {
      /*  if (args.length != 3) {
            System.out.println("Arguments: fastaFile nSamples [label|nolabel]");
            System.out.println(" label will label the target and decoy in header lines");
            System.out.println(" nolabel will not label the target and decoy in header lines");
            return;
        } */

        File fastaFile = new File("uniprot_sprot.txt");
        int nSamples = 100000;
        //boolean label = true;
        
        //int nSamples = Integer.parseInt(args[1]);
        //boolean label = args[2].equals("label");


        GetRandomTrainingSet sampler = new GetRandomTrainingSet(fastaFile);
        List<Sample> sampleList = sampler.randomSamples(nSamples);
        Collections.shuffle(sampleList);
        PrintWriter writer = new PrintWriter("real_protein_training_set.txt", "UTF-8");
        for (int i = 0; i < sampleList.size(); i++) {
            Sample sample = sampleList.get(i);
            /*String header;
            if (label){
                header = sample.isTarget ? ">target" : ">decoy";
            } else {
                header = ">peptide";
            }
            writer.println(header);*/
            writer.println(sample.peptide);
        }
        writer.close();
    }
}
