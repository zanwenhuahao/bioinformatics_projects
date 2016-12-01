import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class SamplerA4 {
    List<String> targetProteins;
    List<String> decoyProteins;
    int[] accumLen;
    Random rand;
    boolean [] isAminoAcid;

    SamplerA4(File file) throws IOException {
        //read sequences
        BufferedReader reader = new BufferedReader(new FileReader(file));
        targetProteins = new ArrayList<String>();
        String line;
        StringBuilder b = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.startsWith(">")) {
                if (b.length() > 0) {
                    targetProteins.add(b.toString());
                    b = new StringBuilder();
                }
            } else {
                b.append(line);
            }
        }
        if (b.length() > 0) {
            targetProteins.add(b.toString());
        }
        reader.close();

        //compute accumlen
        accumLen = new int[targetProteins.size() + 1];
        accumLen[0] = 0;
        for (int i = 0; i < targetProteins.size(); i++) {
            accumLen[i + 1] = accumLen[i] + targetProteins.get(i).length();
        }

        rand = new Random();

        //determine if a letter is an amino acid
        isAminoAcid = new boolean[26];
        String aa = "ARNDCEQGHILKMFPSTWYV";
        for (int i = 0; i < aa.length(); i++) {
            char ch = aa.charAt(i);
            isAminoAcid[ch-'A'] = true;
        }

        //create decoys
        decoyProteins = new ArrayList<String>();
        for (int i = 0; i < targetProteins.size(); i++) {
            decoyProteins.add(shuffle(targetProteins.get(i)));
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
    String randomPeptide(int peptideLength, boolean isTarget) {
        List<String> proteins = isTarget ? this.targetProteins : this.decoyProteins;
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
//                System.out.println(peptide);
            }
        }
    }

    //produce n pairs of samples
    List<Sample> randomTargetSamples(int nSamples){
        List<Sample> sampleList = new ArrayList<Sample>();
        for (int i = 0; i < nSamples; i++) {
            int peptideLength = rand.nextInt(31)+20;
            String target = randomPeptide(peptideLength, true);
            String decoy = randomPeptide(peptideLength, false);
            sampleList.add(new Sample(true, target));
            sampleList.add(new Sample(false, decoy));
        }
        return sampleList;
    }

    public static void main(String[] args) throws IOException {
        /*if (args.length != 3) {
            System.out.println("Arguments: fastaFile nSamples [label|nolabel]");
            System.out.println(" label will label the target and decoy in header lines");
            System.out.println(" nolabel will not label the target and decoy in header lines");
            return;
        }*/

        File fastaFile = new File("uniprot_sprot.txt");
        int nSamples = 50000;
        boolean label = true;


        SamplerA4 sampler = new SamplerA4(fastaFile);
        List<Sample> sampleList = sampler.randomTargetSamples(nSamples);
        Collections.shuffle(sampleList);
        
        PrintWriter writer = new PrintWriter("A4_training_set.txt", "UTF-8");
        for (int i = 0; i < sampleList.size(); i++) {
            Sample sample = sampleList.get(i);
            String header;
            if (label){
                header = sample.isTarget ? ">target" : ">decoy";
            } else {
                header = ">peptide";
            }
            writer.println(header);
            writer.println(sample.peptide);
        }
        writer.close();
    }
}
