package org.example.utils;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

import java.io.*;
import java.util.Map;

public class FastaParser {
    public static String readFastaSequence(String filepath, String sequenceTypeExpected) throws Exception{
        return switch (sequenceTypeExpected) {
            case "d" -> readFastaDNASequence(filepath);
            case "r" -> readFastaRNASequence(filepath);
            default -> throw new Exception("Unknown sequence type: '" + sequenceTypeExpected + "' .. ");
        };
    }

    private static String readFastaDNASequence(String filepath) throws Exception{
        if(!isFastaFile(filepath))
            throw new Exception("No file found at filepath: " + filepath);
        var dnaFile = readFile(filepath);

        try{
            Map<String, DNASequence> sequences = FastaReaderHelper.readFastaDNASequence(dnaFile);

            for(var key : sequences.keySet()){
                var seq = sequences.get(key);
                return seq.getSequenceAsString();
            }
        }
        catch (IOException e){
            System.err.println("Error while reading file: " + filepath + " .. " + e);
        }
        catch (RuntimeException e){
            System.err.println("Wrong type of sequence provided with file: " + filepath + " expected type is DNA .. " + e);
        }

        return "";
    }

    private static String readFastaRNASequence(String filepath) throws Exception{
        if(!isFastaFile(filepath))
            throw new Exception("File format should be in FASTA: " + filepath);
        var rnaFile = readFile(filepath);

        try{
            Map<String, RNASequence> sequences = FastaReaderHelper.readFastaRNASequence(rnaFile);

            for(var key : sequences.keySet()){
                var seq = sequences.get(key);
                return seq.getSequenceAsString();
            }
        }
        catch (IOException e){
            System.err.println("Error while reading file: " + filepath + " .. " + e);
        }
        catch (RuntimeException e){
            System.err.println("Wrong type of sequence provided with file: " + filepath + " expected type is RNA .. " + e);
        }

        return "";
    }

    private static File readFile(String filepath) throws IOException {
        File file = new File(filepath);
        if(!file.exists())
            throw new IOException("File does not exist: " + filepath);

        return file;
    }
    private static Boolean isFastaFile(String filepath) {
        String name = filepath.toLowerCase();

        if(!name.endsWith(".fasta") &&
            !name.endsWith(".fas") &&
            !name.endsWith(".fa") &&
            !name.endsWith(".fna") &&
            !name.endsWith(".ffn") &&
            !name.endsWith(".faa") &&
            !name.endsWith(".mpfa") &&
            !name.endsWith(".frn")) return false;

        try{
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            String firstLine = br.readLine();
            if(firstLine == null || !firstLine.startsWith(">")){
                br.close();
                return false;
            }
            br.close();
        }
        catch (IOException e){
            System.err.println("Error reading file: " + filepath + " .. " + e);
        }

        return true;
    }
}
