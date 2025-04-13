package org.example;

import org.example.algorithms.SlidingWindow;
import org.example.models.AnalysisResult;
import org.example.utils.FastaParser;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.*;

import static java.lang.Math.min;

public class Main {
    public static void main(String[] args) throws Exception {
        if(args.length < 4)
            throw new Exception("Program accepts at least 4 arguments:" +
                    "<sequence type> <recombinant file> <parent 1 file> <parent 2 file> <number of threads - optional> <window size - optional> <step size - optional>\n" +
                    "Accepted sequence types: d - DNA | r - RNA\n" +
                    "Accepted file formats: FASTA");

        long start = System.currentTimeMillis();

        String recombinant = FastaParser.readFastaSequence(args[1], args[0]);
        String parent1 = FastaParser.readFastaSequence(args[2], args[0]);
        String parent2 = FastaParser.readFastaSequence(args[3], args[0]);

        int thCount = 2;
        int windowSize = 20;
        int stepSize = 5;

        if(args.length > 4){
            try{
                thCount = Integer.parseInt(args[4]);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid input for 'number of threads' : " + args[4]);
            }

            if(args.length > 5){
                try{
                    windowSize = Integer.parseInt(args[5]);
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Invalid input for 'window size' : " + args[5]);
                }

                if(args.length > 6){
                    try{
                        stepSize = Integer.parseInt(args[6]);
                    } catch (NumberFormatException e) {
                        throw new IllegalArgumentException("Invalid input for 'step size' : " + args[4]);
                    }
                }
            }
        }

        int length = recombinant.length();
        int chunkSize = (length + thCount - 1) / thCount;

        ExecutorService executor = Executors.newFixedThreadPool(thCount);

        List<Future<List<AnalysisResult>>> futures = new ArrayList<>();

        for(int chunkStartInd = 0; chunkStartInd < length; chunkStartInd += chunkSize){
            int chunkStart = chunkStartInd;
            int chunkEnd = min(chunkStart + chunkSize + windowSize, length);

            String recombinantChunk = recombinant.substring(chunkStart, chunkEnd);

            int finalWindowSize = windowSize;
            int finalStepSize = stepSize;
            futures.add(executor.submit(
                    () -> SlidingWindow.analyzeRecombinantWithParents(
                            recombinantChunk,
                            parent1, parent2,
                            args[0],
                            finalWindowSize, finalStepSize,
                            chunkStart, chunkEnd)
            ));
        }

        executor.shutdown();
        try{
            if(!executor.awaitTermination(60, TimeUnit.SECONDS))
                executor.shutdownNow();
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }

        List<AnalysisResult> finalResults = new ArrayList<>();

        for (Future<List<AnalysisResult>> future : futures) {
            finalResults.addAll(future.get());
        }

        finalResults.sort(Comparator.comparingInt(AnalysisResult::segmentStartLocation));

        for (AnalysisResult result : finalResults) {
            System.out.println(result.toString());
        }

        long end = System.currentTimeMillis();
        long elapsed = end - start;
        System.out.println("Time elapsed to complete task: " + elapsed + "ms.");
    }
}