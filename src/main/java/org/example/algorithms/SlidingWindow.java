package org.example.algorithms;

import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.alignment.Alignments;
import org.example.models.AnalysisResult;
import org.example.utils.enums.RecombinationType;
import org.example.utils.enums.SimilarTo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static java.lang.Math.abs;
import static java.lang.Math.min;

public class SlidingWindow {
    private static final double THRESHOLD = 0.1;

    public static List<AnalysisResult> analyzeRecombinantWithParents(
            String recombinantChunk,
            String parent1, String parent2,
            String sequenceType,
            int window, int step,
            int chunkStartPosition, int fullRecombinantLastPosition) throws Exception {
        if(window <= 0 || step <= 0)
            throw new Exception("Window and step sizes should be greater than zero. Given window size: '" + window + "' and step size: '" + step + "' ..");

        List<AnalysisResult> results = Collections.synchronizedList(new ArrayList<>());

        int length = recombinantChunk.length();

        for(int start = 0; start <= length - window; start += step){
            int end = min((start + window), length);

            int windowStartPosition = chunkStartPosition + start;
            int windowEndPosition = min(chunkStartPosition + end, fullRecombinantLastPosition);

            String recombinantWindow = recombinantChunk.substring(start, end);
            results.add(analyzeWindow(recombinantWindow, parent1, parent2, sequenceType, results, windowStartPosition, windowEndPosition));
        }

        return results;
    }

    private static AnalysisResult analyzeWindow(
        String recombinantWindow,
        String parent1, String parent2,
        String sequenceType,
        List<AnalysisResult> results,
        int windowStartPosition, int windowEndPosition
    )throws Exception
    {
        double similarityWithP1 = SmithWaterman.localAlign(recombinantWindow, parent1, sequenceType);
        double similarityWithP2 = SmithWaterman.localAlign(recombinantWindow, parent2, sequenceType);

        SimilarTo moreSimilarTo;
        if(similarityWithP1 > similarityWithP2) moreSimilarTo = SimilarTo.FirstParent;
        else if (similarityWithP1 < similarityWithP2) moreSimilarTo = SimilarTo.SecondParent;
        else moreSimilarTo = SimilarTo.Equal;

        boolean majorChange = abs(similarityWithP1 - similarityWithP2) >= THRESHOLD;
        
        RecombinationType recombinationType = analyzeStates(results, majorChange, moreSimilarTo, similarityWithP1, similarityWithP2);
        
        return new AnalysisResult(
                windowStartPosition, windowEndPosition,
                moreSimilarTo,
                similarityWithP1, similarityWithP2,
                majorChange,
                recombinationType);
    }
    
    private static RecombinationType analyzeStates(List<AnalysisResult> results, boolean majorChange, SimilarTo currentlySimilarTo,
                                                   double similarityWithP1, double similarityWithP2) {

        if(results.isEmpty()) return RecombinationType.None;

        AnalysisResult last = results.getLast();

        boolean flippedDirectly = currentlySimilarTo != SimilarTo.Equal
                && last.moreSimilarTo() != SimilarTo.Equal
                && currentlySimilarTo != last.moreSimilarTo();

        if(flippedDirectly && majorChange) return RecombinationType.Hard;

        boolean fromEqualToParentWithDrop = currentlySimilarTo != SimilarTo.Equal
                && last.moreSimilarTo() == SimilarTo.Equal
                && (similarityWithP1 <= 0.9 || similarityWithP2 <= 0.9);

        if(fromEqualToParentWithDrop && majorChange) return RecombinationType.Ambiguous;

        if(results.size() > 2){
            AnalysisResult secondLast = results.get(results.size() - 2);

            boolean streakOfTwoBroken = secondLast.moreSimilarTo() == last.moreSimilarTo() && flippedDirectly;

            if(streakOfTwoBroken) return RecombinationType.Ambiguous;

            boolean flippedAfterEqual = secondLast.moreSimilarTo() != SimilarTo.Equal
                    && last.moreSimilarTo() == SimilarTo.Equal
                    && currentlySimilarTo != SimilarTo.Equal
                    && secondLast.moreSimilarTo() != currentlySimilarTo;

            if(flippedAfterEqual && majorChange) return RecombinationType.Hard;
            if(flippedAfterEqual) return RecombinationType.Ambiguous;

            if(results.size() > 3){
                AnalysisResult thirdLast = results.get(results.size() - 3);

                boolean increasingP1 =
                        thirdLast.similarityWithParent1() < secondLast.similarityWithParent1() &&
                                secondLast.similarityWithParent1() < last.similarityWithParent1() &&
                                last.similarityWithParent1() < similarityWithP1;

                boolean decreasingP2 =
                        thirdLast.similarityWithParent2() > secondLast.similarityWithParent2() &&
                                secondLast.similarityWithParent2() > last.similarityWithParent2() &&
                                last.similarityWithParent2() > similarityWithP2;

                boolean decreasingP1 =
                        thirdLast.similarityWithParent1() > secondLast.similarityWithParent1() &&
                                secondLast.similarityWithParent1() > last.similarityWithParent1() &&
                                last.similarityWithParent1() > similarityWithP1;

                boolean increasingP2 =
                        thirdLast.similarityWithParent2() < secondLast.similarityWithParent2() &&
                                secondLast.similarityWithParent2() < last.similarityWithParent2() &&
                                last.similarityWithParent2() < similarityWithP2;

                boolean similaritySwitching = (increasingP1 && decreasingP2) || (increasingP2 && decreasingP1);

                if(similaritySwitching && majorChange) return RecombinationType.Hard;

                boolean suddenChangeFromOneParentToAnother = thirdLast.moreSimilarTo() == secondLast.moreSimilarTo()
                        && secondLast.moreSimilarTo() == last.moreSimilarTo()
                        && last.moreSimilarTo() != currentlySimilarTo
                        && last.moreSimilarTo() != SimilarTo.Equal
                        && currentlySimilarTo != SimilarTo.Equal;

                if(suddenChangeFromOneParentToAnother) return RecombinationType.Hard;
            }
        }

        if(flippedDirectly) return RecombinationType.Ambiguous;

        return RecombinationType.None;
    }
}
