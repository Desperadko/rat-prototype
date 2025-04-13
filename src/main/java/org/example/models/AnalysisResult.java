package org.example.models;

import org.example.utils.enums.RecombinationType;
import org.example.utils.enums.SimilarTo;

public record AnalysisResult(
        int segmentStartLocation,
        int segmentEndLocation,
        SimilarTo moreSimilarTo,
        double similarityWithParent1,
        double similarityWithParent2,
        boolean majorChange,
        RecombinationType recombinationType
) {
    @Override
    public String toString(){
        String similarTo = switch (moreSimilarTo){
            case SimilarTo.FirstParent -> "more similar to the First parent";
            case SimilarTo.SecondParent -> "more similar to the Second parent";
            case SimilarTo.Equal -> "Equal to both";
        };
        return "Analysis result for [" + segmentStartLocation() + " - " + segmentEndLocation() + "]:" +
                " Segment is - " + similarTo +
                ", with similarities -> First parent: " + similarityWithParent1() * 100 + "% | Second parent: " + similarityWithParent2() * 100 + "%." +
                " Has a recombination happened: " + (recombinationType == RecombinationType.None ? "No" : (recombinationType == RecombinationType.Hard ? "Yes" : "Possibly")) + ".";
    }
}
