package org.example.algorithms;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.RNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

public class SmithWaterman {
    public static double localAlign(String window, String parent, String sequenceType) throws Exception {
        return switch (sequenceType){
            case "d" -> localAlignDNA(window, parent);
            case "r" -> localAlignRNA(window, parent);
            default -> throw new Exception("Unknown sequence type: '" + sequenceType + "' ..");
        };
    }

    private static double localAlignDNA(String window, String parent) throws CompoundNotFoundException {
        DNASequence windowSeq = new DNASequence(window);
        DNASequence parentSeq = new DNASequence(parent);

        SequencePair<DNASequence, NucleotideCompound> pair = Alignments.getPairwiseAlignment(
                windowSeq,
                parentSeq,
                Alignments.PairwiseSequenceAlignerType.LOCAL,
                new SimpleGapPenalty(),
                SubstitutionMatrixHelper.getNuc4_4());

        return pair.getPercentageOfIdentity(true);
    }

    private static double localAlignRNA(String window, String parent) throws CompoundNotFoundException {
        RNASequence windowSeq = new RNASequence(window);
        RNASequence parentSeq = new RNASequence(parent);

        SequencePair<RNASequence, NucleotideCompound> pair = Alignments.getPairwiseAlignment(
                windowSeq,
                parentSeq,
                Alignments.PairwiseSequenceAlignerType.LOCAL,
                new SimpleGapPenalty(),
                SubstitutionMatrixHelper.getNuc4_4());

        return pair.getPercentageOfIdentity(true);
    }
}
