'''
This program is designed to perform two different sequence alignments: local and
global. The user can specify the fileName otherwise it defaults to testData.txt.
The sequences can contain up to 300 letters, and each protein must be an
alphabetical literal, as a capital letter. The goal of the program is to score
and display in the results.txt file all of the highest scoring values.

Zachare Lofton
Sequence Alignment Program
Dr. Chung Ng
CSC 440: Design and Analysis of Algorithms
'''
import Sequence_Alignment as seqFile
import SequenceClass as sequence
def promptFileName():
    fileName = str(input("Enter a file name without '.txt'" +
                        "\n(default: testData.txt" +
                        " *Click enter for default*) -> "))
    return fileName + ".txt"

def __main__():
    pass

if __name__ == "__main__":
    fileName = promptFileName()
    sequences = []

    if fileName != ".txt":
        file, seqList = seqFile.grabSequences(fileName)
    else:
        file, seqList = seqFile.grabSequences()

    align = seqFile.determineAlignment(file)
    gap = seqFile.determineGapPenalty(file)
    seqFile.setAlphabet(seqList, file)
    matchMatrix = seqFile.assignBoarders(seqList[0], seqList[1])
    matchMatrix = seqFile.buildMatchMatrix(matchMatrix, file)
    seqMatrix = seqFile.createBoarderSeq(seqList[0], seqList[1], align, gap)
    if align == 'local':
        seqMatrix, maxScore, location = seqFile.fillSequenceMatrixLocal(
                                        seqMatrix, seqList, matchMatrix, gap)
        #seqFile.traceBackLocal(seqMatrix, seqList, location, sequences)
    else:
        seqMatrix, maxScore, location = seqFile.fillSequenceMatrixGlobal(
                                        seqMatrix, seqList, matchMatrix, gap)
    seqFile.traceBack(seqMatrix, seqList, location, sequences)
    seqFile.createOutput(sequences, maxScore)
