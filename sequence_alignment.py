'''
This file contains the code to build the sequence matrices and evaluates
them according to the alignment all functions in this file will be directly
run from the main.py file.
'''

import constants
from SequenceClass import sequence as seq

'''
This function reads the file and stores the data in a sequence array
then returns that value from the function. This function should be
called first also so that the file is open.
    Args:
        -file: fileName that may be specified by the user if not it
            will default to the testData.txt file.
'''
def grabSequences(file=constants.TEST):
    seq1 = seq()
    seq2 = seq()
    f = open(file, 'r')
    temp = str(f.readline()) # Grabs first line of file
    seq1.setSequence(temp)
    temp = str(f.readline()) # Grabs second line of file
    seq2.setSequence(temp)
    
    seqArray = [ seq1, seq2 ]
    
    return f, seqArray

'''
This function reads and interprets the data from file and uses the value as
an index in the alignmentChoice list, then returns what value is at the
specified location.
    Args:
        -file: file name that may be specifed by the user if not it will
            default to testData.txt
    alignmentChoice = [ 'local', 'global' ]
'''
def determineAlignment(file):
    if not file.closed:
        index = int(file.readline()) # Grabs third line of file
        alignment = constants.alignmentChoice[index]
        return alignment
    else:
        print("file closed")

'''

This function reads the file and stores the gapPenalty as a floating point
decimal, then returns that value.
    Args:
        -file: file name that may be declared by the user if not it will
            default to testData.txt
'''
def determineGapPenalty(file):
    if not file.closed:
        penalty = float(file.readline()) # Grabs the 4th line of file
        return penalty
'''
This function sets the alphabet of each respective sequence.
    Args:
        -seqArray: array returned from the grabSequence function
        -file: file name that may be determined by the user if not it will
            default to testData.txt
'''
def setAlphabet(seqArray, file):
    if not file.closed:
        for i in range(2):
            size = int(file.readline())# grabs the alphabet size
            temp = str(file.readline())# grabs the alphabet, will be changed to a list
            temp = temp.strip()
            temp = [x for x in temp] # converts to a list
            seqArray[i].setAlphabet(temp, size)

'''
This function is designed to parse and slice the temp variable and return
specific values based on the switch argument. The switch argument takes 3
possible values:
    'i': integer
    'd': double or float
    'c': for character
These values are used to determine what is outputted from the function.
    Args:
        -temp: a single line from a file which is updated and evaluated
            in this function given an input switch
        -switch: condition used to evaluate the state of the ouput from
            the function
'''
def parseTemp(temp, switch):
    value = ""
    if (switch == 'i'):
        i = 0
        while (temp[i] != ' '):
            i+=1
        value = int(temp[:i])
        temp = temp[i+1:]
    elif switch == 'd':
        value = float(temp)
        temp = ""
    else:
        i = 0
        v = ""

        value = v
        temp = temp[i + 4:]
    return value, temp
'''
This function is designed to build a matrix that contains the scoreing
matrix provided in the given file. The matrix should already have its
boarders created prior to calling this function. Be sure to call the
assignBoarders function on the languages of the given sequences prior
to calling this function.
    Args:
        -2dArray: should be initialized in the assignBoarders function
        -file: file name that may be specified by the user otherwise it
            will default to testData.txt
'''
def buildMatchMatrix(scoreArray, file):
    if not file.closed:
        temp = " "
        while temp:
            temp = file.readline()
            if (not temp):
                break
            temp = temp.strip()

            row, temp = parseTemp(temp, 'i')
            col, temp = parseTemp(temp, 'i')
            throwAway, temp = parseTemp(temp, 'c')

            hold, throwAway = parseTemp(temp, 'd')
            scoreArray[row].append(hold)
    return scoreArray

'''
Creates the boarders for the scoring matrix. This matrix will be indexed
to create the actual sequence matrix.
    Args:
        -seq1, seq2: two instances of the sequence class.
'''

def assignBoarders(seq1, seq2):
    scoreArray = [[0]]
    for letter in seq1.alphabet:
        scoreArray[0].append(letter)
    for letter in seq2.alphabet:
        scoreArray.append([letter])
    return scoreArray

'''
Generates a 2D Matrix given two unique sequences. This matrix will be the
main matrix used throughout the rest of the code.
    Args:
        -seq1, seq2: two instances of the sequence class.
        -align: the current alignment
        -gapPenalty: cost for gap extension or matching to a gap
'''
def createBoarderSeq(seq1, seq2, align, gapPenalty):
    seqArray = [[(0, False, False, False)]]
    if align == "local":
        val = 0
        for letter in seq1.sequence:
            seqArray[0].append((val, False, True, False))
        for letter in seq2.sequence:
            seqArray.append([(val, False, False, True)])
    else:
        for i in range(len(seq1.sequence)):
            val = (i + 1) * gapPenalty
            seqArray[0].append((val, False, True, False))
        for i in range(len(seq2.sequence)):
            val = (i + 1) * gapPenalty
            seqArray.append([(val, False, False, True)])
    return seqArray

'''
Indexes the score matrix given two letters. Performs a linear search of
for each letter from a sequence and indexes the score penalty/profit of
each match.
    Args:
        -letter1, letter2: letters from each sequence.
        -scoreArray: scoreing matrix created from the buildMatchMatrix
            function
'''
def indexScore(letter1, letter2, scoreArray):
    i = j = 0
    while letter1 != scoreArray[0][i]:
        i += 1
    while letter2 != scoreArray[j][0]:
        j += 1
    return scoreArray[i][j]

'''
This function fills the rest of the matrix after initializeMatrix and
createBoarderSequence functions are called. This function is designed
to create a local alignment version of the resulting matrix. This
function represents the Smith-Waterman algorithm for sequence alignment.
This function will also return three values the new seqArray, maxScore,
and the location of maxScore.
    Args:
        -seqArray: the generated matrix after being initialized in local
            alignment.
        -scoreArray: the scoreing matrix used for these two alignments.
        -gapPenalty: the penalty for matching to gaps, mismatching
            letters, and extending gaps.
'''
def fillSequenceMatrixLocal(seqMatrix, seqArray, scoreArray, gapPenalty):
    maxScore = 0
    location = [(1, 1)]

    for row in range(1, len(seqMatrix)):
        y = seqArray[1].sequence[row - 1]
        for col in range(1, len(seqMatrix[0])):
            diag = left = up = False
            x = seqArray[0].sequence[col - 1]
        # Diagonal
        m1 = seqMatrix[row - 1][col - 1][0] + indexScore(x, y, scoreArray)
        m1 = round(m1, 1)
        # Up
        m2 = seqMatrix[row - 1][col][0] + gapPenalty
        m2 = round(m2, 1)
        # Left
        m3 = seqMatrix[row][col - 1][0] + gapPenalty
        m3 = round(m3, 1)
        val = max(m1, m2, m3, 0)
        if val == m1:
            diag = True
        if val == m2:
            up = True
        if val == m3:
            left = True
        
        seqMatrix[row].append((val, diag, left, up))

        if maxScore <= seqMatrix[row][col][0]:
            maxScore = seqMatrix[row][col][0]
            location = (row, col)
    return seqMatrix, maxScore, location

'''
This function fills the rest of the matrix provided an initialized matrix in
global alignment. Thus this should be called after the initializeMatrix and
createBoarderSequence functions are called. This function represents the
Needleman-Wunsch algorithm. This function returns 3 values: the resulting
matrix, maxScore, and the location of the maxScore.
    Args:
        -seqArray: an array initialized in global alignment.
        -scoreArray: an array containing the match and mismatch scores
        -gapPenalty: penalty for a gap extension or matching to a gap
'''
def fillSequenceMatrixGlobal(seqMatrix, seqArray, scoreArray, gapPenalty):
    location = (len(seqMatrix) - 1, len(seqMatrix[0]) - 1)
    
    for row in range(1, len(seqMatrix)):
        y = seqArray[1].sequence[row - 1]
        for col in range(1, len(seqMatrix[0])):
            diag = left = up = False
            x = seqArray[0].sequence[col - 1]
            # Diagonal
            m1 = seqMatrix[row - 1][col - 1][0] + indexScore(x, y, scoreArray)
            #m1 = round(m1, 1)
            # Up
            m2 = seqMatrix[row - 1][col][0] + gapPenalty
            #m2 = round(m2, 1)
            # Left
            m3 = seqMatrix[row][col - 1][0] + gapPenalty
            #m3 = round(m3, 1)
            val = max(m1, m2, m3)
            
            if val == m1:
                diag = True
            if val == m2:
                up = True
            if val == m3:
                left = True
            
            seqMatrix[row].append((val, diag, left, up))

    row, col = location
    maxScore = seqMatrix[row][col][0]

    return seqMatrix, maxScore, location
'''
This traceback follows the sequences from the local alginment matrix. It
stores the strings formed by the matrix in a list as a tuple in the form
of (sequence1, sequence2).
    Args:
        -seqMatrix: the local alignment matrix formed from the
            fillSequenceMatrixLocal function
        -seqArray: the array containing the sequences used in this file
        -location: the location of maximum or current score
        -seqMatch: array that will be used to store each alignment
        -seq: a tuple containing the current states of the alignments
'''
def traceBack(seqMatrix, seqArray, location, seqMatch, seq=('','')):
    i, j = location
    # val = (score, diag, left, up)
    # diag, left, up are bool type
    val = seqMatrix[i][j]
    seq1, seq2 = seq
    if val[0] == 0 or (i == 0 and j == 0):
        seqMatch.append(seq)

    if val[1] and val[0] != 0:
        # Handles diagonal
        seq1a = seqArray[0].sequence[j - 1] + seq1
        seq2a = seqArray[1].sequence[i - 1] + seq2
        updated = (seq1a, seq2a)
        traceBack(seqMatrix, seqArray,(i - 1, j - 1), seqMatch, updated)
    if val[2] and val[0] != 0:
        # Handles left
        seq1b = seqArray[0].sequence[j - 1] + seq1
        seq2b = '_' + seq2
        updated = (seq1b, seq2b)
        traceBack(seqMatrix, seqArray,(i, j - 1), seqMatch, updated)
    if val[3] and val[0] != 0:
        # Handles upward node
        seq1c = '_' + seq1
        seq2c = seqArray[1].sequence[i - 1] + seq2
        updated = (seq1c, seq2c)
        traceBack(seqMatrix, seqArray,(i - 1, j), seqMatch, updated)


'''
May not need this function....

'''
def traceBackLocal(seqMatrix, seqList, location, sequences, seq=("", "")):
    pass

'''
This writes the sequences found from this program to a text file along with
the maximum score from the sequence matrix.
    Args:
        -sequences: a list containing tuples of all the sequences with a
            maxScore alignment.
        -maxScore: the maximum score from the sequences.
'''
def createOutput(sequences, maxScore):
    output = open("Output.txt", 'w')
    output.write(str(maxScore) + '\n')
    for seq1, seq2 in sequences:
        output.write(str(seq1) + '\n')
        output.write(str(seq2) + '\n')
        output.write('\n')

'''
This function is strictly for test purposes. It was so that I could make check
that the arrays are filling properly.
    -Args:
        -seqMatrix: the dynamic matrix created from two sequences. Each node
            is in the form of (score, diag, left, up). This allows for an
            easier traceback
'''
def displayMatrix(seqMatrix):
    size = len(seqMatrix[0])
    top = [ x for x in range(size)]
    top = [' '] + top
    
    print(top)
    for row in range(len(seqMatrix)):
        r = "["
        print(row, end=' ')
        for col in range(len(seqMatrix[0])):
            if (col + 1) != len(seqMatrix[0]):
                r += (' ' + str(seqMatrix[row][col][0]) + ',')
            else:
                r += (' ' + str(seqMatrix[row][col][0]) + ' ]')
        print(r)
