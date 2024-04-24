'''
This file contains the class Sequence which contains the following
data:
    -Sequence: a string that contains a DNA strand
    -alphabet: an array of the language of the alphabet for this
        strand of DNA
    -sizeOfSequence, sizeOfAlphaber: length of the above components
```` Try a list of dictionaries for the sequences instead of a class ````
'''
class sequence(object):
    def __init__(self):
        self.sequence = ""
        self.alphabet = list()
        self.sizeOfAlphabet = 0
        self.sizeOfSequence = 0
    def __getSizeOfSequence(self):
        return len(self.sequence)
    def __getSizeOfAlphabet(self, size):
        self.sizeOfAlphabet = size
    def setSequence(self, seq):
        self.sequence = seq.strip()
        self.sizeOfSequence = self.__getSizeOfSequence()
    def setAlphabet(self, alph, size):
        self.alphabet = alph
        self.__getSizeOfAlphabet(size)
