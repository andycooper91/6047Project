#!/usr/bin/env python

import sys

PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3

def readCodons(seq1, seq2, offset1, offset2):
    aa1 = ""
    for i in range(offset1, len(seq1)-2, 3):
        aa1 += codons[seq1[i:i+3]]

    aa2 = ""
    for j in range(offset2, len(seq2)-2, 3):
        aa2 += codons[seq2[j:j+3]]

    return aa1, aa2

def seqalignDP(aa1, aa2):
    """return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
       Note: gap_pen should be positive (it is subtracted)
    """
    
    F = [[0 for j in range(len(aa2)+1)] for i in range(len(aa1)+1)]
    TB = [[PTR_NONE for j in range(len(aa2)+1)] for i in range(len(aa1)+1)]

    # initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
    for i in range(1,len(aa1)+1):
        F[i][0] = 0
        TB[i][0] = PTR_NONE
    for j in range(1,len(aa2)+1):
        F[0][j] = 0
        TB[0][j] = PTR_NONE


    # Fill in the dynamic programming tables F and TB, starting at [1][1]
    # Hints: The first row and first column of the table F[i][0] and F[0][j] are dummies
    #        Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and seq2[j-1].
    #        Use the dictionary base_idx to convert from the character to an index to
    #         look up entries of the substitution matrix.
    
    for i in range(1, len(aa1)+1):
        for j in range(1, len(aa2)+1):
            match = F[i-1][j-1] + AA[aa_idx[aa1[i-1]]][aa_idx[aa2[j-1]]]
            delete = F[i-1][j] - gap_pen
            insert = F[i][j-1] - gap_pen
            
            F[i][j] = max(match, insert, delete, 0)
            if F[i][j] == match:
                TB[i][j] = PTR_BASE
            elif F[i][j] == delete:
                TB[i][j] = PTR_GAP2
            elif F[i][j] == insert:
                TB[i][j] = PTR_GAP1
            else:
                TB[i][j] = PTR_NONE

    return F[len(aa1)][len(aa2)], F, TB

def traceback(seq1, seq2, aa1, aa2, offset1, offset2, F, TB):
    s1 = ""
    s2 = ""

    a1 = ""
    a2 = ""

    i = len(aa1)
    j = len(aa2)

    high = 0
    for x in range(len(aa1)):
        for y in range(len(aa2)):
            if F[x][y] >= high:
                high = F[x][y]
                i = x
                j = y
    print F[i][j]

    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[3*(i-1)+offset1 : 3*i+offset1] + s1
            s2 = seq2[3*(j-1)+offset2 : 3*j+offset2] + s2
            a1 = aa1[i-1:i] + a1
            a2 = aa2[j-1:j] + a2
            i = i-1
            j = j-1
        elif TB[i][j] == PTR_GAP1:
            s1 = '---' + s1
            s2 = seq2[3*(j-1)+offset2 : 3*j+offset2] + s2
            a1 = '-' + a1
            a2 = aa2[j-1:j] + a2
            j = j-1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[3*(i-1)+offset1 : 3*i+offset1] + s1
            s2 = '---' + s2
            a1 = aa1[i-1:i] + a1
            a2 = '-' + a2
            i = i-1
        else: assert False

    return s1, s2, a1, a2

def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)

def maxScore(aa1, aa2):
    # for each elt in the shorter sequence, look up the perfect match
    # subtract gap penalties for difference in length
    if len(aa1) == len(aa2):
        score1 = calcMax(aa1)
        score2 = calcMax(aa2)
        print "MAX: Lengths equal. Scores " + str(score1) + " and " + str(score2)
        return max(score1, score2)
    else:
        diffL = abs(len(aa1) - len(aa2))
        if len(aa1) < len(aa2):
            score = calcMax(aa1) - gap_pen * diffL
            print "MAX: aa1 shorter. Score " + str(score)
        else:
            score = calcMax(aa2) - gap_pen * diffL
            print "MAX: aa2 shorter. Score " + str(score)

        return score

def calcMax(seq):
    matchSum = 0

    for i in range(len(seq)):
        matchSum += AA[aa_idx[seq[i]]][aa_idx[seq[i]]]

    return matchSum

def minScore(aa1, aa2):
    # for each elt in the shorter sequence, look up the worst match
    # subtract gap penalties for difference in length
    if len(aa1) == len(aa2):
        score1 = calcMin(aa1)
        score2 = calcMin(aa2)
        print "MIN: Lengths equal. Scores " + str(score1) + " and " + str(score2)
        return min(score1, score2)
    else:
        diffL = abs(len(aa1) - len(aa2))
        if len(aa1) < len(aa2):
            score = calcMin(aa1) - gap_pen * diffL
            print "MIN: aa1 shorter. Score " + str(score)
        else:
            score = calcMin(aa2) - gap_pen * diffL
            print "MIN: aa2 shorter. Score " + str(score)

        return score

def calcMin(seq):
    matchSum = 0

    for i in range(len(seq)):
        matchSum += min(AA[aa_idx[seq[i]]][:20])

    return matchSum

## Amino Acid Scores BLOSUM
AA = [
    # Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val STP
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -10], #Ala/A
    [ -1,  5,  0, -2, -3,  1,  0, -1,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -10], #Arg/R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -10], #Asn/N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -10], #Asp/D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -10], #Cys/C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -10], #Gln/Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -2, -1, -2, -10], #Glu/E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -10], #Gly/G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -10], #His/H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -2, -1,  1, -10], #Ile/I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -10], #Leu/L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -10], #Lys/K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -10], #Met/M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -10], #Phe/F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -10], #Pro/P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -10], #Ser/S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -10], #Thr/T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -10], #Trp/W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -10], #Tyr/Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -10], #Val/V
    [-10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, 15] #STP/X
    ]

gap_pen = 12

aa_idx = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4,
          'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9,
          'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
          'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'X': 20}

codons = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
          'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
          'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X',
          'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W',
          
          'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
          'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
          'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
          'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
          
          'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
          'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
          'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
          'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
          
          'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
          'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
          'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
          'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

def main(arg1, arg2):
    seq1 = readSeq(arg1)
    seq2 = readSeq(arg2)
    
    for off1 in range(3):
        for off2 in range(3):
            aa1, aa2 = readCodons(seq1, seq2, off1, off2)

            score, F, TB = seqalignDP(aa1, aa2)

            print "Offset 1: " + str(off1) + ", Offset2: " + str(off2)
            print "Score: " + str(score)

            normalized = float(score) / (maxScore(aa1, aa2))
            print "Normalized Score: " + str(normalized)

            s1, s2, a1, a2 = traceback(seq1, seq2, aa1, aa2, off1, off2, F, TB)
            print s1
            print s2
            print "------------"
            print a1
            print a2
            print ""

if __name__ == "__main__":
##    main("mouse_HoxA13.fa", "mouse_HoxD13.fa")
##    main("human_HoxA13.fa", "mouse_HoxA13.fa")
##    main("human_HoxA13.fa", "human_HoxD13.fa")
    main("human_HoxD13down900.fa", "mouse_HoxD13down900.fa")
##    main("human_HoxD13v2.fa", "mouse_HoxD13v2.fa")
##    main("human_HoxD13messy.fa", "mouse_HoxD13messy.fa")
