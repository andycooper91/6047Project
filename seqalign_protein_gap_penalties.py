#!/usr/bin/env python

import sys
import math

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
PTR_NONE, PTR_GAP1_1, PTR_GAP2_1, PTR_BASE_1, PTR_GAP1_2, PTR_GAP2_2, PTR_BASE_2,PTR_GAP1_3, PTR_GAP2_3, PTR_BASE_3   = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

def seqalignDPV2(seq1,seq2):
    break_frame_pen = 15
    continue_frame_pen = 4
    fix_frame_pen = -10

    """return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
       Note: gap_pen should be positive (it is subtracted)
    """
    F1 = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    F2 = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    F3 = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    TB1 = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB2 = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB3 = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]




    # initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
    for i in range(1,len(seq1)+1):
        F1[i][0] = -5000
        F2[i][0] = -5000
        F3[i][0] = -5000
        TB1[i][0] = PTR_GAP2_1 # indicates a gap in seq2
        TB2[i][0] = PTR_GAP2_2 # indicates a gap in seq2
        TB3[i][0] = PTR_GAP2_3 # indicates a gap in seq2

    for j in range(1,len(seq2)+1):
        F1[0][j] = -5000
        F2[0][j] = -5000
        F3[0][j] = -5000
        TB1[0][j] = PTR_GAP1_1 # indicates a gap in seq1
        TB2[0][j] = PTR_GAP1_2 # indicates a gap in seq1
        TB3[0][j] = PTR_GAP1_3 # indicates a gap in seq1



    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            match1 = F1[i-1][j-1] + S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
            match2 = F2[i-1][j-1] + S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]
            match3 = F3[i-1][j-1] + S[base_idx[seq1[i-1]]][base_idx[seq2[j-1]]]

            insert1 = F3[i][j-1] - fix_frame_pen
            insert2 = F1[i][j-1] - break_frame_pen
            insert3 = F2[i][j-1] - continue_frame_pen

            delete1 = F2[i-1][j] - fix_frame_pen
            delete2 = F3[i-1][j] - continue_frame_pen
            delete3 = F1[i-1][j] - break_frame_pen

            
            F1[i][j] = max(match1, insert1, delete1)
            F2[i][j] = max(match2, insert2, delete2)
            F3[i][j] = max(match3, insert3, delete3)

            
            if F1[i][j] == match1:
                TB1[i][j] = PTR_BASE_1
            elif F1[i][j] == delete1:
                TB1[i][j] = PTR_GAP2_2
            else:
                TB1[i][j] = PTR_GAP1_3

            if F2[i][j] == match2:
                TB2[i][j] = PTR_BASE_2
            elif F2[i][j] == delete2:
                TB2[i][j] = PTR_GAP2_3
            else:
                TB2[i][j] = PTR_GAP1_1

            if F3[i][j] == match3:
                TB3[i][j] = PTR_BASE_3
            elif F3[i][j] == delete3:
                TB3[i][j] = PTR_GAP2_1
            else:
                TB3[i][j] = PTR_GAP1_2


##    return F[len(seq1)][len(seq2)], F, TB
    return F1[len(seq1)][len(seq2)],F2[len(seq1)][len(seq2)],F3[len(seq1)][len(seq2)], F1, F2, F3, TB1, TB2, TB3

def traceback(seq1,seq2,TB1,TB2,TB3,initial):
    s1 = ""
    s2 = ""
    current = initial

    
    i = len(seq1)
    j = len(seq2)

    while TB1[i][j] != PTR_NONE and TB2[i][j] != PTR_NONE and TB3[i][j] != PTR_NONE :
        if current[i][j] == PTR_BASE_1:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i=i-1
            j=j-1
        elif current[i][j] == PTR_BASE_2:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i=i-1
            j=j-1
        elif current[i][j] == PTR_BASE_3:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i=i-1
            j=j-1
            
        elif current[i][j] == PTR_GAP1_1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j=j-1
            current = TB1
        elif current[i][j] == PTR_GAP1_2:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j=j-1
            current = TB2
        elif current[i][j] == PTR_GAP1_3:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j=j-1
            current = TB3
            
        elif current[i][j] == PTR_GAP2_1:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i=i-1
            current = TB1
        elif current[i][j] == PTR_GAP2_2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i=i-1
            current = TB2
        elif current[i][j] == PTR_GAP2_3:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i=i-1
            current = TB3
        else:
            print current[i][j]
            assert False

    return s1,s2

def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)

S = [
	# A  G   C   T
	[3, -1, -2, -2], # A
	[-1, 3, -2, -2], # G
	[-2, -2, 3, -1], # C
	[-2, -2, -1, 3]  # T
	]

gap_pen = 4

def main(arg1, arg2):
    # parse commandline
##    if len(sys.argv) < 3:
##        print "you must call program as: python ps1-seqalign.py <FASTA 1> <FASTA 2>"
##        sys.exit(1)
##
##    file1 = sys.argv[1]
##    file2 = sys.argv[2]

    seq1 = readSeq(arg1)
    seq2 = readSeq(arg2)
##    seq1 = "AGGTGAT"
##    seq2 = "AGTAA"
##    seq1 = "AAAAAAAAAAAAAAAA"
##    seq2 = "TTTTTTTTTTTTTTTT"
##    seq1 = arg1
##    seq2 = arg2

    ##score, F, TB = seqalignDP(seq1,seq2,S,gap_pen)

    score1, score2, score3, F1, F2, F3, TB1, TB2, TB3 = seqalignDPV2(seq1,seq2)
    score = max(score1,score2,score3)
    print >> sys.stderr, score

    if score1 >= score2 and score1 >= score3:
        s1, s2 = traceback(seq1,seq2,TB1, TB2, TB3, TB1)

    elif score2 >= score3:
        s1, s2 = traceback(seq1,seq2,TB1, TB2, TB3, TB2)

    else:
        s1, s2 = traceback(seq1,seq2,TB1, TB2, TB3, TB3)

    max_score = (min(len(seq1), len(seq2)) * 3) - ((abs(len(seq1)-len(seq2))) * 3)
    min_score = (min(len(seq1), len(seq2)) * -2) - ((abs(len(seq1)-len(seq2))) * 3)
    max_score += float(abs(min_score))
    score += float(abs(min_score))
    print "Normalized score: " + str(score/max_score)
    print s1
    print s2

if __name__ == "__main__":
##    main("mouse_HoxA13.fa", "mouse_HoxD13.fa")
##    main("human_HoxA13.fa", "mouse_HoxA13.fa")
##    main("human_HoxA13.fa", "human_HoxD13.fa")
    main("human_HoxD13.fa", "mouse_HoxD13.fa")
##    main("AGAGT", "AGGC")
##    main("AGGT", "AGGA")
