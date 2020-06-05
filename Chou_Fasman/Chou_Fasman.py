"""
Chou Fasman Algorithm Implementation 

by Gurpinder Singh
Date: 26 May 2020 09:44 PM 

Note: The implemenation has not been tested for accuracy yet.
It was just created to understand the alogrithm better and for a bit of fun.

Citations:
Chou PY, Fasman GD (1974)
'Prediction of protein conformation'
Biochemistry. 13 (2): 222–245

Table values from: http://prowl.rockefeller.edu/aainfo/chou.htm


Colouring library: Sty 1.0.0b12
Check this awesome library at https://pypi.org/project/sty/

"""
import numpy as np
from sty import fg, rs

    
def finder(Seq, Table, Helix= True):
    """
    Finds alpha helices and beta sheets from the given Sequence using Table values.
    For Helix = True, the window size is set to 6
    For Helix = False, the window size is set to 5
    
    Output: Returns a list of tuples containing start + stop indices.
    """
    if Helix: #parameters for finding Helices
        windowSize = 6
        minReq = 4
        ind = 0
    else: #parameters for finding sheets
        windowSize = 5
        minReq = 3
        ind = 1
        
    results = []
    for i in range(len(Seq)-windowSize):
        window = Seq[i:i+6]
        count = 0
        for aa in window:
            if Table[aa][ind] > 100:
                count += 1
        if count > minReq: #min req for respective helix/sheet
            results.append((i, i+6))
             
    return results


def extender(Seq, Table, start, end, Helix = True):
    """
    Extends the helix/sheet region on both sides
    Inputs:
        Seq: Protein sequnce
        Table: parameter Table
        start = start point
        end = end point
        Helix: True for helix parameters/ False for sheet parameters
        
    Returns:
        Tuple: New corresponding start and end indices
    
    """
    if Helix:
        ind = 0
    else:
        ind = 1
    #while the sum of four residues, starting from end-3 to end+1 is > 100 
    while (sum([Table[aa][ind] for aa in Seq[end-3 : end+1]])/4 > 100) and (end < len(Seq) - 1):
            end += 1
    #while the sum of four residues, starting from start to start +4 is > 100
    while (sum([Table[aa][ind] for aa in Seq[start : start+4]])/4 > 100) and (start > 0):
            start -= 1
    return (start, end)


def turn_finder(Seq, Table):
    """
    Find the residues that might correspond to turns in a sequence.
    Conditions are as per the algorithm.
    Returns a list of indices with predicted turns.
    """
    results = []
    for i in range(len(Seq)-4):
        window = Seq[i:i+4]
        
        #Three Conditions
        Pt = Table[window[0]][3] * Table[window[1]][4] * Table[window[2]][5] * Table[window[3]][6] > 0.000075
        Av_PTurn = sum([Table[aa][2] for aa in window])/4 > 1.0
        PTurn = sum([Table[aa][2] for aa in window]) > max([sum([Table[aa][0] for aa in window]), sum([Table[aa][1] for aa in window])])
        
        #if all thee conditions are true
        if Pt and Av_PTurn and PTurn:
            results.append(i)
    return results

    

def isoverlap(tup1, tup2):
    """
    Determine if two intervals overlap.
    """
    return (tup2[0] <= tup1[0] <= tup2[1]) or (tup1[0] <= tup2[0] <= tup1[1])


def intersect(tup1, tup2):
    """
    Return the common interval between tup1 and tup2. 
    """
    return (max(tup1[0], tup2[0]), min(tup1[1], tup2[0]))



def overlap_remover(tup1, tup2):
    """
    To find the interval in tup1 that does not lie in tup2.
    """
    #if tup1 starts before tup2 and ends within tup2
    if tup1[0] <= tup2[0] and tup1[1] <= tup2[1]:
        #the end part of it needs to be chopped off.
        return [(tup1[0], tup2[0]-1)]
    #if tup1 starts within tup2 and ends after tup2 has ended
    elif tup2[0] <= tup1[0] <=tup2[1] and tup1[1] >= tup2[1]:
        # the starting part needs to be chopped off
        return [(tup2[1]+1, tup1[1])]
    #if tup2 lies entirely within tup1, then we need to output two tuples
    elif tup1[0] <= tup2[0] <= tup1[1] and tup1[0] <= tup2[1] <= tup1[1]:
        return [(tup1[0], tup2[0]-1), (tup2[1]+1, tup1[1])]
    else:
        return ()


def overlap_resolver(Seq, Table, alphas, betas):
    """
    Resolves conflicts between alpha helix and beta sheet sites.
    Input:
        Seq= protein sequene
        Table: parameter table
        alphas: List of predicted alpha indices
        betas: list of predicted beta indices
        
    Returns:
        new lists of alpha helix and beta sheet indices.
    """
    alphas_copy = alphas[:] #create copies of lists for interation
    betas_copy = betas[:]
    
    for alpha in alphas_copy:
        for beta in betas_copy:
            
            if isoverlap(alpha, beta): #if the indices overlap
                inter = intersect(alpha, beta) #check their overlapping region
                Palpha = sum([Table[aa][0] for aa in Seq[inter[0]: inter[1]+1]])
                Pbetas = sum([Table[aa][1] for aa in Seq[inter[0] : inter[1]+1]])
                if Palpha > Pbetas:
                    #it means the region is alpha helix        
                    b = overlap_remover(beta, alpha) #remove the overalpping region from beta
                    for i in b: #overlap remover outputs index tuples inside of lists
                        if i[1] - i[0] > 4: #if the length is atleast 5
                           betas.append(i) #append it to betas
                           try:
                               betas.remove(beta) #remove the onld entry, if it exists.
                           except:
                               continue
                # same as above but for alpha helix this time.           
                if Pbetas > Palpha:
                    #it means the region is beta sheet
                    
                    a = overlap_remover(alpha, beta)
                    for i in a:
                        if i[1] - i[0] > 4:
                            alphas.append(i)
                            try:
                                alphas.remove(alpha)
                            except:
                                continue
    return alphas, betas
                           

def chou_fasman(Seq, Table):
    """
    Final Implementation of Chou Fasman Algorithm
    Inputs:
        Seq: protein sequence
        Table: parameter table
    Returns:
        Lists with alpha helix indices, beta sheet indices and turn indices.
    """
    
    #Find the alphas and betas
    alphas = finder(Seq, Table)
    betas = finder(Seq, Table, Helix=False)
    
    #extend the alphas and betas
    extended_alphas = [extender(Seq, Table, interval[0], interval[1]) for interval in alphas]
    extended_betas = [extender(Seq, Table, interval[0], interval[1], Helix=False) for interval in betas]
    
    #resolving alphas and betas overlaps
    final_alphas, final_betas = overlap_resolver(Seq, Table, extended_alphas, extended_betas)    
    
    #Find Turns in the Sequence
    Turns = turn_finder(Seq, Table)
    
    return final_alphas, final_betas, Turns



def OutputDisplay(Seq, final_alphas, final_betas, Turns):
    """
    Output Code
        Inputs:
            Seq: Protein Sequence
            final_alphas, final_betas, Turns: lists of indices from chou fasman algorithm
        Returns: None        
    """
    #adding blank spaces corresponding to the length of sequence
    
    Helix = np.array([' ' for i in range(len(Seq))])
    Sheets = np.array([' ' for i in range(len(Seq))])
    Turn = np.array([' ' for i in range(len(Seq))])
    
    #adding appropriate symbols at the specified indices
    for interval in final_alphas:
        Helix[np.arange(interval[0], interval[1]+1)] = "H"
    for interval in final_betas:
        Sheets[np.arange(interval[0], interval[1]+1)] = "B"
    Turn[Turns] = "T"
    
    #creating start points to print the sequences
    starts = [60 * i for i in range(len(Seq)) if 60*i < len(Seq)]
    starts.append(len(Seq))
    
    #actual printing
    for i in range(len(starts)-1):
        seq_string = Seq[starts[i] : starts[i+1]]
        Helix_string = fg.blue + str(''.join(list(Helix)[starts[i] : starts[i+1]])) + rs.fg
        Sheet_string = fg.yellow + str(''.join(list(Sheets)[starts[i] : starts[i+1]])) + rs.fg
        Turn_string = fg.green +str(''.join(list(Turn)[starts[i] : starts[i+1]])) + rs.fg
        print("Seq:  \t", seq_string)
        print("Helix:\t", Helix_string)
        print("Sheet:\t", Sheet_string)
        print("Turns:\t", Turn_string, "\n\n")
    
    #small statistics
    print("\nStatistics:")
    print("Length of Protein: {} aa".format(len(protein)))
    print("Helices: {}\tSheets: {}\tTurns: {}\t".format(list(Helix).count("H"), list(Sheets).count("B"), list(Turn).count("T")))
    print("Helices: {}\tSheets: {}\tTurns: {}\t".format(round(list(Helix).count("H")/len(protein), 2),\
                                                        round(list(Sheets).count("B")/len(protein), 2),\
                                                            round(list(Turn).count("T")/len(protein), 2)))

#Table of form AA: [P(a), P(b), P(t), f(i), f(i+1), f(i+2), f(i+3)]

Table = {'A': [142, 83, 66, 0.06, 0.076, 0.035, 0.058],
         'R': [98, 93, 95, 0.070, 0.106, 0.099, 0.085],
         'N': [101, 54, 146, 0.147, 0.110, 0.179, 0.081],
         'D': [67, 89, 156, 0.161, 0.083, 0.191, 0.091],
         'C': [70, 119, 119, 0.149, 0.050, 0.117, 0.128],
         'E': [151, 37, 74, 0.056, 0.060, 0.077, 0.064],
         'Q': [111, 110, 98, 0.074, 0.098, 0.037, 0.098],
         'G': [57, 75, 156, 0.102, 0.085, 0.190, 0.152],
         'H': [100, 87, 95, 0.140, 0.047, 0.093, 0.054],
         'I': [108, 160, 47, 0.043, 0.034, 0.013, 0.056],
         'L': [121, 130, 59, 0.061, 0.025, 0.036, 0.070],
         'K': [114, 74, 101, 0.055, 0.115, 0.072, 0.095],
         'M': [145, 105, 60, 0.068, 0.082, 0.014, 0.055],
         'F': [113, 138, 60, 0.059, 0.041, 0.065, 0.065],
         'P': [57, 55, 152, 0.102, 0.301, 0.034, 0.068],
         'S': [77, 75, 143, 0.120, 0.139, 0.125, 0.106],
         'T': [83, 119, 96, 0.086, 0.108, 0.065, 0.079],
         'W': [108, 137, 96, 0.077, 0.013, 0.064, 0.167],
         'Y': [69, 147, 114, 0.082, 0.065, 0.114, 0.125],
         'V': [106, 170, 50, 0.062, 0.048, 0.028, 0.053]
         }

        
#Main Code
        
print("""Welcome to Chou Fasman Algorithm Implementation
Created on a fine wednesday afternoon, 27th May 2020.
Author: Gurpinder Singh\n""")

print("-"*30)
print("\nNote: The implemenation has not been tested for accuracy yet.\nIt was just created to understand the alogrithm better and for a bit of fun.\n")
print("Citations:\nChou PY, Fasman GD (1974)\n'Prediction of protein conformation'\nBiochemistry. 13 (2): 222–245\n\nTable values from: http://prowl.rockefeller.edu/aainfo/chou.htm\n")
print("Colouring library: Sty 1.0.0b12\nCheck this awesome library at https://pypi.org/project/sty/\n")
print("-"*30)

file = input("Enter the name of the file containing protein Sequence (in FASTA format): ")
meta = ""
protein = ""
with open(file, 'r') as fh:
    for line in fh:
        if ">" in line:
            meta = line
        else:
           protein = protein + line.replace("\n", '')
           
Alphas, Sheets, Turns = chou_fasman(protein, Table)

print("\n\nChouFasman Plot for @", meta.split(' ')[0][1:], "\n")    
OutputDisplay(protein, Alphas, Sheets, Turns)




