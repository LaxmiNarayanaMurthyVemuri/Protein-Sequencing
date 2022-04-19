"""
Protein Sequencing Project
Name:
Roll Number:
"""

from tkinter.ttk import LabeledScale
import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file=open(filename,"r")
    content=file.read()
    #print(type(content))
    string_without_line_breaks = ""
    for line in content:
        stripped_line = line.rstrip()
        string_without_line_breaks+= stripped_line
    file.close()
    #print(string_without_line_breaks)
    return string_without_line_breaks
    # #splitlines() will reaturn a list. 
'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    list=[]
    n=3
    Dna=dna.replace("T","U")
    split_Dna=[Dna[i:i+n] for i in range(startIndex,len(Dna),n)]
    for i in range(len(split_Dna)):
        #print(split_Dna[i])
        list.append(split_Dna[i])
        if (split_Dna[i]=="UAA" or split_Dna[i]=="UAG" or split_Dna[i]=="UGA"):
            break
    #print(list)
    return list
'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    dict={}
    import json
    f=open(filename)
    data=json.load(f)
    for key,value in data.items():
        for j in range(len(value)):
            x=value[j].replace("T","U")
            dict[x]=key
    #print(dict)
    return dict
'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    list=["Start"]
    for i in range(1,len(codons)):
        list.append(codonD[codons[i]])
    #print(list)
    return list
'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    proteins_list=[]
    Data=readFile(dnaFilename)
    count=0
    i=0
    while(i<len(Data)):
        if(Data[i:i+3]=="ATG"):
            RNA=dnaToRna(Data,i)
            #print(RNA)
            Dict=makeCodonDictionary(codonFilename)
            #print(Dict)
            proteins=generateProtein(RNA,Dict)
            proteins_list.append(proteins)
            i=i+3*len(RNA)
        else:
            i=i+1
            count+=1
    return proteins_list

def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###
'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    common=[]
    for i in range(len(proteinList1)):
        if proteinList1[i] in proteinList2:  
            if proteinList1[i] not in common:
                common.append(proteinList1[i])
    return common 
'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    list=[]
    for i in range(len(proteinList)):
        for j in range(len(proteinList[i])):
            list.append(proteinList[i][j])
    #print(list)
    return list
'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dict={}
    for i in range(len(aaList)):
        if aaList[i] not in dict:
            dict[aaList[i]]=1
        else:
            dict[aaList[i]]+=1
    return dict
'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    list=[]
    x=aminoAcidDictionary(combineProteins(proteinList1))
    y=aminoAcidDictionary(combineProteins(proteinList2))
    for key,value in x.items():
        x[key]=value/len(combineProteins(proteinList1))
        if key not in y:
            y[key]=0 
    for key,value in y.items():
        y[key]=value/len(combineProteins(proteinList2)) 
        if key not in x:
            x[key]=0    
        if abs(x[key]-y[key])>cutoff:
            #ignore=["Start","Stop"]
            if key!= "Start" and key!= "Stop" :
                list.append([key, x[key], y[key]])
    return list
'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    #print(commonalities)
    lst = sorted(commonalities)
    print("common proteins are:")
    for i in range(len(lst)):
        for j in range(len(lst[i])):
            if lst[i][j] == "Start" or lst[i][j] == "Stop":
                continue
            else:
                print(lst[i][j], end=' ')
        print( )       
    #print(differences)        
    print("amino acid occured at modt differences are:")
    for i in range(len(differences)):
        print(differences[i][0],":",round(differences[i][1]*100,2),"% in seq1,",round(differences[i][2]*100,2),"% in seq2")
        print( )
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    common=[]
    x=combineProteins(proteinList1)
    y=combineProteins(proteinList2)
    x.extend(y)
    #print(x)
    for i in range(len(x)):
        if x[i] not in common:
            common.append(x[i])
    #print(sorted(common))
    return sorted(common) 


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    list=[]
    #print(labels)
    #print(proteinList)
    com=combineProteins(proteinList)
    D=aminoAcidDictionary(com)
    #print(D)
    for key in labels:
        if key in D:
            list.append(D[key]/len(com))
        else:
            list.append(0)
    #print(list)
    return list


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList):
    import matplotlib.pyplot as plt
    #xLabels=makeAminoAcidLabels(label1,label2)
    w = 0.35  # the width of the bars

    plt.bar(xLabels, freqList1, width=-w, align='edge', label=label1, edgecolor=edgeList)
    plt.bar(xLabels, freqList2, width= w, align='edge', label=label2, edgecolor=edgeList)

    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("Graph")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    #print(biggestDiffs)
    list = []
    edge = []
    for i in biggestDiffs:
        list.append(i[0])
    for j in labels:
        if j in list:
            edge.append("black")
        else:
            edge.append("white")
    return edge


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")

    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()

