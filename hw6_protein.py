"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this
import numpy
import matplotlib


### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file=open(filename,"r")
    res=file.read()
    res=res.replace("\n","")
    return res


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    dna=dna.replace("T","U")
    lst=[]
    for i in range(startIndex,len(dna),3):
        str1=dna[i:i+3]
        if(str1=="UAA" or str1=="UAG" or str1=="UGA"):
            lst.append(str1)
            break
        else:
            lst.append(str1)
    return lst


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    file=open(filename,"r")
    obj=json.load(file)
    dict={}
    for each in obj:
        for word in obj[each]:
            word=word.replace("T","U")
            dict[word]=each
    return dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    lst=[]
    for each in codons:
        if each=="AUG" and lst==[]:
            lst.append("Start")
        elif(each=="UAA" or each=="UAG" or each=="UGA"):
            lst.append("Stop")
            break
        else:
            lst.append(codonD[each])
    return lst


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna=readFile(dnaFilename)
    codondict=makeCodonDictionary(codonFilename)
    counter=0
    proteinlst=[]
    i=0
    while(i!=len(dna)):
        str2=dna[i:i+3]
        if(str2=="ATG"):
            RNA=dnaToRna(dna,i)
            lst=generateProtein(RNA,codondict)
            proteinlst.append(lst)
            i=i+3*len(RNA)
        else:
            counter=counter+1
            i=i+1
    #print("Total Bases:",len(dna)/3)
    #print("Unused Bases:",counter)
    #print(len(proteinlst))
    return proteinlst


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
    resultlst=[]
    for each in proteinList1:
        if each in proteinList2 and each not in resultlst:
            resultlst.append(each)
    return resultlst


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    result=[]
    for protein in proteinList:
        for amino in protein:
            result.append(amino)
    return result


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dict={}
    for each in aaList:
        if each not in dict:
            dict[each]=0
        dict[each]+=1
    return dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    amino1=combineProteins(proteinList1)
    aminodict1=aminoAcidDictionary(amino1)
    amino2=combineProteins(proteinList2)
    aminodict2=aminoAcidDictionary(amino2)
    len1=len(amino1)
    len2=len(amino2)
    dict1={}
    dict2={}
    res=[]
    for each in aminodict1:
        dict1[each]=aminodict1[each]/len1
    for each in aminodict2:
        dict2[each]=aminodict2[each]/len2
    for each in dict1:
        if each not in dict2:
            if (dict1[each]>cutoff):
                res.append([each,dict1[each],0])
        elif(each!="Start" and  each!="Stop"): 
            sub=dict1[each]-dict2[each]
            if(abs(sub)>cutoff):
                res.append([each,dict1[each],dict2[each]])
    for each in dict2:
        if each not in dict1 and dict2[each]>cutoff:
            if(each!="Start" and each!="Stop"):
                res.append([each,0,dict2[each]])
    return res



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    lst1=commonalities[:]
    lst2=differences[:]
    lst3=[]
    for each in lst1:
        str=""
        for word in each:
            if(word!="Start" and word!="Stop"):
                if(len(each)>4 and each.index(word)!=len(each)-2):
                    str=str+word+"-"
                elif(len(each)>4):
                    str=str+word
                else:
                    str=word
        lst3.append(str)
    lst4=lst3.sort()
    for i in lst3:
        print(i)
    print("The following amino acids occurred at very different rates in the two DNA sequence")
    for each in lst2:
        r1=round(each[1]*100,2)
        r2=round(each[2]*100,2)
        print(each[0],":",r1, "% in  seq1,",r2, "% in seq2") 

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
    Aminolst1=combineProteins(proteinList1)
    Aminolst2=combineProteins(proteinList2)
    Aminodict1=aminoAcidDictionary(Aminolst1)
    Aminodict2=aminoAcidDictionary(Aminolst2)
    #print(Aminodict1)
    #print(Aminodict2)
    lst1=[]
    key1=list(Aminodict1.keys())
    key2=list(Aminodict2.keys())
    #print(key1)
    #print(key2)
    lst1=key1[:]
    #print("sample:",lst1)
    for each in key2:
        if each not in key1:
            lst1.append(each)
    #print("sample1:",lst1)
    lst1.sort()
    #print("sample3:",lst1)
    return lst1


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    Aminolst=combineProteins(proteinList)
    Aminodict=aminoAcidDictionary(Aminolst)
    len1=len(Aminolst)
    data=[]
    dict={}
    for each in Aminodict:
        dict[each]=Aminodict[each]/len1
    for each in labels:
        if each in dict:
            data.append(dict[each])
        else:
            data.append(0)
    return data


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


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

