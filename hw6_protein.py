"""
Protein Sequencing Project
Name:
Roll Number:
"""

from os import read
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
   open_=open(filename,"r").read().splitlines()
   string="".join(open_)
   return string 


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    res=dna[startIndex:]
    list_=[]
    string=""#atgaug
    ignore=["UAG", "UAA","UGA"]
    for letter in range(len(res)):
        if len(string)!=3:
            string+=res[letter]
        if len(string)==3:
            x= string.replace("T","U")
            if x in ignore:
                list_.append(x)
                return list_
            else:
                list_.append(x)
                string=""
    return list_


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f = open(filename)
    read = json.load(f)
    d={}
    for x,y in read.items():
        for i in y:
            d[i.replace('T','U')]=x
    return d


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    list_=[]
    if codons[0] =="AUG":
        list_.append("Start")
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            list_.append(codonD[codons[i]])
    return list_


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    reading_file = readFile(dnaFilename) #returs string
    making_dic = makeCodonDictionary(codonFilename) # returns dict
    i=0
    count=0
    temp=[]
    while i < len(reading_file):
        if reading_file[i:i+3] == "ATG":
            dna_list= dnaToRna(reading_file,i) #returns list
            prot = generateProtein(dna_list, making_dic) #list of strings
            temp.append(prot)
            i = i+3*len(dna_list)
        else:
            i+=1
            count+=1
    return temp


# def runWeek1():
#     print("Human DNA")
#     humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
#     print("Elephant DNA")
#     elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    temp=[]
    for arow in proteinList1:
        for brow in proteinList2:
            if arow==brow and arow not in temp:
                temp.append(arow)

    return temp


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    _list=[]
    for i in proteinList:
        for word in i:
            _list.append(word)
    # print(_list)
    return _list


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    result={}
    for x in aaList:
        if x not in result:
            result[x]=1
        else:
            result[x]+=1
    return result


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    combine1,combine2 = combineProteins(proteinList1), combineProteins(proteinList2) #returns 1d list
    dict1,dict2 = aminoAcidDictionary(combine1),aminoAcidDictionary(combine2) #returns dictinory
    temp,result=[],[]             #holds acids,#final result list
    freq_dict1,freq_dict2={},{}   #proteinList1,#proteindict2   
    for i in dict1:
        freq_dict1[i] = dict1[i]/len(combine1)
        if i not in temp and i !="Start" and i!="Stop":
            temp.append(i)
    for a in dict2:
        freq_dict2[a] = dict2[a]/len(combine2)
        if a not in temp and a !="Start" and a!="Stop":
            temp.append(a)
    # print(freq_dict1,len(freq_dict1),"f1")
    # print(freq_dict2,len(freq_dict2),"f2")
    for ac in temp:
        freq1,freq2=0,0
        if ac in freq_dict1:
            freq1= freq_dict1[ac]
        if ac in freq_dict2:
            freq2= freq_dict2[ac]
        difference = freq2-freq1
        if difference < -cutoff or difference > cutoff  :
            result.append([ac , freq1, freq2])
    return result


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("Printing the Commnalities!")
    for i in sorted(commonalities):
        commonProteins = ""
        let = i[1:len(i)-1] # to remove start and stop [start, ser, pro, leu,stop]
        count=0
        for j in let:
            commonProteins+=j   #ser,pro,leu
            count+=1            #1,2,3
            if count !=len(let): #not equal
                commonProteins+="-" #ser-pro-leu
        print(commonProteins)
        '''
        for i in commonalities:
            if i =="Start" or i =="Stop":
                pass
            else:
                print(i)
        '''

    print("Printing DNA sequences!")
    for item in differences:
        print(item[0],":",round(item[1]*100,2),"% in seq1,",round(item[2]*100,2),"% in seq2")
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
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


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
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()

    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
    #test.testReadFile()
    #test.testDnaToRna()
    #test.testMakeCodonDictionary()
    #test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    test.runWeek2()

