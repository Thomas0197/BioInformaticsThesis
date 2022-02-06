# -*- coding: utf-8 -*-
"""
Created on Mon May 24 18:07:42 2021

@author: Thomas
"""
import sys

#Retracing the methylation genes using the UC.x markers to the differential methylation sites they are affected from
# and placing combining the information into a single file

def some(filepath, filepath2, filepath3):
    count = 0
    FullList = []
    
    Genes = open(filepath, "r")
    Input = open(filepath2, "r")
    meth = open(filepath3, "r")
    
    f = open(r"./siteGeneComparisons.tsv" ,"x")
    f.write("chr" + "\t" + "start" + "\t" + "end" + "\t" + "strand" + "\t" + "pvalue" + "\t" + "qvalue" + "\t" + "meth.diff" + "\t" +"ID" +"\t" + "gene" + "\t" + "distance" + "\n")
            
            
            
    for line in Genes:
        linea = line.rsplit()
        #print(line)
        #ID = linea[1]
        Gene = linea[0]
        #Dist = linea[2]
        
        #WholeLine = line
        Input = open(filepath2, "r")
        for line in Input:
             line = line.rstrip()
             LineList = line.split()
             UC = LineList[3]
             #PosStart = LineList[1]
             j=0
             for i in linea:
                 j=j+1
                 if (UC == i):
                     #print(j)
                     count = count + 1
                     LineList.append(Gene)
                     LineList.append(linea[j])
                     FullList.append(LineList)
        
    for line in meth:
        line = line.rstrip()
        List = line.split()
        myDiff = List[1]
        #print(line)
        for i in FullList:
            #print(i[1])
            #print (myDiff)
            if myDiff == i[1]:
                #print("ok")
                print(str(List[0]) +"\t" +str(List[1]) +"\t" + str(List[2])+"\t" +str(List[3]) + "\t" +str(List[4]) + "\t" + str(List[5]) + "\t" + str(List[6]) + "\t" + str(i[3]) +  "\t" + str(i[4]) +  "\t" + str(i[5]) +"\n")
                f.write(str(List[0]) +"\t" +str(List[1]) +"\t" + str(List[2])+"\t" +str(List[3]) + "\t" +str(List[4]) + "\t" + str(List[5]) + "\t" + str(List[6]) + "\t" + str(i[3]) +  "\t" + str(i[4]) + "\t" + str(i[5]) + "\n")
    
    #print(str(List[0]) +"\t" +str(List[1]) +"\t" + str(List[2])+"\t" +str(List[3]) + "\t" +str(List[4]) + "\t" + str(List[5]) + "\t" + str(List[6]) + "\t" + str(i[3]) +  "\t" + str(i[4]))
    #print(FullList)
    
            
    
    
#some(r"C:\Users\Thomas\Downloads\geneComparisons.txt", r"C:\Users\Thomas\Downloads\GREAT_Input_40.tsv", r"C:\Users\Thomas\Downloads\myDiff40_mWT-Df1-Sup-AAD_VS_mDf1-Df1-Sup-AAD.tsv")

if __name__ == "__main__":
  filename=sys.argv[1]
  filename2=sys.argv[2]
  filename3=sys.argv[3]
  some(filename, filename2,filename3)
