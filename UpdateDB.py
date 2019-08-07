"""
Program: IRSpectrum.py
Programmed by: Josh Ellis, Josh Hollingsworth, Aaron Kruger, Alex Matthews, and
    Joseph Sneddon
Description: This program will recieve an IR Spectrograph of an unknown
    molecule and use our algorithm to compare that graph to a stored database of
    known molecules and their IR Spectrographs. This program will then return a
    list of the closest Spectrographs matches as determined by our algorithm.
UpdateDB.py: This part of the program imports all pdf files from */IR_samples
    and updates the database (IR.db) with each new compound found.
"""
#---------------------------------Imports--------------------------------------
import sys
import sqlite3
import os
from PIL import Image
from shutil import copyfile
import multiprocessing as mp
import time
from IR_Functions import *
#------------------------------------------------------------------------------

#---------------------------------Classes/Functions----------------------------
def initializeDB():
    #If IR.db somehow gets deleted then re-create it.
    if not os.path.exists("IR.db"):
        file = open('IR.db', 'w+')
        file.close()

    sqlData = "CREATE TABLE IF NOT EXISTS `IR_Data` ( `CAS_Num` TEXT, `Type` \
                TEXT, `Wavelength` NUMERIC, `Value` NUMERIC )"
    sqlInfo = "CREATE TABLE IF NOT EXISTS `IR_Info` ( `Spectrum_ID` TEXT, \
                `CAS_Num` TEXT, `Formula` TEXT, `Compound_Name` TEXT, \
                PRIMARY KEY(`Spectrum_ID`) )"

    myIRDB = IRDB()
    myIRDB.writeIRDB(sqlData)
    myIRDB.writeIRDB(sqlInfo)
    myIRDB.commitIRDB()

def tryWork(Jobs,comparisonTypes):
    try:
        file = Jobs.get()

        """ Open the source image """
        images = PullImages(file)

        fname = file.split("\\")[-1]
        casNum = fname.split(".")[0]

        """ is this file already in the database? """
        myIRDB = IRDB()
        sqlQ = "SELECT CAS_Num FROM IR_Info WHERE CAS_Num='"+casNum+"'"
        sqlData = "SELECT CAS_Num FROM IR_Data WHERE CAS_Num='"+casNum+"'"
        sqlInfo = "INSERT INTO IR_Info(Spectrum_ID, CAS_Num, Formula, \
                                        Compound_Name) VALUES (?, ?, ?, ?)"

        myIRDB.writeIRDB(sqlQ)
        myIRDB.writeIRDB(sqlData)
        qData = myIRDB.fetchallIRDB()

        """ if not in the database set the flag to add it """
        if len(qData)==0:

            copyfile(images[0],"public\\images\\"+casNum+".jpg")

            structure=PullStructure(file)[0]
            CleanStructure(structure)
            copyfile(structure,"public\\info\\"+structure.split("\\")[-1])
            os.remove(structure)

            values=PullText(file)
            #Save compound data into the database
            dbvalues = (list(values.values())[0], casNum,
                        list(values.values())[2], list(values.values())[3])

            myIRDB.writeIRDB(sqlInfo, dbvalues)
            myIRDB.commitIRDB()

            f=open("public\\info\\"+casNum+".json",'w')
            f.write(str(values).replace("'",'"'))
            f.close()
        else:
            os.remove(images[0])
            return casNum+" already in DB"

        data = ReadGraph(images[0])  #ReadGraph() from IR_Functions.py
        os.remove(images[0])

        #calculate each transformation
        comparisonDict={}
        for cType in comparisonTypes:
            comparisonDict[cType]=Convert(data,cType)

        sqlQ = "INSERT INTO IR_Data(CAS_Num, Type, Wavelength, Value) \
                    VALUES (?, ?, ?, ?)"
        #save each transformation to file
        for cType in comparisonDict:
            d=[]
            for row in comparisonDict[cType]:
                d+=[str(row[0])+','+str(row[1])]
                dbvalues = (casNum, cType, row[0], row[1])
                myIRDB.writeIRDB(sqlQ, dbvalues)
                #save data

        myIRDB.commitIRDB()
        return casNum+" added to DB"

    except Exception as e:
        print('\nERROR!:')
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print('%s' % e)
        print("\n"+str(exc_tb.tb_lineno)+" "+str(exc_obj)+" "+str(exc_tb),"\n")
        return False
#------------------------------------------------------------------------------

#----------------------------Multiprocessing functions-------------------------
def worker(Jobs,workerNo,NofWorkers,JobsDoneQ,NofJobs,comparisonTypes):
    working=True
    while working:
        message=tryWork(Jobs,comparisonTypes)
        if message:
            jobNo=JobsDoneQ.get()
            print("[Worker No. "+str(workerNo)+"] "+str(jobNo)+" of "
                    +str(NofJobs)+" "+message)
            if NofJobs-jobNo <= NofWorkers-1:
                working = False
        else:
            working=False

def multiProcessUpdater(comparisonTypes):
    filedir=[os.path.join("IR_samples",file) for file in
                os.listdir("IR_samples") if file.endswith(".pdf")]

    Jobs=mp.Queue()
    JobsDoneQ=mp.Queue()
    for i in range(len(filedir)):
        Jobs.put(filedir[i])
        JobsDoneQ.put(i+1)

    CORES = min(mp.cpu_count(),len(filedir))
    p={}
    print("Starting")
    start=time.time()
    for core in range(CORES):
        p[core] = mp.Process(target = worker, args=[Jobs,core,CORES,JobsDoneQ,len(filedir),
                                                    comparisonTypes])
        p[core].start()
    for core in range(CORES):
        p[core].join()
    print("Done and Done "+str(time.time()-start))
#------------------------------------------------------------------------------

#---------------------------------Program Main---------------------------------
def main():

    comparisonTypes=ReadComparisonKeys()

    #Edits comparisonTypes to include only a single raw
    #comparisons with the raw argument will be calculated in the future.
    raws=[]
    for icomp in range(len(comparisonTypes)-1,-1,-1):
        if 'raw' in comparisonTypes[icomp]:
            raws+=[comparisonTypes.pop(icomp)]
    if len(raws)>0:
        comparisonTypes+=['raw']

    initializeDB()

    multiProcessUpdater(comparisonTypes)

if __name__ == "__main__":
    main()
#---------------------------------End of Program-------------------------------
