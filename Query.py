"""
Program: IRSpectrum.py
Programmed by: Josh Ellis, Josh Hollingsworth, Aaron Kruger, Alex Matthews, and
    Joseph Sneddon
Description: This program will recieve an IR Spectrograph of an unknown
    molecule and use our algorithm to compare that graph to a stored database of
    known molecules and their IR Spectrographs. This program will then return a
    list of the closest Spectrographs matches as determined by our algorithm.
Query.py: This part of the program recieves the file location of a query IR
    Spectrograph downloaded by main.js. formatQueryData() then formats the query
    data and returns a dictionary, queryDict, of the formated query data.
    compareQueryToDB() then takes that dictionary and compares it against all of
    the IR spectrographs imported from our IR spectrum database (IR.db).
    compareQueryToDB() then sends a string back to main.js of the closest IR
    spectrographs found in the IR.db.
"""
#---------------------------------Imports--------------------------------------
import sys
import os
from IR_Functions import *
import multiprocessing as mp
from shutil import copyfile
#------------------------------------------------------------------------------

#---------------------------------Classes/Functions----------------------------
class FormatQueryData:
    def __new__(self, queryPath, comparisonTypes, filename):
        """
        Creates class object and initializes class variables, then returns a
        dictionary of the formated query data.
        """
        self.queryPath = queryPath
        self.comparisonTypes = comparisonTypes
        self.filename = filename

        return self.formatQueryData(self)

    def timeStamp(self, f):
        return int(f.split('.')[0].split('_')[-1])

    def cleanupQueryData(self, images):
        """Removes all generated query data that is more than 5 min old."""
        currentTime = self.timeStamp(self, self.filename)
        holdTime=2*60*1000
        for each in [file for file in os.listdir("public\\uploads")
                        if file.endswith(".jpg")]:
            try:
                if self.timeStamp(self, each)<currentTime-holdTime:
                    os.remove("public\\uploads\\"+each)
            except:
                pass

        #Deletes the temp file downloaded by main.js
        os.remove(images[0])
        if 'temp' in self.queryPath:
            os.remove(self.queryPath)

    def formatQueryData(self):
        #Open the source image
        images = PullImages(self.queryPath)  #PullImages() from IR_Functions.py
        IR_Data = ReadGraph(images[0])  #ReadGraph() from IR_Functions.py

        copyfile(images[0], "public\\uploads\\" + self.filename)

        #Cleans up temp data from queries.
        self.cleanupQueryData(self, images)

        #Calculate each transformation. ConvertQuery() from IR_Functions.py
        queryDict=ConvertQuery(IR_Data,self.comparisonTypes)

        return queryDict
#------------------------------------------------------------------------------

#----------------------------Multiprocessing functions-------------------------
def work(DataQ, ReturnQ, query, comparisonTypes):
    try:
        casNum, dataDict = DataQ.get()

        differenceDict = {}
        for cType in comparisonTypes:
            differenceDict[cType] = []

            dif=Compare(cType,dataDict[cType],query[cType])

            differenceDict[cType] += [(dif, casNum)]
        ReturnQ.put(differenceDict)
        return True
    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        if int(exc_tb.tb_lineno)==83:
            #error due to active update
            ReturnQ.put(None)
            return False
        print('\nERROR!:')
        print('%s' % e)
        print("\n"+str(exc_tb.tb_lineno)+" "+str(exc_obj)+" "+str(exc_tb),"\n")
        return False

def worker(workerNo, JobsDoneQ, NofJobs, NofWorkers, ReturnQ, DataQ, query,
            comparisonTypes):
    # Worker loop
    working = True
    while working:
        jobNo = JobsDoneQ.get()
        work(DataQ, ReturnQ, query, comparisonTypes)
        if NofJobs-jobNo <= NofWorkers-1:
            working = False

def multiProcessController(formatedQueryData,comparisonTypes,IR_Info,dataDict,differenceDict):
    CORES = min(mp.cpu_count(),len(IR_Info))
    
    JobsDoneQ=mp.Queue()
    ReturnQ=mp.Queue()
    ReadRequestQ=mp.Queue()
    DataQ=mp.Queue()
    DataBuffer=min(CORES*2,len(IR_Info))

    for iCompound in range(len(IR_Info)):
        JobsDoneQ.put(iCompound+1)
        ReadRequestQ.put(1)
    for iCompound in range(DataBuffer):
        DataQ.put((IR_Info[iCompound][0], dataDict[IR_Info[iCompound][0]]))
        ReadRequestQ.get()
        ReadRequestQ.put(0)

    p = {}
    for core in range(CORES):
        p[core] = mp.Process(target=worker,
                          args=[core, JobsDoneQ, len(IR_Info), CORES, ReturnQ, DataQ,
                                formatedQueryData, comparisonTypes])
        p[core].start()

    #Read returned data from workers, add new read reqests
    for iCompound in range(DataBuffer, len(IR_Info)+DataBuffer):
        retDict = ReturnQ.get()
        if retDict:
            for cType in comparisonTypes:
                differenceDict[cType] += retDict[cType]
        else:#not found due to active update
            for cType in comparisonTypes:
                differenceDict[cType] += [(0,)]
        if ReadRequestQ.get():
            DataQ.put((IR_Info[iCompound][0], dataDict[IR_Info[iCompound][0]]))

    for core in range(CORES):
        p[core].join()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def importDB():
    try:
        myIRDB = IRDB()
        IR_Info = myIRDB.searchIRDB("SELECT CAS_Num FROM IR_Info GROUP BY CAS_Num")
        IR_Data=myIRDB.searchIRDB("SELECT CAS_Num,Type,Wavelength,Value FROM IR_Data")

        return IR_Info, IR_Data
    except:
        return None

def generateDataDict(IR_Info, IR_Data, comparisonTypes):
    dataDict = {}
    for iCompound in range(len(IR_Info)):
        dataDict[IR_Info[iCompound][0]] = {}
        for cType in comparisonTypes:
            dataDict[IR_Info[iCompound][0]][cType] = []
    for iDBrow in range(len(IR_Data)):
        if 'raw'!=IR_Data[iDBrow][1]:
            dataDict[IR_Data[iDBrow][0]][IR_Data[iDBrow][1]]+=[IR_Data[iDBrow][2:]]
        else:
            for cType in comparisonTypes:
                if 'raw' in cType:
                    dataDict[IR_Data[iDBrow][0]][cType]+=[IR_Data[iDBrow][2:]]
    return dataDict

def generateDifDict(comparisonTypes):
    differenceDict = {}
    for cType in comparisonTypes:
        differenceDict[cType]=[]
    return differenceDict

def compareQueryToDB(formatedQueryData,comparisonTypes):
    IR_Info, IR_Data = importDB()

    dataDict = generateDataDict(IR_Info, IR_Data, comparisonTypes)

    differenceDict = generateDifDict(comparisonTypes)

    multiProcessController(formatedQueryData,comparisonTypes,IR_Info,dataDict,differenceDict)

    #Sort compounds by difference. SmartSortResults() from IR_Functions.py
    results=SmartSortResults(differenceDict,[a[0] for a in IR_Info])[:min(20,len(IR_Info))]
    retString=""

    #Save list of compound differences to file
    for iResult in range(len(results)):
        retString+=results[iResult][1]+" "

    #Gives sorted list of Output to main.js
    return retString.strip()
#------------------------------------------------------------------------------

#---------------------------------Program Main---------------------------------
def main(queryPath, filename):

    if importDB():
        #get comparison types from file
        comparisonTypes=ReadComparisonKeys()
        
        formatedQueryData = FormatQueryData(queryPath,comparisonTypes,filename)

        results = compareQueryToDB(formatedQueryData,comparisonTypes)
        print(results)

        sys.stdout.flush()
    else:
        print("DB_Not_Found")

        sys.stdout.flush()
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
#---------------------------------End of Program-------------------------------
