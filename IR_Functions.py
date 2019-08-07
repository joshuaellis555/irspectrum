"""
Program: IRSpectrum.py
Programmed by: Josh Ellis, Josh Hollingsworth, Aaron Kruger, Alex Matthews, and
    Joseph Sneddon
Description: This program will recieve an IR Spectrograph of an unknown
    molecule and use our algorithm to compare that graph to a stored database of
    known molecules and their IR Spectrographs. This program will then return a
    list of the closest Spectrographs matches as determined by our algorithm.
IR_Functions.py: This part of the program contains most of the functions used by
    Query.py and UpdatedDB.py.
"""
#---------------------------------Imports--------------------------------------
import PyPDF2
import sqlite3
from PIL import Image
import sys
import warnings
import os

warnings.filterwarnings("ignore")
#------------------------------------------------------------------------------

#---------------------------------Variables------------------------------------

#------------------------------------------------------------------------------

#---------------------------------Classes/Functions----------------------------
def PullImages(filename):
    '''
    Pull graph image from first page of PDF
    '''
    file = PyPDF2.PdfFileReader(open(filename, "rb"))
    xObject = file.getPage(0)


    xObject = xObject['/Resources']['/XObject'].getObject()

    images=[]

    for obj in xObject:

        if xObject[obj]['/Subtype'] == '/Image':
            size = (xObject[obj]['/Width'], xObject[obj]['/Height'])
            data = xObject[obj]._data
            if xObject[obj]['/ColorSpace'] == '/DeviceRGB':
                mode = "RGB"
            else:
                mode = "P"

            if xObject[obj]['/Filter'] == '/FlateDecode':
                img = Image.frombytes(mode, size, data)
                img.save(filename + ".png")
                images+=[filename + ".png"]
            elif xObject[obj]['/Filter'] == '/DCTDecode':
                img = open(filename + ".jpg", "wb")
                img.write(data)
                img.close()
                images+=[filename + ".jpg"]
            elif xObject[obj]['/Filter'] == '/JPXDecode':
                img = open(filename + ".jp2", "wb")
                img.write(data)
                img.close()
                images+=[filename + ".jp2"]
    return images

def PullStructure(filename):
    '''
    Pulls the image of the molecular structure from page 2 as a png
    '''
    file = PyPDF2.PdfFileReader(open(filename, "rb"))
    xObject = file.getPage(1)


    xObject = xObject['/Resources']['/XObject'].getObject()

    images=[]

    for obj in xObject:
        if xObject[obj]['/Subtype'] == '/Image':
            size = (xObject[obj]['/Width'], xObject[obj]['/Height'])
            data = xObject[obj].getData()
            if xObject[obj]['/Filter'] == '/FlateDecode':
                img = Image.frombytes("P", size, data)
                img.save(filename.split('.')[0] + ".png")
                images+=[filename.split('.')[0] + ".png"]
    return images

def PullText(filename):
    '''
    Pull text from the first page of a PDF
    returns an array containing:
    [ SpectrumID, CAS Number, Molecular Formula, Compound Name ]
    '''
    specID = ""
    cas = ""
    formula = ""
    name = ""

    try:
        file = PyPDF2.PdfFileReader(open(filename, "rb"))
        page = file.getPage(0)

        page_content = page.extractText()

        idIndex = page_content.find("Spectrum ID")
        casIndex = page_content.find("CAS Registry Number")
        formulaIndex = page_content.find("Formula")
        nameIndex = page_content.find("CAS Index Name")
        sourceIndex = page_content.find("Source")
        startIndex = casIndex

        begin = idIndex + 11
        end = casIndex
        while begin != end:
            specID += page_content[begin]
            begin += 1

        begin = casIndex + 19
        end = formulaIndex
        while begin != end:
            cas += page_content[begin]
            begin += 1

        begin = formulaIndex + 7
        end = nameIndex
        while begin != end:
            formula += page_content[begin]
            begin += 1

        begin = nameIndex + 14
        end = sourceIndex
        while begin != end:
            name += page_content[begin]
            begin += 1
    except:
        print("There was an error extracting text from the PDF")

    #return [specID, cas, formula, name]
    return {"spectrumID":specID, "cas":cas, "formula":formula, "name":name}

def CleanStructure(filename):
    '''
    Changes all of the brightest pixels in a compound structure image
    to full alpha
    '''
    img = Image.open(filename)
    imgdata=list(img.getdata())#the pixels from the image

    img = Image.new('RGBA', (img.size[0],img.size[1]))

    imgdata=[(i,i,i,255)  if i<31 else (i,i,i,0) for i in imgdata]

    img.putdata(imgdata)
    img.save(filename)

def ReadComparisonKeys():
    f=open("public\\types.keys",'r')
    transformTypes=f.readlines()
    
    f.close()
    transformTypes=[line for line in \
                    [lines.strip() for lines in transformTypes] \
                    if len(line)]
    
    return transformTypes
    

class ReadGraph:
    '''
    Reads each datapoint in the graph and converts it to an x,y coordinate
    Each datapoint gets added to a list and returned
    '''
    def __new__(self, image):
        '''area of image scanned for data'''
        self.image = image
        self.xMin=200
        self.xMax=4100
        self.xRange=self.xMax-self.xMin #the x-range of the graph.
        self.yMin=1.02
        self.yMax=-0.05
        self.yRange=self.yMax-self.yMin #the y-range of the graph.
        #This is the width and height standard for all IR samples
        self.width=1024
        self.height=768
        #the area of each image that we want (the graph)
        self.targetRect=(113,978,29,724) #(left,right,top,bottom)

        return self.readGraph(self)

    #copies pixels from the source image within the targetRect
    def cropRect(self, source):
        left,right,top,bottom=self.targetRect
        newImg=[]
        for y in range(top,bottom+1):
            for x in range(left,right+1):
                newImg+=[source[y*self.width+x]]
        return newImg

    #checks if the pixel at x,y is black
    def pix(self, graph,x,y):
        r,g,b=graph[y*self.width+x]
        if r+g+b>=100:
            return False#not black
        else:
            return True#black

    #These two functions convert graph x,y into scientific x,y
    def convertx(self, x):
        return self.xMin+self.xRange*(x/self.width)
    def converty(self, y):
        return self.yMin+self.yRange*(y/self.height)

    def convertGraph(self, graph):
        """
        Creates a graphData list by finding each black pixel on the x axis. For each
        x get the y range over which the graph has black pixels or None if the graph
        is empty at that x value. It stores the min and max y values in the
        graphData list. Then returns the filled graphData List.
        """
        graphData=[]#to be filled with values from graph
        #For each x get the y range over which the graph has black pixels
        # or None if the graph is empty at that x value
        for x in range(0,self.width):
            graphData+=[None]
            foundPix=False#have you found a pixel while looping through the column
            for y in range(0,self.height):
                p=self.pix(self,graph,x,y)#is the pixel black
                if p and not foundPix:
                    #record the first black pixels y value
                    foundPix=True
                    maxVal=y
                elif not p and foundPix:
                    #record the last black pixels y value
                    minVal=y
                    graphData[-1]=(minVal,maxVal)#write these values to data
                    break#next x

        return graphData

    #convert graph into datapoints
    def cleanData(self, graphData):
        data=[]
        for x in range(len(graphData)):
            #Points in format x,y
            if graphData[x]:
                data+=[(self.convertx(self,x),self.converty(self,graphData[x][1]))]

        return data

    def readGraph(self,):
        #Crops the image
        img = Image.open(self.image)
        imgdata=list(img.getdata())#the pixels from the image

        #The graph is cut out of the larger image
        graph=self.cropRect(self,imgdata)

        #width and height of out cropped graph
        self.width=self.targetRect[1]-self.targetRect[0]+1
        self.height=self.targetRect[3]-self.targetRect[2]+1

        #Fills graphData with values from 'graph'
        graphData = self.convertGraph(self,graph)

        #return only x,maxy and skip any none values
        data = self.cleanData(self,graphData)
        return data

def ConvertQuery(l,comparisonTypes):
    '''for each type being processed, convert the query
    and add the result to a dictionary to be returned'''
    queryDict={}
    for cType in comparisonTypes:
        queryDict[cType]=[]
        queryDict[cType]+=Convert(l,cType)
    return queryDict

class Convert():
    '''
    takes the raw data and converts it into a format that can be compared later
    '''
    def __new__(self,raw,cType):
        f=open('debugConvert.txt','w+')
        f.write(str([cType,raw[:10]]))
        f.close()
        if "raw" == cType:
            #if the comparison is by raw
            return raw
        elif '.' in cType:
            #convert the raw data into the appropiate format
            if cType.split('.')[0] == "AbsoluteROC":
                return self.AbsoluteROC(self,raw,int(cType.split('.')[-1]))
        raise ValueError("Convert type not found: "+str(cType))

    def AbsoluteROC(self,raw,scanrange):
        '''
        The absolute value of the slope of the curve in list l
        Note: this method may not be useful for matching compounds
        '''
        retlist=[]
        for i in range(len(raw)-1):
            retlist+=[(raw[i][0], abs(raw[i+1][1]-raw[i][1]) )]

        return self.Cumulative(self,retlist,scanrange)

class Compare():
    '''
    Compares a query to a subject in the database
    Converts the subject first if needed
    '''
    def __new__(self,cType,subject,query):
        f=open('debug.txt','w+')
        f.write(str([cType,subject[:10],query[:10]]))
        f.close()
        if not "raw" in cType or "raw" == cType:
            #if the subject doesn't need to be converted
            return self.CorCompare(self,subject,query)
        elif "." in cType:
            #else the subject need to be converted
            if cType.split('.')[0] in ["AbsoluteROC"]:
                return self.CorCompare(self, Convert(subject,cType) ,query)
        raise ValueError("Compare type not found: "+str(cType))

    def CorCompare(self,raw1,raw2):
        
        avg1=sum([r[1] for r in raw1])/len(raw1)
        avg2=sum([r[1] for r in raw2])/len(raw2)

        A=0
        B=0
        C=0
        
        x2=0
        x1=0
        while x1<len(raw1) and x2<len(raw2):
            if abs(raw1[x1][0]-raw2[x2][0]) < 7:
                A+=(raw1[x1][1]-avg1)*(raw2[x2][1]-avg2)
                B+=(raw1[x1][1]-avg1)**2
                C+=(raw2[x2][1]-avg2)**2
                
                x2+=1
                x1+=1
            else:
                if raw1[x1][0]>raw2[x2][0]:
                    x2+=1
                elif raw1[x1][0]<raw2[x2][0]:
                    x1+=1

        return 999-999*(A/((B*C)**.5)+1)/2

def AddSortResults(differenceDict,casNums):
    '''
    Take a dictionary with casNums as keys filled with dictionaries with types as keys
    Add the differences of the types for each casnum together and return a sorted list
    where each element is the compound's difference from the query followed by the casnum
    '''
    comparisonTypes=list(differenceDict.keys())[:]

    differenceList=[]
    for i in range(len(casNums)):
        dif=0
        for cType in comparisonTypes:
            if differenceDict[cType][i]:
                dif+=differenceDict[cType][i][0]
        differenceList+=[(dif,differenceDict[cType][i][1])]
    differenceList.sort()

    return differenceList

def SmartSortResults(differenceDict,casNums):
    '''
    Take a dictionary with casNums as keys filled with dictionaries with types as keys
    Return a sorted list where each element is the compound's difference from
    the query followed by the casnum
    
    The compounds are sorted by first seperating compounds by type and then sorting each list
    Each list adds its top result to the bestDict, then any compounds that have been paced
    in the bestDict by the majority of the comparison types are added to the bottom of the
    difference list
    '''
    comparisonTypes=list(differenceDict.keys())[:]

    for cType in comparisonTypes:
        differenceDict[cType].sort()
    differenceList=[]

    bestDict={}
    for i in range(len(casNums)):#casNum
        bestDict[casNums[i]]=[]

    for i in range(len(casNums)):
        tempList=[]
        for cType in comparisonTypes:
            if differenceDict[cType][i]!=(0,):#not found due to active update
                if bestDict[differenceDict[cType][i][1]]!="Done":
                    bestDict[differenceDict[cType][i][1]]+=[(differenceDict[cType][i][0],cType)]   
        for casNum in list(bestDict.keys()):
            if bestDict[casNum]!="Done":
                if len(bestDict[casNum])>=max(1,len(comparisonTypes)//2+1):
                    dif=0
                    for comp in bestDict[casNum]:
                        dif=max(dif,comp[0])
                    tempList+=[(dif,casNum)]
                    bestDict[casNum]="Done"
        if tempList:
            tempList.sort()
            differenceList+=tempList

    return differenceList

class IRDB:
    def __init__(self):
        self.conn = sqlite3.connect(os.path.realpath("IR.db"))
        self.cur = self.conn.cursor()

    def searchIRDB(self, sqlQuery):
        self.cur.execute(sqlQuery)
        return self.cur.fetchall()

    def writeIRDB(self, sqlWrite, dbValues=None):
        try:
            if dbValues:
                self.cur.execute(sqlWrite, dbValues)
            else:
                self.cur.execute(sqlWrite)
            return True
        except Exception as e:
            return False

    def commitIRDB(self):
        try:
            self.conn.commit()
            return True
        except Exception as e:
            return False

    def fetchallIRDB(self):
        return self.cur.fetchall()
#------------------------------------------------------------------------------
