# -*- coding: utf-8 -*-                                                                                                            
""" 
DESCRIPTION

The script reads the ODF (TRK-2-18) radio science data file and fills out the velocity template. 

INPUTS:                                                                                                         
    1. template      : Path to the reference velocity template.
    2. inputFile     : Path to the input data file.
    3. nameNoExe     : File name with no extension.
    4. creTime       : Creation or current time.
    5. size          : Size of the input file.
    6. verbose       : Boolean, displays the progress bar if set to True.

RETURN:
    1. Filled velocity template.

AUTHOR  
        
   Ashok Verma (ashokverma@ucla.edu)
   Department of Earth, Planetary, and Space Sciences
   University of California, Los Angeles

"""
#-------------------------------------------------------------------------------------------------------------------
import os, sys, traceback, optparse, ast
import struct, shutil
import subprocess as sp
import collections
from time import localtime, strftime
from Data import Data
import datetime as dt
from astropy.time import Time
import numpy as np
import functions as fn
from dsnContext import dsnContext

#-------------------------------------------------------------------------------------------------------------------
class bcolors:
    '''Colors class:
    reset all colors with colors.reset
    two subclasses fg for foreground and bg for background.
    use as colors.subclass.colorname.
    i.e. colors.fg.red or colors.bg.green
    also, the generic bold, disable, underline, reverse, strikethrough,
    and invisible work with the main class
    i.e. colors.bold
    '''
    reset='\033[0m'
    bold='\033[01m'
    disable='\033[02m'
    underline='\033[04m'
    reverse='\033[07m'
    strikethrough='\033[09m'
    invisible='\033[08m'
    class fg:
        black='\033[30m'
        red='\033[31m'
        green='\033[32m'
        orange='\033[33m'
        blue='\033[34m'
        purple='\033[35m'
        cyan='\033[36m'
        lightgrey='\033[37m'
        darkgrey='\033[90m'
        lightred='\033[91m'
        lightgreen='\033[92m'
        yellow='\033[93m'
        lightblue='\033[94m'
        pink='\033[95m'
        lightcyan='\033[96m'
    class bg:
        black='\033[40m'
        red='\033[41m'
        green='\033[42m'
        orange='\033[43m'
        blue='\033[44m'
        purple='\033[45m'
        cyan='\033[46m'
        lightgrey='\033[47m'


#-------------------------------------------------------------------------------------------------------------------
class byteFmt:
      # define data type format
      # Ns  == string of N bytes
      byteOrder = '>' # > == network (big-endian) byte order
      uint1     = "B" # B == 1 byte unsigned integer (b==signed) 
      uint2     = "H" # H == 2 byte unsigned integer (h==signed)
      int4      = "I" # I == 4 byte unsigned integer (i==signed)
      int4      = "i" # i == 4 byte signed integer (i==signed)
      uint8     = "Q" # Q == 8 byte unsigned long long (q==signed)
      float4    = "f" # f == 4 byte float
      float8    = "d" # d == 8 byte float (double)


#-------------------------------------------------------------------------------------------------------------------
def readFile(data):                                                                                                                      
     # read input file
     with open(data, "rb") as f: byte = f.read()
     f.close()
     return byte

#-------------------------------------------------------------------------------------------------------------------                 
def getPastSec2Date(sec):

    # convert second past from J1950 to calendar time     
    time = str(dt.datetime(1950, 1, 1 ) + dt.timedelta(seconds=sec))
    time = time.split()[0]+"T"+time.split()[1]

    return time

#-------------------------------------------------------------------------------------------------------------------
def getFormat(setup, byteLength):
    # define format of the data chunk
    dataFormat = byteFmt.byteOrder[:]

    fieldNames = []
    for f in setup:
        dataFormat += f[1]
        fieldNames.append( f[0] )
    
    format = struct.Struct( dataFormat )

    if format.size != byteLength:
       print("Error while formating")
       sys.exit()

    return (fieldNames, format)

#-------------------------------------------------------------------------------------------------------------------
def getStat(fileName):
    # get outfile file's statistics
    byte = readFile(fileName)
    size = len(byte)
    obj, end = getHeader(byte, 0)
    rcrdLen = obj["sfdu_length"]+20
    rcrdSize = round(size/rcrdLen)

    return rcrdLen, rcrdSize

#-------------------------------------------------------------------------------------------------------------------
def dataType(i):
    # data type format ids
    formatIds = {
         11 : "DSN One-Way Doppler",    
         12 : "DSN Two-Way Doppler",    
         13 : "DSN Three-Way Doppler",    
         37 : "DSN Two-Way Range",    
         }
     
    return formatIds[i]

#-------------------------------------------------------------------------------------------------------------------
def unpack( fields, format, data, offset=0):
    # Unpack the data.
    elems = format.unpack_from( data, offset )
   
    obj = {}
    for i in range( len( fields ) ):
        obj.update({fields[i]:elems[i]})

    return obj


#-------------------------------------------------------------------------------------------------------------------
def TableHeader (data, start, checkPrimaryKeys=True):
    byteLength = 16
    end = start + byteLength  + 20

    setup = [
       ( "Primary_Key", byteFmt.int4 ),
       ( "Secondary_Key", byteFmt.int4 ),
       ( "Logical_Record_Length", byteFmt.int4 ),
       ( "Group_Start_Packet_Number", byteFmt.int4 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    primaryKeys = [101, 107, 109, 2030, 2040, 105, -1]

    if checkPrimaryKeys:
       if obj['Primary_Key'] not in primaryKeys:
          err = "Unexpected Primary key: %s \n" \
                "Expeted Primary keys  : %s \n" 
          fn.raiseErr (err %(obj['Primary_Key'], primaryKeys))
           
    return obj, end

#-------------------------------------------------------------------------------------------------------------------
def getSCID (data, start):
    byteLength = 4 
    start = start + 16
    end = start + byteLength  + 16


    setup = [("scID", byteFmt.int4)]

    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    return obj['scID'], end
#-------------------------------------------------------------------------------------------------------------------
def getBits (items):
    itms = []
    for i in range(len(items)):
        if sys.version_info[0] == 2:
           byte = bin(ord(items[i]))[2:].rjust(8,'0')
        else:      
           byte = format(items[i], 'b').rjust(8,'0')

        for bit in byte:
            itms += [str(bit)]

    r  = ("".join(itms))
    del itms
    return r

#-------------------------------------------------------------------------------------------------------------------
def twos(val, sign=False):
    if sign:
       """compute the 2's complement of int value val"""
       bits = len(val)
       val = int(val, 2)
       if (val & (1 << (bits - 1))) != 0:
          val = val - (1 << bits)
    else: val = int(val,2)

    return val


def getItems2To3 (items):
    r = getBits (items)
    obj = {}
    obj.update({"Item02": twos(r[00:10])})
    obj.update({"Item03": twos(r[10:32])})

    return obj

def getItems6To19 (items, start):
    r = getBits (items)
    obj = {}
    obj.update({"Item06": twos(r[00:3])})
    obj.update({"Item07": twos(r[3:10])})
    obj.update({"Item08": twos(r[10:17])})
    obj.update({"Item09": twos(r[17:19])})
    obj.update({"Item10": twos(r[19:25])})
    obj.update({"Item11": twos(r[25:27])})
    obj.update({"Item12": twos(r[27:29])})
    obj.update({"Item13": twos(r[29:31])})
    obj.update({"Item14": twos(r[31:32])})
    obj.update({"Item15": twos(r[32:39])})
    obj.update({"Item16": twos(r[39:49])})
    obj.update({"Item17": twos(r[49:50])})
    obj.update({"Item18": twos(r[50:72])})
    obj.update({"Item19": twos(r[72:96])})

    msg = "\nUnexpected Format ID: %s\n"\
          "Expected Format ID  : 2\n\n"\
          "If Format ID: 1,  \n" \
          "The ODF was created on or before 1997-04-14 \n"\
          "and will not be accurately described by \n"\
          "this script. \n \n"\
          "If Format ID: 2, see: \n"\
          "JPL/DSN Document 820-13; Rev A \n"\
          "   DSN System Requirements \n"\
          "   Detail Interface Design \n"\
          "          TRK-2-18  \n"\
          "DSN Tracking System Interfaces\n"\
          "  Orbit Data File Interface \n"\
          "          Mark IVA \n"\
          " Effective Date: May 15, 1984 \n"
    if obj[ "Item06"] != 2 and start <= 180:
       fn.raiseErr (msg%(obj[ "Item06"]))

    return obj

def getItems20To22 (items):
    r = getBits (items)
    obj = {}
    obj.update({"Item20": twos(r[00:20], True)})
    obj.update({"Item21": twos(r[20:42])})
    obj.update({"Item22": twos(r[42:64])})

    return obj


def getObservable (obj):

    # SPACECRAFT ID
    scID    = obj["Item16"]

    # OBSERVATION TIME FROM 1950 
    timeTag = obj["Item01"] + obj["Item02"]*1e-3 

    data = Data (scID=scID, timeTag=timeTag, rcvr=str(obj["Item07"]), groupStart=int(obj["Item05"])) 

    return data 


def TableOrbitData (data, start):

    byteLength = 36
    end = start + byteLength

    setup = [

       ( "Item01", byteFmt.int4 ),
       ( "Item 2-3", '4s' ),
       ( "Item04", byteFmt.int4 ),
       ( "Item05", byteFmt.int4 ),
       ( "Item 6-19", "12s" ),
       ( "Item 20-22", "8s" ),

       ]
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start)

    objs = {}
    objs.update({"Item01":obj["Item01"]})
    objs.update(getItems2To3 (obj['Item 2-3']))
    objs.update({"Item04":obj["Item04"]})
    objs.update({"Item05":obj["Item05"]})
    objs.update(getItems6To19(obj['Item 6-19'], start))
    objs.update(getItems20To22(obj['Item 20-22']))
    Observable = getObservable(objs)

    return Observable, end
    
#-------------------------------------------------------------------------------------------------------------------
def getRampItems5To6 (items):
    r = getBits (items)
    obj = {}
    obj.update({"Item05": twos(r[00:22])})
    obj.update({"Item06": twos(r[22:32])})

    return obj

def TableRampData (data, start):
    
    byteLength = 36
    end = start + byteLength
    
    setup = [
       ( "Item01", byteFmt.int4 ),
       ( "Item02", byteFmt.int4 ),
       ( "Item03",  byteFmt.int4 ),
       ( "Item04",  byteFmt.int4 ),
       ( "Item 5-6", "4s" ),
       ( "Item07", byteFmt.int4 ),
       ( "Item08", byteFmt.int4 ),
       ( "Item09", byteFmt.int4 ),
       ( "Item10", byteFmt.int4 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    objs = unpack( _fields, _format, data, start)

    objs.update(getRampItems5To6 (objs['Item 5-6']))
    rampData = Data(stnID = objs["Item06"])

    return rampData, end

#-------------------------------------------------------------------------------------------------------------------
def getByteRecs (byteSum, table, start, end):
    recLen = int((end-start)/36)
    if recLen != 0: byteSum.write("%20s, %20s, %20s \n" %(table, start, recLen))

#-------------------------------------------------------------------------------------------------------------------
def writeVelocityTemp (template, inputFile, nameNoExe, creTime, md5, size, dsn):
    # get size, md5, and creation date
    inFile = os.path.basename(inputFile) 

    # filled and generate the velocity templete
    newTab, start, stop = getTabelXML (os.path.dirname(template), "/tmp/_temp_%s_odf.sum"%(nameNoExe))

    # filled the file Area
    fileArea = fillFileArea (template, inFile, nameNoExe, creTime, start, stop, size, md5, dsn)

    # write velocity templete
    xml = fileArea.strip()+"\n             "+newTab.strip()+"\n      </File_Area_Observational>\n"+"</Product_Observational>"
    tempFile = "_temp_%s.xml"%(nameNoExe)
    outVM = open("/tmp/"+tempFile, "w")
    outVM.write(xml)
    outVM.close()

    return tempFile 

#-------------------------------------------------------------------------------------------------------------------
def docLid (tt):
    tt = Time(tt)
    if tt < Time("1988-10-15"): return "dsn.trk-2-18:1988-01-15"
    if tt >= Time("1988-10-15") and tt < Time("1996-08-15"): return "dsn.trk-2-18:1988-10-15"
    if tt >= Time("1996-08-15") and tt < Time("2000-06-15"): return "dsn.trk-2-18:1996-08-15"
    if tt >= Time("2000-06-15") and tt < Time("2008-02-29"): return "dsn.trk-2-18:2000-06-15"
    if tt >= Time("2008-02-29"): return "dsn.trk-2-18:2008-02-29"

#----------------------------------------------------------------------------------------------------------------

def fillFileArea (template, inputFile, nameNoExe, creTime, start, stop, size, md5, dsn):

    year = int(start[:4])
    _dsn =  list(map(int, dsn))
    _dsn = str(_dsn)[1:-1] 
    
    readFile = open(template, "r")
    fileData = readFile.read()
    readFile.close()

    newArea = fileData.replace("FILEWITHNOEXE", nameNoExe)
    newArea = newArea.replace("CREATION_DAY", creTime[:10])
    newArea = newArea.replace("START_TIME", start)
    newArea = newArea.replace("STOP_TIME", stop)
    newArea = newArea.replace("FILE_NAME", inputFile)
    newArea = newArea.replace("CREATION_TIME", creTime[:10])
    newArea = newArea.replace("FILE_SIZE", str(size))
    newArea = newArea.replace("MD5", md5)
    newArea = newArea.replace("DSSUSED", str(_dsn))
    newArea = newArea.replace("DSNOBSSYS", dsnContext(dsn, year, inputFile))
    newArea = newArea.replace("DOCLID", "urn:nasa:pds:radiosci.documentation:"+docLid(start))

    return newArea

#----------------------------------------------------------------------------------------------------------------
def fillTable(templeteDir, table, offSet, rcrdSize, stn):

    fileName = templeteDir+"/"+table+".xml"
    readTab = open(fileName, "r")
    tabData = readTab.read()
    readTab.close()

    newTab = tabData.replace("STNID", str(stn).strip())
    newTab = newTab.replace("OFFSET", str(offSet).strip())
    newTab = newTab.replace("RECORDS", str(rcrdSize).strip())

    return newTab

#-------------------------------------------------------------------------------------------------------------------
def getTabelXML (templeteDir, fileName):
    file = open(fileName, "r")
    newTab = ""
    stn = "NAN"
    for l in file:
        if "#" not in l:
            getAll = l.strip().split(",")
            t = getAll[0]; o = getAll[1]; r = getAll[2]
            if "Table" in t:
               if "Table-3-4" in t: stn = t[-2:]; t = t[:-2]
               newTab = newTab + fillTable(templeteDir, t, o, r, stn)
            else:
               start = o; stop = r

    return newTab, start.strip(), stop.strip()


#-------------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def getTemplet(template, inputFile, nameNoExe, creTime, size, verbose):
    

    # setup update bar, if requested
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Reading and unpacking ODF file")


    # odf summary file
    sumFile = open("/tmp/_temp_%s_odf.sum"%(nameNoExe), "w")
    sumFile.write("%20s, %20s, %20s \n" %("# Table", "Offset", "Record Length"))

    # read input file 
    byte = readFile(inputFile)
    size = len(byte)

    timeAll = []
    start = 0
    # Table 3-1a. ODF File Label Group Header Format
    obj, end = TableHeader (byte, start)
    if obj['Primary_Key'] == 101: 
       getByteRecs(sumFile, "Table-3-1a", start, end)
       # Table 3-1b. ODF File Label Group Data Format
       scID, start = getSCID (byte, end)
       getByteRecs(sumFile, "Table-3-1b", end, start)
    
    # Table 3-2a. ODF Identifier Group Header Format
    obj, end = TableHeader (byte, start)
    if obj['Primary_Key'] == 107: 
       getByteRecs(sumFile, "Table-3-2a", start, end)
       # Table 3-2b. ODF Identifier Group Data Format
       start = end + 36
       getByteRecs(sumFile, "Table-3-2b", end, start)
    
    # Table 3-3a. ODF Orbit Data Group Header Format
    obj, end = TableHeader (byte, start)
    dss = []
    if obj['Primary_Key'] == 109: 
       getByteRecs(sumFile, "Table-3-3a", start, end)
       start = end
       # Table 3-3b. ODF Orbit Data Group Data Format
       while start < size:
             Observable, start = TableOrbitData (byte, start)
             if scID != Observable.scID and Observable.groupStart+1 == int(start/36) :
                start = start - 36
                getByteRecs(sumFile, "Table-3-3b", end, start)
                break 
             else:
                if scID == Observable.scID:             
                   timeAll.append(Observable.timeTag)                
                   dss.append(Observable.rcvr)
                   if verbose: prog.setStep(0.95*start/size)
                else: 
                   msg = "ERROR: While attempting to read Orbit Data Group \n"\
                         "File Name: %s \n"\
                         "Expected Spacecraft ID: %s; instead got %s \n"
                   fn.raiseErr(msg%(inputFile, scID, Observable.scID))                 
    
    # Table 3-4a. ODF Ramp Groups Header Format
    for i in range(100):
        obj, end = TableHeader (byte, start)
        if obj['Primary_Key'] == 2030:
           stnID = obj['Secondary_Key']
           getByteRecs(sumFile, "Table-3-4a"+str(stnID), start, end)
           start = end
           # Table 3-4b. ODF Ramp Groups Data Format
           while start < size:
                rampData, start = TableRampData (byte, start)
                if rampData.stnID != stnID:
                   start = start - 36
                   getByteRecs(sumFile, "Table-3-4b"+str(stnID), end, start)
                   break
                else: 
                   if verbose: prog.setStep(0.95*start/size)    
        else: break 
    
    # Table 3-5a. ODF Clock Offsets Group Header Format
    obj, end = TableHeader (byte, start)
    if obj['Primary_Key'] == 2040:
       # Table 3-5b. ODF Clock Offsets Group Data Format
       getByteRecs(sumFile, "Table-3-5a", start, end)
       start = end + 36
       getByteRecs(sumFile, "Table-3-5b", end, start)
       if verbose: prog.setStep(0.95*start/size) 
    
    # Table 3-7a. ODF Data Summary Group Header Format
    obj, end = TableHeader (byte, start)
    if obj['Primary_Key'] == 105:
       # Table 3-7b. ODF Data Summary Group Data Formatt
       getByteRecs(sumFile, "Table-3-7a", start, end)
       start = end 
       while start < size:
             obj, start = TableHeader (byte, start, False)           
             if obj['Primary_Key'] == -1:
                start = start - 36                 
                getByteRecs(sumFile, "Table-3-7b", end, start)                
                break        
             else: 
                if verbose: prog.setStep(0.95*start/size)                 
    
    # Table 3-8. ODF End-of-File Group Format
    obj, end = TableHeader (byte, start)
    if obj['Primary_Key'] == -1:
       getByteRecs(sumFile, "Table-3-8a", start, end)
       if verbose: prog.setStep(0.95*start/size) 
   
    # get start and stop time of the output file
    sumFile.write(("%20s, %20s, %20s \n")%("Start-Stop Time",  getPastSec2Date(min(timeAll)), getPastSec2Date(max(timeAll))))
    sumFile.close()
    del timeAll, byte

    # get stations
    dsn = np.unique(dss)
    if verbose: prog.stop()

    # update template 
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Updating velocity template")
    
    md5 = fn.calculateMD5(inputFile)
    tempOut=writeVelocityTemp (template, inputFile, nameNoExe, creTime, md5, size, dsn)

    # remove tempOut file
    os.remove("/tmp/_temp_%s_odf.sum"%(nameNoExe))
    if verbose: prog.stop()

    return tempOut

#-------------------------------------------------------------------------------------------------------------------
