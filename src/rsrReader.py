""" 
DESCRIPTION                                                                                                           

The script reads the open-loop RSR radio science data file and fills out the velocity template. 

INPUTS:                                                                                                         
    1. template      : Path to the reference velocity template.
    2. inputFile     : Path to the input data file.
    3. outFile       : Path to the output file.
    4. nameNoExe     : File name with no extension.
    5. creTime       : Creation or current time.
    6. size          : Size of the input file.
    7. verbose       : Boolean, displays the progress bar if set to True.

RETURN:
    1. Filled velocity template.

AUTHOR  
        
   Ashok Verma (ashokverma@ucla.edu)
   Department of Earth, Planetary, and Space Sciences
   University of California, Los Angeles

"""
#-------------------------------------------------------------------------------------------------------------------
import os, sys, traceback, optparse, ast
import struct, glob, shutil
import subprocess as sp
import hashlib
from time import localtime, strftime 
import datetime as dt
import functions as fn
from astropy.time import Time
import numpy as np
from dsnContext import dsnContext

#-------------------------------------------------------------------------------------------------------------------
class byteFmt:
      # define data type format
      # Ns  == string of N bytes
      byteOrder = '>' # > == network (big-endian) byte order
      uint1     = "B" # B == 1 byte unsigned integer (b==signed) 
      uint2     = "H" # H == 2 byte unsigned integer (h==signed)
      uint4     = "I" # I == 4 byte unsigned integer (i==signed)
      uint8     = "Q" # Q == 8 byte unsigned long long (q==signed)
      float4    = "f" # f == 4 byte float
      float8    = "d" # d == 8 byte float (double)
      
#-------------------------------------------------------------------------------------------------------------------
def getDate (y, doy, sec):

    date = (dt.datetime(y, 1, 1) + dt.timedelta(days=doy - 1, seconds=sec))
    date = date.strftime('%Y-%m-%dT%H:%M:%S')

    return date

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
       fn.colorTxt("\nError while formating\n", 'lightred')
       sys.exit()

    return (fieldNames, format)

#-------------------------------------------------------------------------------------------------------------------
def unpack( fields, format, data, offset=0):

    # Unpack the data.
    elems = format.unpack_from( data, offset )
   
    obj = {}
    for i in range( len( fields ) ):
        obj.update({fields[i]:elems[i]})

    return obj

#-------------------------------------------------------------------------------------------------------------------
def getHeader (data, start):
    byteLength = 20
    end = start + byteLength

    setup = [
       ( "control_auth_id", "4s" ),
       ( "sfdu_version_id", "1s" ),
       ( "sfdu_class_id", "1s" ),
       ( "reserve2", "2s" ),
       ( "data_description_id", "4s" ),
       ( "sfdu_length", byteFmt.uint8 ),
       ]

    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    if obj["control_auth_id"] != b'NJPL':
       errMsg = "Error while reading SFDU Label: \n"\
                "Expected control_auth_id : NJPL; got %s instead \n"
       fn.raiseErr (errMsg %(obj["control_auth_id"]))

    if obj["data_description_id"] != b'C997':
       errMsg = "Error while reading SFDU Label: \n"\
                "Expected data_description_id : C997; got %s instead \n"
       fn.raiseErr (errMsg %(obj["data_description_id"][1:]))


    return obj, end 
    
#-------------------------------------------------------------------------------------------------------------------
def getAggregate (data, start):
    byteLength = 4 
    end = start + byteLength

    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),                   
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    if obj["chdo_type"] != 1:
       errMsg = "Error while reading Aggregation CHDO Label: \n"\
                "Expected chdo_type : 1; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_type"]))                     
    
    if obj["chdo_length"] != 232:
       errMsg = "Error while reading Aggregation CHDO Label: \n"\
                "Expected chdo_length : 232; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_length"]))                     


    return obj, end 

#-------------------------------------------------------------------------------------------------------------------
def getPrimary (data, start):
    byteLength = 8 
    end = start + byteLength

    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),
       ( "mjr_data_class", byteFmt.uint1 ),
       ( "mnr_data_class", byteFmt.uint1 ),
       ( "mission_id", byteFmt.uint1 ),
       ( "format_code", byteFmt.uint1 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    if obj["chdo_type"] != 2:
       errMsg = "Error while reading Primary Header CHDO: \n"\
                "Expected chdo_type : 2; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_type"]))                     
    
    return obj, end 

#-------------------------------------------------------------------------------------------------------------------
def getSecondary (data, start ):
    byteLength = 56
    end = start + 224
    success = False

    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),
       ( "not_in_use", "7s" ),
       ( "rcvr", byteFmt.uint1 ),
       ( "not_in_use", "32s" ),
       ( "year", byteFmt.uint2 ),
       ( "doy", byteFmt.uint2 ),
       ( "sod", byteFmt.float8 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    if obj["chdo_type"] != 104:
       errMsg = "Error while reading Secondary Header CHDO: \n"\
                "Expected chdo_type : 104; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_type"]))                     
  
    if obj["chdo_length"] != 220:
       errMsg = "Error while reading Secondary Header CHDO: \n"\
                "Expected chdo_length : 220; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_length"]))                     

    date = getDate(obj["year"], obj["doy"], obj["sod"])
    rcvr = obj["rcvr"]

    return obj, end, date, rcvr

#-------------------------------------------------------------------------------------------------------------------
def getDataCHOD (data, start ):
    byteLength = 4 

    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    if obj["chdo_type"] != 10:
       errMsg = "Error while reading DATA CHOD: \n"\
                "Expected CHOD Type : 10; got %s instead \n"
       fn.raiseErr (errMsg %(obj["chdo_type"]))                

    return obj


#-------------------------------------------------------------------------------------------------------------------
def readFile(data):
     # read input file
     with open(data, "rb") as f: byte = f.read()
     f.close()
     return byte

#-------------------------------------------------------------------------------------------------------------------
def findStartByte(data):
    
    # find start byte
    dataLen = len( data )
    start = 0
    while start < dataLen:
          i = data.find(b"N")
          sfdu_class_id = struct.unpack('s', data[i+5:i+6])[0]
          if sfdu_class_id == b"I":
             startByte = i
             return startByte
          elif sfdu_class_id == b"K":
             i = data.find(b"CCSD$$MARKER$T-2-34$", i)
             startByte = i + 20 + 20
             return startByte
          else:
             print("The file '%s' doesn't contain the appropriate header" %(inputFile))
             sys.exit()

#----------------------------------------------------------------------------------------------------------------

def docLid (tt):
    tt = Time(tt)
    if tt < Time("2008-02-29"): return "dsn.0159-science:2001-02-28"
    if tt >= Time("2008-02-29") and tt < Time("2009-05-18"): return "dsn.0159-science:2008-02-29"
    if tt >= Time("2009-05-18"): return "dsn.0159-science:2009-05-18"

#----------------------------------------------------------------------------------------------------------------
def fillFileArea (inputFile, velTemp, nameNoExe, dsn, start, stop, creTime, size, recod, lenth, md5, rep, grpLen):

    year = int(start[:4])
    _dsn =  list(map(int, dsn))
    _dsn = str(_dsn)[1:-1] 
    
    readFile = open(velTemp, "r")
    fileData = readFile.read()
    readFile.close()

    newArea = fileData.replace("FILE_NAME", inputFile.split("/")[-1])
    newArea = newArea.replace("FILEWITHNOEXE", nameNoExe.lower())
    newArea = newArea.replace("START_TIME", start)
    newArea = newArea.replace("STOP_TIME", stop)
    newArea = newArea.replace("CREATION_DAY", creTime[:10])
    newArea = newArea.replace("CREATION_TIME", creTime[:10])
    newArea = newArea.replace("FILE_SIZE", str(size))
    newArea = newArea.replace("RECORDS", str(recod))
    newArea = newArea.replace("LENGTH", str(lenth))
    newArea = newArea.replace("MD5", str(md5))
    newArea = newArea.replace("DSSUSED", str(_dsn))
    newArea = newArea.replace("REPETITION", str(rep))
    newArea = newArea.replace("GROUP", str(grpLen))
    newArea = newArea.replace("DSNOBSSYS", dsnContext(dsn, year, inputFile))
    newArea = newArea.replace("DOCLID", "urn:nasa:pds:radiosci.documentation:"+docLid(start))

    tempFile = "_temp_%s.xml"%(nameNoExe)
    outFile = open("/tmp/"+tempFile, "w")
    outFile.write(newArea)
    outFile.close()

    return tempFile

#----------------------------------------------------------------------------------------------------------------
def fixedPartialRecords (chunkSize, size, frac, inputFile, outputFile):

    msg = "WARNING:\n \
    Partial record found at the end of file: \n \
    File Name: %s \n \
    Partial records: %s Byte(s) \n \
    Fix partial record: %s\n"

    endByte = int(frac)*chunkSize
    partialRecords = size-endByte

    fn.colorTxt (msg%(inputFile, partialRecords, fixParRecords), 'purple')

    endByte = endByte-chunkSize
    with open(outputFile, 'rb+') as f:
         f.seek(0,2)
         f.truncate(size-partialRecords)
    f.close()         

    return endByte 

#-------------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def getTemplet(velTemp, inputFile, outFile, nameNoExe, creTime, size, verbose):
    
    # setup update bar, if requested
    if verbose:
       prog = fn.progressBar()                                  
       prog.start(1, "Reading and unpacking RSR file")

    # read input file
    byte=readFile(inputFile)
    size = len(byte)
    
    # initialize arrays and variable
    outfile = [0]*18 
    nsfdu = [0]*18  
    secPast = []
    orec = 0
    
    # get start byte
    start = findStartByte(byte)
    
    # start loop
    i = 0
    coverage = []
    dss = []
    while start < size:
          
          # read tracking SFDU label 
          objH, end = getHeader(byte, start)

          # read aggregation CHDO label
          objA, end = getAggregate (byte, end)
    
          # read primary CHDO data
          objP, end = getPrimary(byte, end)
    
          # read secondary CHDO data
          objS, end, date, rcvr = getSecondary (byte, end)
          coverage.append(date)
          dss.append(str(rcvr))

          # read data CHOD
          objD = getDataCHOD (byte, end)
    
          # update chunkSize
          chunkSize = objH['sfdu_length'] + 20

          # getting last chunk
          i = i + 1
          if i == 2: break
          frac = size/chunkSize

          if not frac.is_integer():
             start = fixedPartialRecords (chunkSize, size, frac, inputFile, outFile) 
          else:
             start = size - chunkSize
          
          if verbose: prog.setStep(start*1.0/size)
          
    #  get station    
    dsn = np.unique(dss)
    prog.stop()

    # fill template
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Updating velocity template")
    
    md5 = fn.calculateMD5(outFile)
    tempOut=fillFileArea (outFile, velTemp, nameNoExe, dsn, coverage[0], coverage[1], creTime, int(size), 
    int(size/chunkSize), int(chunkSize), md5, int(objD['chdo_length']/4), int(objD['chdo_length']))

    # empty memory
    del byte
    if verbose: prog.stop()
    return tempOut
    
#-------------------------------------------------------------------------------------------------------------------
