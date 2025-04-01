# -*- coding: utf-8 -*-                                                                                                            
""" 
 DESCRIPTION                                                                                                           
 
 The script reads the radio science Archival Tracking Data File (ATDF) and fills out the velocity template. 

 INPUTS:                                                                                                         
     1. template      : Path to the reference velocity template.
     2. inputFile     : Path to the input data file.
     3. nameNoExe     : File name with no extension.
     4. creTime       : Creation or current time.
     5. dsn           : List of Deep Space Network (dsn).
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
import sys, os, traceback, optparse, ast
import struct, shutil
import numpy as np
import subprocess as sp
import collections
from Data import Data
from time import localtime, strftime
import datetime as dt
from astropy.time import Time
from dsnContext import dsnContext
import functions as fn

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
def getDate (year, doy, hr, mn, sec):

    date = (dt.datetime(year, 1, 1) + dt.timedelta(days=doy - 1, hours=hr, minutes=mn, seconds=sec))
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
       fn.raiseErr("Error while formating")

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
def readFile(data):
     # read input file
     with open(data, "rb") as f: byte = f.read()
     f.close()
     return byte

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
def getFLoat (hp, lp):
    return hp*10000.0 + lp/1000.0   

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


#-------------------------------------------------------------------------------------------------------------------
def TableUnpack (data):
    byteLength =  288 

    setup = [
       
       ( "clm", '288s' ),
       ]
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data)
    
    return obj

#-------------------------------------------------------------------------------------------------------------------
def getOBJ (format, r):
    
    obj = {}
    i = 0
    bs = 0
    for fmt in format:
        f = False
        i = i + 1
        item = "Item"+str("%03d" %i)
        fmt = fmt.strip()
        if fmt[0] == 'S': f=True
        be = bs + int(fmt[1:])

        obj.update({item : twos(r[bs:be], f)})
        bs = be

    return obj

#-------------------------------------------------------------------------------------------------------------------
def getTableItems (items, tableFormat):
    r = getBits (items)
    obj = getOBJ (tableFormat, r)

    return obj

#-------------------------------------------------------------------------------------------------------------------
def bits2string(b=None):
    return ''.join([chr(int(x, 2)) for x in b])

#-------------------------------------------------------------------------------------------------------------------
def getTable1Data (obj):
    items = {}
    if obj["Item002"] != 128:
       msg = "Item #2 in Table 3-1, ATDF File Identification Logical Record Format, \n" \
             "should be 128, got %s instead \n"
       fn.raiseErr (msg %(obj["Item002"]))
    
    atdf = bits2string([bin(obj["Item015"])[2:]])+bits2string([bin(obj["Item016"])[2:]]) + \
           bits2string([bin(obj["Item017"])[2:]])+bits2string([bin(obj["Item018"])[2:]])

    if atdf != "ATDF":
       msg = "Item #15-18 in Table 3-1, ATDF File Identification Logical Record Format, \n" \
              "should be ATDF, got %s instead \n"
       fn.raiseErr (msg %(atdf))              

    return items

#-------------------------------------------------------------------------------------------------------------------
def getTable2Data (obj):
    items = {}
    if obj["Item002"] != 128:
       msg = "Item #2 in Table 3-2, ATDF Transponder Logical Record Format, \n" \
             "should be 128, got %s instead \n"
       fn.raiseErr (msg %(obj["Item002"]))
     
    if obj["Item003"] != 30:
       msg = "Item #2 in Table 3-2, ATDF Transponder Logical Record Format, \n" \
             "should be 30, got %s instead \n"
       fn.raiseErr (msg %(obj["Item003"]))
    
#-------------------------------------------------------------------------------------------------------------------
def docLid (tt):
    tt = Time(tt)
    if tt < Time("1986-01-21"): return "dsn.trk-2-25:1977-05-19"
    if tt >= Time("1986-01-21") and tt < Time("1988-10-15"): return "dsn.trk-2-25:1986-01-21"
    if tt >= Time("1988-10-15") and tt < Time("1996-07-31"): return "dsn.trk-2-25:1988-10-15"
    if tt >= Time("1996-07-31"): return "dsn.trk-2-25:1996-07-31"

#-------------------------------------------------------------------------------------------------------------------
def writeVelocityTemp (template, inputFile, nameNoExe, start, stop, currentDate, md5, size, records, dsn, padding):

    tempFile = "_temp_%s.xml"%(nameNoExe)
    readFile = open(template, "r")
    fileDate = readFile.read()
    readFile.close()
    year = int(start[:4])
    _dsn =  list(map(int, dsn))
    _dsn = str(_dsn)[1:-1] 

    newArea = fileDate.replace("FILEWITHNOEXE", nameNoExe)
    newArea = newArea.replace("CREATION_DAY", currentDate[:10])
    newArea = newArea.replace("START_TIME", start)
    newArea = newArea.replace("STOP_TIME", stop)
    newArea = newArea.replace("FILE_NAME", os.path.basename(inputFile))
    newArea = newArea.replace("FILE_SIZE", str(size))
    newArea = newArea.replace("CREATION_TIME", currentDate)
    newArea = newArea.replace("MD5", md5)
    newArea = newArea.replace("DSSUSED", str(_dsn))
    newArea = newArea.replace("RECORDS", str(records))
    newArea = newArea.replace("PADDINGOFFSET", str(int(size-padding*288.0)))
    newArea = newArea.replace("PADDING", str(padding))
    newArea = newArea.replace("DSNOBSSYS", dsnContext(dsn, year, inputFile))
    newArea = newArea.replace("DOCLID", "urn:nasa:pds:radiosci.documentation:"+docLid(start))
    newArea = newArea.replace("""\t""", "        ")

    outFile=open("/tmp/"+tempFile, "w")
    outFile.write(newArea)
    outFile.close()

    return tempFile 

#-------------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def getTemplet(template, inputFile, nameNoExe, creTime, dsn, size, verbose):

    # setup update bar, if requested
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Reading and unpacking ATDF file")

    # ATDF format
    T3_1 = [
           "I32", "I8 ", "I32", "I12", "I16", "I8 ", "I12", "I8 ", "I12", "I16", "I8 ", "I8", "I8 ", "I12", "I16", "I8 ", "I12", "I8 ", "I16",\
           "I4 ", "I2048 "
           ]
    
    T3_2 = [
           "I32", "I8 ", "I32", "I12", "I16", "I8 ", "I12", "I8 ", "I12", "I16", "I8 ", "I8 ", "I8 ", "I12", "I16", "I8 ", "I12", "I8 ", "I16",\
           "I12", "I24", "I12", "I24", "I28", "I952"
           ]
    
    T3_3 = [
            "I32", "I8 ", "I32", "I12", "I16", "I8 ", "I8 ", "I8 ", "I20", "I10", "I8 ", "I6 ", "I4 ", "I4 ", "I16", "I8 ", "I8 ", "I8 ", "I1 ",\
            "S18", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I6 ", "I6 ", "I4 ", "I32", "I24", "I24", "I24", "I24", "I24", "I24", "I8 ", "I28", "I24",\
            "I24", "I24", "S24", "S24", "I32", "I32", "S32", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24",\
            "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "I24", "S4 ", "S32", "S4 ", "S32",\
            "S18", "S18", "I8 ", "I4 ", "I2 ", "I1 ", "I1 ", "I1 ", "I1 ", "I8 ", "I10", "S18", "S18", "I24", "I24", "I1 ", "I1 ", "I1 ", "I1 ",\
            "I1 ", "I1",  "I1 ", "I1 ", "I1 ", "I4 ", "I1 ", "I10", "I24", "S12", "S4 ", "S32", "S4 ", "S32", "I4 ", "I32", "S22", "I14", "I23",\
            "I1 ", "I1 ", "I1 ", "I10", "I8 ", "S32", "S32", "I4 ", "I32", "I4 ", "I32", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ",\
            "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I1 ", "I28", "I30", "I32", "I32", "I32", "I32", "I32", "I32", "I32", "I32", "I32"
            ]
   
    # read all bytes
    f = open(inputFile, "rb")
    byte = f.read()

    # Unpack Table 3-1 
    obj = TableUnpack (byte[:288])
    obj = getTableItems (obj["clm"], T3_1)
    getTable1Data (obj)

    # Unpack Table 3-2
    obj = TableUnpack (byte[288:576])
    obj = getTableItems (obj["clm"], T3_2)
    getTable2Data (obj)

    # Unpack Table 3-2
    obj = TableUnpack (byte[576:864])
    obj = getTableItems (obj["clm"], T3_3)
    if obj["Item001"] !=8 : 
       msg = "Item #1 in Table 3-3, ATDF Tracking Data Logical Records Format, \n" \
             "should be 8, got %s instead \n"
       fn.raiseErr (msg %(obj["Item001"]))
    else:
        Year = obj["Item004"] + 1900
        stime = getDate (Year, obj["Item005"], obj["Item006"], obj["Item007"], obj["Item008"])

    if verbose: prog.setStep(0.5)
    getPadding = size
    i = 0 
    etime=''
    while getPadding>288:
          obj = TableUnpack (byte[getPadding-288:getPadding])
          obj = getTableItems (obj["clm"], T3_3)
          Item1 = obj["Item001"]
          if Item1 == 8:     
             Year = obj["Item004"] + 1900
             etime = getDate (Year, obj["Item005"], obj["Item006"], obj["Item007"], obj["Item008"])
             break
          i = i + 1
          getPadding = getPadding-288

    if not etime: fn.raiseErr ("Unable to find the end time of the sample. \nMay be a corrupted file!!!")
    
    records = int(size/288.0-i-2)
    padding = int(i)
    if verbose: prog.stop()

    # updating velocity template
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Updating velocity template")
    
    md5 = fn.calculateMD5(inputFile) 
    tempOut=writeVelocityTemp (template, inputFile, nameNoExe, stime, etime, creTime, md5, size, records, dsn, padding)
    prog.stop()

    return tempOut
    
#-------------------------------------------------------------------------------------------------------------------
