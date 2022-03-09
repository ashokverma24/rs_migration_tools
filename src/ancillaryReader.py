""" 
DESCRIPTION

The script reads radio science ancillary data (e.g. ckf, spk, eop, ion, tro)
and fills out the velocity template. 

INPUTS:
    1. dType         : Data type.
    2. velTemp       : Path to the reference velocity template.
    3. inputFile     : Path to the input data file.
    4. nameNoExe     : File name with no extension.
    5. creTime       : Creation or current time.
    6. dsn           : List of Deep Space Network (dsn), if any.
    7. size          : Size of the input file.
    8. verbose       : Boolean, displays the progress bar if set to True.

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
import numpy as np
import functions as fn
import julian
import hashlib
from time import localtime, strftime 
import datetime as dt
from itertools import (takewhile,repeat)  

#-------------------------------------------------------------------------------------------------------------------
def delFiles (wildCard):
    for f in glob.glob(wildCard):
        if os.path.exists(f): os.remove(f)

#-------------------------------------------------------------------------------------------------------------------
def runCommands(commandArgs):
    # run shell commond within python
    proc = sp.Popen(commandArgs,shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()

    return err

#-------------------------------------------------------------------------------------------------------------------
def fillFileArea (inputFile, velTemp, fileNoEXE, startTime, stopTime, creTime, 
                  size, md5, headerSize, records, recordSize):

    tempFile = "/tmp/_temp_%s.xml"%(fileNoEXE)
    os.system('cp %s %s'%(velTemp, tempFile)) 
    readFile = open(tempFile, "r")
    fileDate = readFile.read()
    readFile.close()

    newArea = fileDate.replace("FILENAME", inputFile.split("/")[-1])
    newArea = newArea.replace("FILEWITHNOEXE", fileNoEXE)
    newArea = newArea.replace("STARTTIME", startTime.upper())
    newArea = newArea.replace("STOPTIME", stopTime.upper())
    newArea = newArea.replace("CREATIONDATE", str(creTime))
    newArea = newArea.replace("FILESIZE", str(size))
    newArea = newArea.replace("MD5", str(md5))
    newArea = newArea.replace("HEADER_SIZE", str(headerSize))
    newArea = newArea.replace("RECORDS", str(records))
    newArea = newArea.replace("RECORD_SIZE", str(recordSize))

    outFile=open(tempFile, "w")
    outFile.write(newArea)
    outFile.close()

    return "_temp_%s.xml"%(fileNoEXE)

#----------------------------------------------------------------------------------------------------------------
def dataTypeInfo (velTemp):
    # read velocity templete
    f = open(velTemp, 'r')
    dType = 'nan'
    for l in f:
        if l.startswith('#set ($dataType'):
           dType = l.split("=")[-1].split(")")[0].strip()[1:-1]
           break
    f.close()
    return dType
           
#----------------------------------------------------------------------------------------------------------------
def dos2unix(inp, out):
    content = ''
    outsize = 0
    with open(inp, 'rb') as infile:
      content = infile.read()
    with open(out, 'wb') as output:
      for line in content.splitlines():
        outsize += len(line) + 1
        output.write(line + b'\n')

#----------------------------------------------------------------------------------------------------------------
def spkInfo (fileName):
    # load supporting data
    lsk = "./../spiceData/lsk/naif0012.tls"
    sclk = "./../spiceData/sclk/cas00172.tsc"

    # out and temporary files
    spacitIn = fileName+"_.spacit.in"
    spacitOut= fileName+"_.bsp"
    spacitlog= fileName+"_.log"
    delFiles(fileName+"_.*")

    # dos2unix
    dos2unix (fileName, fileName+".temp")
    shutil.move(fileName+".temp", fileName)

    # set variable
    startTime = None
    stopTime = None

    # remove old outputs
    delFiles(fileName+"_.*")

    # convert transfer file to binary file 
    try:
       spacitCmd = open(spacitIn, "w")
       spacitCmd.write("T\n%s\n%s\nQ\nN"%(fileName, spacitOut))
       spacitCmd.close()
       cmd = "spacit < %s > %s"%(spacitIn, spacitlog)
       os.system(cmd)
    except Exception as e:
      print (e.message, e.args)

    # if successful, dump summary to the log file
    if os.path.exists(spacitOut):
       cmd = "brief -t -utc %s %s >> %s"%(lsk, spacitOut, spacitlog)
       os.system(cmd)
    
    # read log file and extract info
    if os.path.exists(spacitOut):
       log = open(spacitlog, "r")
       lines = np.array([l.strip() for l in log])
       log.close()

       # let's check if dump was successful
       index = np.where(lines =='A traceback follows.  The name of the highest level module is first.')[0]
       if index:
          fn.raiseErr("\n".join(lines))

       # let's read further if successful          
       index = np.where(lines=='Bodies                Start of Interval (UTC)         End of Interval (UTC)')[0][0]+2
       info = lines[index].split()
       startTime = julian.cal_dt(info[2]+"T"+info[3])
       stopTime =  julian.cal_dt(info[4]+"T"+info[5])
   
    delFiles(fileName+"_.*")
    return startTime, stopTime, None, None, None 

#----------------------------------------------------------------------------------------------------------------
def ckfInfo (fileName):
    # load supporting data
    lsk = "./../spiceData/lsk/naif0012.tls"
    sclk = "./../spiceData/sclk/cas00172.tsc"

    # out and temporary files
    spacitIn = fileName+"_.spacit.in"
    spacitOut= fileName+"_.bc"
    spacitlog= fileName+"_.log"
    delFiles(fileName+"_.*")

    # dos2unix
    dos2unix (fileName, fileName+".temp")
    shutil.move(fileName+".temp", fileName)        

    # set variable
    startTime = None
    stopTime = None

    # remove old outputs
    delFiles(fileName+"_.*")

    # convert transfer file to binary file 
    try:
       spacitCmd = open(spacitIn, "w")
       spacitCmd.write("T\n%s\n%s\nQ\nN"%(fileName, spacitOut))
       spacitCmd.close()
       cmd = "spacit < %s > %s"%(spacitIn, spacitlog)
       os.system(cmd)
    except Exception as e:
      print (e.message, e.args)

    # if successful, dump summary to the log file
    if os.path.exists(spacitOut):
       cmd = "ckbrief -utc %s %s %s>> %s"%(spacitOut, lsk, sclk, spacitlog)
       os.system(cmd)
    
    # read log file and extract info
    if os.path.exists(spacitOut):
       log = open(spacitlog, "r")
       lines = np.array([l.strip() for l in log])
       log.close()

       # let's check if dump was successful
       index = np.where(lines =='A traceback follows.  The name of the highest level module is first.')[0]
       if index:
          fn.raiseErr("\n".join(lines))
       
       # let's read further if successful
       index = np.where(lines=='Interval Begin UTC       Interval End UTC         AV')[0][0]+2
       info = lines[index].split()
       startTime = julian.cal_dt(info[0]+"T"+info[1])
       stopTime =  julian.cal_dt(info[2]+"T"+info[3])

    delFiles(fileName+"_.*")
    return startTime, stopTime, None, None, None 

#-------------------------------------------------------------------------------------------------------------------
def getEopRecords(fileName):
    # read EOP file and get records infor
    f = open(fileName, "r")
    s = []
    headerSize=None
    dataRecords=None
    recordSize=None
    headerInfo=None
    for t in f:
        if t.strip(): s += [len(t+"\n")]
        if t.startswith(" EOP="):
           headerInfo = [sum(s), len(s)]
           s = []
    f.close()

    if headerInfo is not None:
       headerSize = headerInfo[0]
       dataRecords = len(s)
       recordSize = s[0]

    return headerSize, dataRecords, recordSize

#----------------------------------------------------------------------------------------------------------------
def eopInfo (fileName):
    
    # read eop file
    f = open(fileName, 'r')
    getAll = np.array([d.strip() for d in f])
    f.close()

    # get index of the first data point
    index = np.where(getAll=='EOP=')[0][0]+1
    start = getAll[index].split(",")[0]
    start = julian.from_jd(float(start), fmt='mjd') 

    # get end time
    stop = getAll[-1].split(",")[0]
    stop = julian.from_jd(float(stop), fmt='mjd') 

    # get header info
    headerSize, dataRecords, recordSize = getEopRecords (fileName)

    return start, stop, headerSize, dataRecords, recordSize

#----------------------------------------------------------------------------------------------------------------
def mediaTimeFmt(time):
    
    # get output time formate
    time = time.split(',')
    time = time[0].replace("/","-")+"T"+time[-1]
    if int(time[:2]) > 60: time="19"+time
    else: time="20"+time

    return time

#----------------------------------------------------------------------------------------------------------------
def mediaInfo (fileName):
    
    # read ion file
    f = open(fileName, 'r')
    getAll = np.array([d.strip() for d in f if 'FROM' in d])
    f.close()

    # get start time
    start = getAll[0].split('TO')[0].split('FROM')[-1][1:-1] 
    start = mediaTimeFmt(start)

    # get stop time
    stop = getAll[-1].split('TO')[-1].split('DSN')[0][1:-1]
    stop = mediaTimeFmt(stop)

    return start, stop, None, None, None
    
#----------------------------------------------------------------------------------------------------------------
def getDataInfo (dType, fileName):
    
    if dType == 'spk':
       start, stop, headerSize, dataRecords, recordSize = spkInfo(fileName)
    
    if dType == 'ckf':
       start, stop, headerSize, dataRecords, recordSize = ckfInfo(fileName)

    if dType == 'eop':
       start, stop, headerSize, dataRecords, recordSize = eopInfo(fileName)

    if dType == 'ion' or dType == 'tro':   
       start, stop, headerSize, dataRecords, recordSize = mediaInfo(fileName)

    if start is None or stop is None:
       msg = "\nUnable to extract the file info!\nFile Name: %s\n"
       fn.raiseErr (msg%(fileName))
    
    return start, stop, headerSize, dataRecords, recordSize

#----------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def getTemplet(dType, velTemp, inputFile, nameNoExe, creTime, dsn, size, verbose):
    
    # setup update bar, if requested
    if verbose:
       prog = fn.progressBar()                                
       prog.start(1, "Reading and unpacking %s file"%(dType))

    # start and stop time
    start, stop, headerSize, records, recordSize = getDataInfo(dType, inputFile)
    if verbose: prog.stop()


    # updating velocity template
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Updating velocity template")
    
    md5 = fn.calculateMD5(inputFile)
    tempOut=fillFileArea (inputFile, velTemp, nameNoExe, start, stop, creTime, 
                          int(size), md5, headerSize, records, recordSize)
    
    prog.stop() 
    return tempOut
    
#-------------------------------------------------------------------------------------------------------------------
