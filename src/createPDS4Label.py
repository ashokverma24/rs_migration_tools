""" 
USAGE    

    python createPDS4Label.py [-h,--help] 
    python createPDS4Label.py -i inputFile -v velocityTemplate [options...]
    
DESCRIPTION

    This script generates PDS4 XML labels for radio science ancillary files. 
    Since the DATA_SET_ID, PRODUCT_ID, and TARGET are used in the migration process, 
    it is imperative that the original PDS3 ancillary data file and its label are placed in the same folder. 

INPUTS:
    1.  -i   (inputFile)    :-  Path to the input data file.                        
    2.  -v   (velTemp)      :-  Path to the reference velocity template.                        
    3. [-u]  (usrDir)           [Optional] The path to the user-defined options directory. The directory should contain:
                                           1.  options.txt:  a file containing all variables of the run.
                                           2.  targLid.txt:  a file containing the list of target LIDs.
                                           3.  targType.txt: a file containing the list of the target types.
                                Default directory is: ../usrInputs/ 
    4. [-o]  (outDir)       :- [Optional] Path to the output directory. Default is set to the home directory of this script.

OUTPUT:
    1. Generate PDS4 XML label. 
   
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
import ancillaryReader
import tdfReader
import rsrReader
import tnfReader
import odfReader
import hashlib
from time import localtime, strftime 
import datetime as dt
from itertools import (takewhile,repeat)  

#-------------------------------------------------------------------------------------------------------------------
def delFiles (wildCard):
    for f in glob.glob(wildCard):
        if os.path.exists(f): os.remove(f)

#-------------------------------------------------------------------------------------------------------------------
def del_oldOutFile():
    # delete old output files, if exists
    for i in range(19):
        outFile = "stnfpy%s.out"%("%02d"%i)
        if os.path.exists(outFile):
           os.remove(outFile)

#-------------------------------------------------------------------------------------------------------------------
def closeOutFile(outfile):
    # close output files
    for i in range(len(outfile)):
        if hasattr(outfile[i], 'close'):
           try:
             outfile[i].close()
           except AttributeError:
             raise 

#-------------------------------------------------------------------------------------------------------------------
def runCommands(commandArgs):
    # run shell commond within python
    proc = sp.Popen(commandArgs,shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()

    return err

#-------------------------------------------------------------------------------------------------------------------
def is_binary_file(fileName):
    is_binary = False
    try:
        with open(fileName, "r") as f:
             for l in f: break
    except UnicodeDecodeError:
       is_binary = True
       pass
    
    return is_binary

#-------------------------------------------------------------------------------------------------------------------
def getOutFileName (fileName, outDir):
    # define the name of output files 
    split = fileName.split("/")[-1]
    exe = split.split(".")[-1]

    if os.path.exists(fileName[:-(len(exe)+1)]+".lbl"): lblFile = fileName[:-(len(exe)+1)]+".lbl"
    elif os.path.exists(fileName[:-(len(exe)+1)]+".LBL"): lblFile = fileName[:-(len(exe)+1)]+".LBL"
    else:
       msg = fn.colorTxt("Searched for %s or %s supporting lbl file, but none of them presents."\
                      %(fileName[:-(len(exe)+1)]+".lbl", fileName[:-(len(exe)+1)]+".LBL"), 'red')
       fn.raiseErr (msg)

    nameNoExe = split[:-(len(exe)+1)].lower()

    # check if input file is binary or text
    is_binary = is_binary_file(fileName)
    if is_binary: inputFile = outDir+"/"+nameNoExe+".dat"
    else: inputFile = outDir+"/"+nameNoExe+".tab"

    return is_binary, inputFile, nameNoExe, lblFile        

#-------------------------------------------------------------------------------------------------------------------
def readLbl(dType, filename, lblFile, fileNoExe, usrDir):
    tmplbl = "/tmp/_temp_%s.lbl"%(fileNoExe)
    out = open(tmplbl, "w")
    lbl = open(lblFile, "r")
    dsn=None
    target=None
    for m in lbl:
        msg = m.strip().lower()
        if "description" in msg: break
        if "dsn_station" in msg:
           if "{" in msg: dsn = msg.split("=")[-1].strip()[1:-1].split(",")
           else: dsn = [msg.split("=")[-1].strip()]
        if "start_time" in msg: start = msg.split("=")[-1].strip()
        if "stop_time" in msg: stop = msg.split("=")[-1].strip()
        if "target_name" in msg: 
           target = msg.split("=")[-1].strip()
           if '''"''' == target[0]: target = target[1:-1]
           if """'""" == target[0]: target = target[1:-1]
      
        out.write(m)
    lbl.close()
    out.close()
    
    # load all target types
    allTargs = {}
    for l in open(usrDir+"/targType.txt").readlines():
        l = l.split("=")
        allTargs.update({l[0].strip().lower():l[1].strip()})

    # check if target is available in the PDS3 label
    if target is None:
       msg = "\n=====>  The TARGET name is missing in the PDS3 label. <=====\n" \
             "Label Name: %s \n"
       if os.path.exists(filename): os.remove(filename)
       fn.raiseErr (msg%(lblFile))
    else:
       # check if we got the valid target!
       try:
           targType = allTargs[target]
       except KeyError as e:
          msg = "\n=====>  Invalid TARGET ['%s'] in the PDS3 label. <=====\n" \
                "Label Name: %s \n\n" \
                "Please check or edit the PDS3 label. Alternatively you can add this target to the existing list at %s/targType.txt\n"
          if os.path.exists(filename): os.remove(filename)
          fn.raiseErr (msg%(target, lblFile, usrDir))

    if target == 'unk':
       msg="\nWarning:\nFileName: %s\nThe target name on the PDS3 label is 'UNK'.\nMake sure this is the correct target type.\n"
       print (fn.colorTxt(msg%(lblFile), 'purple'))

    if dType in ['rsr', 'tnf', 'odf']: return dsn, tmplbl
    if dsn is not None:
       dsn = list(map(str.strip, dsn))
       for _dsn in dsn:
          if not _dsn:
             msg = "\n=====>  The DSN station(s) is(are) missing in the PDS3 label. <=====\n" \
                   "Label Name: %s \n"
             if os.path.exists(filename): os.remove(filename)
             fn.raiseErr (msg%(lblFile))

    return dsn, tmplbl

#-------------------------------------------------------------------------------------------------------------------
def updateDate (date):

    ymd = date.split("T")[0]
    hms = date.split("T")[1]

    ymd = ymd.split("-")

    if len(ymd) == 2:
       y   = int(ymd[0])
       doy = int(ymd[1])
       ymd = (str(dt.datetime(y, 1, 1) + dt.timedelta(days=doy - 1))[:10])
       date = ymd+"T"+hms

    return date

#-------------------------------------------------------------------------------------------------------------------
def getTime (inputFile):
    
    with open(inputFile, "rb") as f:
         first = f.readline()
         for last in f: pass
    f.close()

    startTime = updateDate(first.split(b",")[-1].decode().strip())
    stopTime  = updateDate(last.split(b",")[-1].decode().strip())

    return startTime, stopTime

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
def runDocgen(tempOut, inputFile, outDir, nameNoExe, lblFile, usrDir):                                                                          

    cmd = "docgen -t /tmp %s "\
          "-o %s/%s.xml  "\
          "opt:%s/options.txt  "\
          "lbl:%s  "\
          "targ_type:%s/targType.txt  "\
          "targ_lid:%s/targLid.txt;  "\
          "xml-finish.150714.pl  "\
          "%s/%s.xml"

    try:
       error=runCommands(cmd %(tempOut, outDir, nameNoExe, usrDir, lblFile, usrDir, usrDir, outDir, nameNoExe))
       if error:
          os.remove(inputFile)
          print(fn.colorTxt("\nDOCGEN COMMAND: %s\n"%(cmd%(tempOut, outDir, nameNoExe, usrDir, lblFile, usrDir, usrDir, outDir, nameNoExe)), 'blue'))
          fn.raiseErr(fn.colorTxt(str(error), 'red'))
    except Exception as err:
          fn.raiseErr(fn.colorTxt(str(error), 'red'))

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
def carriageReturn (is_binary, inputFile, outFile): 

    if is_binary:
       shutil.copy(inputFile, outFile)
    else:   
       # read input and write Carriage-Return Line-Feed record delimiters output
       inpF = open(inputFile, "r")
       outF = open(outFile, "w")

       for txt in inpF:
           if txt.strip(): outF.write(txt.rstrip()+"\r\n")

       inpF.close()
       outF.close()

    return outFile

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
def getDataInfo (dType, velTemp, inputFile, nameNoExe, creTime, dsn, size, progBar):
    
    # read all options
    global options
    rawFile  = options.inputFile
    tempOut = None

    # Raw data type
    if dType == 'tdf':
       tempOut = tdfReader.getTemplet(velTemp, inputFile, nameNoExe, creTime, dsn, size, progBar) 

    if dType == 'odf':
       tempOut = odfReader.getTemplet(velTemp, inputFile, nameNoExe, creTime, size, progBar) 

    if dType == 'tnf':
       tempOut = tnfReader.getTemplet(velTemp, rawFile, inputFile, nameNoExe, creTime, size, progBar) 
    
    if dType == 'rsr':
       tempOut = rsrReader.getTemplet(velTemp, rawFile, inputFile, nameNoExe, creTime, size, progBar) 
    
    # Ancillary data type
    if dType in ['spk', 'ckf', 'ion', 'tro', 'eop']:
       tempOut = ancillaryReader.getTemplet(dType, velTemp, inputFile, nameNoExe, creTime, dsn, size, progBar) 
    
    if tempOut is None:
       msg = "An unrecognized data type, %s, was encountered."
       fn.raiseErr(msg)

    return tempOut 

#----------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def main():
    
    # read all options
    global options
    inputFile  = options.inputFile
    velTemp    = options.velTemp
    outDir     = options.outDir
    usrDir     = options.usrDir
    progBar    = options.progBar

    # get local running time
    t1 = localtime()
    
    # get data type info
    dType = dataTypeInfo (velTemp)

    # define output file name based on input file
    is_binary, outFile, nameNoExe, lblFile = getOutFileName (inputFile, outDir)  

    # write Carriage-Return Line-Feed record delimiters output
    inputFile = carriageReturn(is_binary, inputFile, outFile)
    
    # find size of the input file
    size = os.path.getsize(inputFile)

    # delete old outputs, if exists
    del_oldOutFile()
    
    # read PDS3 label file 
    dsn, lblFile = readLbl(dType, inputFile, lblFile, nameNoExe, usrDir)

    # creation date
    creTime  = strftime("%Y-%m-%d", localtime())
    
    # start and stop tim
    tempOut = getDataInfo(dType, velTemp, inputFile, nameNoExe, creTime, dsn, size, progBar) 

    # write xml label
    if progBar: 
       prog = fn.progressBar()
       prog.start(1, "Writing PDS4 Label")
    runDocgen(tempOut, inputFile, outDir, nameNoExe, lblFile, usrDir)

    # remove tempOut file
    os.remove("/tmp/%s"%(tempOut))
    os.remove(lblFile)
    if progBar: prog.stop()
    
#-------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'])
        parser.add_option ('-i', '--inputFile',  action='store', default='', help='Path to the input data file.')
        parser.add_option ('-v', '--velTemp',  action='store', default='', help='Path to the reference velocity template.')
        parser.add_option ('-u', '--usrDir', action='store', default='./../usrInputs/', help='The path to the user-defined options directory.')
        parser.add_option ('-o', '--outDir', action='store', default='./output', help='Path to the output directory. Default is set to the home directory of this script.')
        parser.add_option ('-p', '--progBar', action='store', default=True, help='Show progress bar.')
        (options, args) = parser.parse_args()

        if options.inputFile:
           if not os.path.exists (options.inputFile): 
              print (fn.colorTxt("\nThe path specified for the file %s is not valid.\n" %options.inputFile, 'red'))
              sys.exit(0)
        else:
           print (fn.colorTxt("\nPlease '-i' option to provide path to the input file (use -h for help).\n", 'red'))
           sys.exit(0)

        if options.velTemp:  
           if not os.path.exists (options.velTemp): 
              print (fn.colorTxt("\nThe path specified for the file %s is not valid.\n" %options.velTemp, 'red'))
              sys.exit(0)
        else:
           print (fn.colorTxt("\nPlease use '-v' option to provide path to the reference velocity template (use -h for help).\n", 'red'))
           sys.exit(0)

        if not os.path.exists(options.usrDir):
           print (fn.colorTxt("\nThe user defined directory '%s' not found.\n"%(options.usrDir), 'red'))
           sys.exit(0)
        else:
           expectedFiles = ['options.txt', 'targLid.txt', 'targType.txt']
           for ef in expectedFiles:
               if not os.path.exists(options.usrDir+"/"+ef):
                  print (fn.colorTxt("\nUnder user defined directory: '%s' file was expected, but not found.\n"%(options.usrDir+"/"+ef), 'red'))
                  sys.exit(0)

        if not os.path.exists (options.outDir): 
              try:
                 os.mkdir(options.outDir)        
              except Exception as err:
                   print(fn.colorTxt("\n Error while creating output directory.\n", 'red'))
                   sys.exit(0)

        main()
        sys.exit()
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
        sys.modules[__name__].__dict__.clear()
    except SystemExit as e: # sys.exit()
        raise e
        sys.modules[__name__].__dict__.clear()
    except Exception as e:
        msg = "ERROR, UNEXPECTED EXCEPTION \nFile: %s" %(options.inputFile)
        print (fn.bcolors.fg.red+'\n'+msg+'\n')
        traceback.print_exc()
        print ('\n '+fn.bcolors.reset)
        sys.modules[__name__].__dict__.clear()

