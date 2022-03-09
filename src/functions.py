 # -*- coding: utf-8 -*-
"""

Miscellaneous functions

AUTHOR  
        
   Ashok Verma (ashokverma@ucla.edu)
   Department of Earth, Planetary, and Space Sciences
   University of California, Los Angeles
"""

#-------------------------------------------------------------------------------------------------------------------
import os, glob, sys, traceback
import datetime as dt
import shutil as sh
from inspect import currentframe, getframeinfo
import signal, importlib
cf = currentframe()
sn = getframeinfo(cf).filename.split("/")[-1]

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
def getLineNumber():
    '''
    Obtain line number of the code.

    '''

    cf = currentframe()
    return cf.f_back.f_lineno

#-------------------------------------------------------------------------------------------------------------------
def raiseErr (msg):

    '''

    Raise exception when an error occurs.

    INPUTS:

    --msg  Error message 

    '''

    rmPYC()
    msg = "\nError:\n" + str(msg)
    print((colorTxt(msg, "red")))
    os.killpg(0, signal.SIGKILL)
    raise Exception("") 

#-------------------------------------------------------------------------------------------------------------------
def raiseWarning (msg):
    
    '''
    Raise warning.

    INPUTS:

    --msg  Warning message

    '''
    print((bcolors.fg.pink+'\n' + msg + '\n'+bcolors.reset))

#-------------------------------------------------------------------------------------------------------------------
def colorTxt (txt, clr):                                   

    '''
    Change the color of the text.

    INPUTS:

    --txt  Input text
    --clr  Color of the text

    '''
    mod=getattr(bcolors.fg, clr)
    return mod+txt+bcolors.reset

#-------------------------------------------------------------------------------------------------------------------
def rmDir(path, isPathDir=True, createDir=False):

    '''

    Remove directories (parent and child) or file. 

    INPUTS:

    --path  Path of the directory that needs to be removed.
    --isPathDir  Flag that describes the type of the path.
    --createDir  Flag to create the empty copy of the deleted directory. 

    '''
    if os.path.exists(path):
         if isPathDir: 
            sh.rmtree(path)
         else:
            os.remove(path)
    
    if createDir:
       os.mkdir(path)
       os.chmod(path, 0o774)

#-------------------------------------------------------------------------------------------------------------------
def calculateMD5(filename, block_size=2**20):
        """Returns MD% checksum for given file.
        """
        import hashlib

        md5 = hashlib.md5()
        try:
                file = open(filename, 'rb')
                while True:
                        data = file.read(block_size)
                        if not data:
                                break
                        md5.update(data)
        except IOError:
                colorTxt('\nFile \'' + filename + '\' not found!\n', 'lightred')
                if os.path.exists(filename): os.remove(filename)
                return None
        except:
                if os.path.exists(filename): os.remove(filename)
                return None
        return md5.hexdigest()

#-------------------------------------------------------------------------------------------------------------------
def replaceElement (_list, str2replace, str2replaceWith):
    
    '''

    Replace element of the list.

    INPUTS:

    --_list  Path of the list
    --str2replace  String that needs to be replaced.
    --str2replaceWith  String that will be replaced with.


    '''
    _bool = [str2replace==x for x in _list]
    return any(_bool), [str2replaceWith if x==str2replace else x for x in _list]

#-------------------------------------------------------------------------------------------------------------------
def checkRunStatus (runFlag, cmd, script="", lineNumber=""):
    
    '''

    Check run status.

    INPUTS:

    --runFlag  Error message of the run, if any
    --cmd  The command that is running. 
    --script  The name of the script that is running. 
    --lineNumber  Line number from where the error evoke.

    '''
    if "error" in str(runFlag).lower() and len(str(runFlag)):
       rmPYC()

       msg = "="*100+"\n\n"
       msg += "%s \n Error while '%s' \n Script Name: %s \n Line Number: %s\n\n" %(runFlag, cmd, script, lineNumber)
       msg += "="*100+"\n\n"
       raiseErr(msg)

#-------------------------------------------------------------------------------------------------------------------
def subDirPath (d):
    
    '''

    Get absolute path of all child directories.

    INPUTS:

    --d  The path to the parent directory.

    '''
    return list(filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)]))

#-------------------------------------------------------------------------------------------------------------------
def move (src, des):
 
    '''

    Rename or move the file.

    INPUTS:

    --src  Source file name.
    --des  Destination file name.

    '''

    sh.move(src, des)

#-------------------------------------------------------------------------------------------------------------------
def copy (src, des, isDesDir=False, exclude='NONE'):

    '''

    Copy files or directores.

    INPUTS:

    --src  Source file/directory name.
    --des  Destination file/directory name.
    --isDesDir  Flag indicating whether destination is a directory.
    --exclude  Files to exclude.
    
    '''

    if isDesDir:
       rmDir (des, createDir=True)
    for fileORdir in glob.glob(src):
        if exclude != fileORdir:
           cmd = "cp -r '%s' '%s'" %(fileORdir, des)
           cpFlag = runSubProcess (cmd)
           checkRunStatus(cpFlag, cmd, sn, getLineNumber())

#-------------------------------------------------------------------------------------------------------------------
def pathcheck(pathList, script="", lineNumeber=""):
    
    '''

    Check path of the input files.

    INPUTS:

    --pathList  The path of the input file or the list of files.
    --script  Name of the script that runs this function. 
    --lineNumeber  Line number from where the error evoke.

    '''

    if type(pathList) is list:
       for path in pathList:
           path = os.path.abspath(path)
           if not os.path.exists(path):
              checkRunStatus("Error: No such file or directory", "path check: '%s'" %(path), script, lineNumeber)
    else:
       path = os.path.abspath(pathList)
       if not os.path.exists(path):
          checkRunStatus("Error: No such file or directory", "path check: '%s'" %(pathList), script, lineNumeber)

#-------------------------------------------------------------------------------------------------------------------
def rmFiles(des, exclude=[]):
    
    '''

    Remove multiple files.

    INPUTS:

    --des  The path of the input file or directory or wildcard (``yourpath/*.pyc``).
    --exclude  List of files that needs to be excluded.

    '''

    getAll = glob.glob(des)

    if len(exclude) > 0:
       getAll = []
       for excl in exclude:
           files = []
           for trash in glob.glob(des):
               if excl not in trash:
                  files.append(trash)
           
           getAll.append (files)
       getAll = list(set.intersection(*list(map(set, getAll))))           

    for trash in getAll:
        if os.path.exists(trash):os.remove(trash)

#-------------------------------------------------------------------------------------------------------------------
def rmPYC():

    '''

    Remove Python interpreter files.

    '''
    for gin in glob.glob("./*.pyc"):
        os.remove(gin)

#-------------------------------------------------------------------------------------------------------------------
def time_now ():

    '''

    Get local time of the run.

    '''
    return dt.datetime.now()

#-------------------------------------------------------------------------------------------------------------------
def time_diff (time_start, time_end):
    
    '''

    Get the difference bwtween two times.

    INPUTS:

    --tStart  Start time.
    --tEnd  End time.
    
    '''
    sec = (time_end - time_start).total_seconds()

    hr = sec / (60 * 60)
    mn = (hr - int(hr)) * 60
    sec = (mn - int(mn)) * 60

    return str("Time taken: %s Hr: %s Min: %s Sec \n" % (int(hr), int(mn), int(sec)))

#-------------------------------------------------------------------------------------------------------------------
class progressBar:
      
      '''
      
      Track progress for the running operations.
      
      '''

      def __init__( self, maxLabelLen=36, blockSize=20 ):
         
         '''

         INPUTS:

         --maxLabelLen  Length of the longest job label.
         --blockSize  The size of the progress bar.         

         '''
         self.colLen = maxLabelLen
         self.step = None
         self.totalSteps = None
         self.lastOutputPct = None
         self.blockSize = blockSize
         self.st = time_now()
      
      def start( self, numStep=100, label=""):
         
         '''
         Start a job.

         This will print a job label and start the counter.

         INPUTS:

         --numStep  The number of steps in the job.
         --label  The label for the job. If this is longer than maxLabelLen, it will be truncated.

         '''

         self.step = 0
         self.label = "%40s" %label.replace("Postfit measurements","Filter Solution").title()
         self.totalSteps = numStep
         self.lastOutputPct = 0
         self.localStep = 0
      
         lbl = label
         if len( lbl ) > self.colLen:
            lbl = label[:self.colLen]
      
         clr = bcolors.fg.red
         status = ""
         text = "\r"+self.label+": [{0}] {1}% {2}".format( clr + "█"*self.localStep + \
                "-"*(self.blockSize-self.localStep) + bcolors.reset, 0, status)

         sys.stdout.write(text)
         sys.stdout.flush()
      
      def setLabel( self, label ):
          self.label = label

      def increment( self ):

         '''

         Increment the current step by 1.

         '''

         self.step += 1
         self.setStep( self.step )
      
      def setStep( self, step ):
         
         '''

         Set the current step to the input value.

         INPUTS:

         --step  The current step to the input value.

         '''
         self.step = step
      
         x = float( self.step ) / float( self.totalSteps )
         pct = float( self.step ) / float( self.totalSteps ) * 100
         delta = pct - self.lastOutputPct
      
         clr = bcolors.fg.red
         status = ""
         chrSize = 100/self.blockSize
         while delta >= chrSize:
            self.localStep += 1
            text = "\r"+self.label+": [{0}] {1}% {2}".format( clr + "█"*self.localStep + \
                   "-"*(self.blockSize-self.localStep) + bcolors.reset, int(pct), status)
            sys.stdout.write(text)
      
            self.lastOutputPct += chrSize 
            delta -= chrSize 
      
         sys.stdout.flush()
      
      def stop( self ):

         '''

         Stop the job.

         '''

         clr = bcolors.fg.green
         status = time_diff(self.st, time_now())
         text = "\r"+self.label+": [{0}] {1}% {2}\n".format( clr + "█"*self.blockSize +\
                bcolors.reset, 100, status)
         sys.stdout.write( text )
         sys.stdout.flush()
      
#-------------------------------------------------------------------------------------------------------------------
