# -*- coding: utf-8 -*-
""" 
The script reads the TNF (TRK-2-34) radio science data file and fills out the velocity template. 

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
import numpy as np
import datetime as dt
import functions as fn
from time import localtime, strftime 
from astropy.time import Time
from dsnContext import dsnContext

#-------------------------------------------------------------------------------------------------------------------
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
def getFormat(setup, byteLength):
    # define format of the data chunk
    dataFormat = byteFmt.byteOrder[:]

    fieldNames = []
    for f in setup:
        dataFormat += f[1]
        fieldNames.append( f[0] )
    
    format = struct.Struct( dataFormat )

    if format.size != byteLength:
       print(fn.colorTxt("\nError while formating\n", 'lightred'))
       sys.exit()

    return (fieldNames, format)

#-------------------------------------------------------------------------------------------------------------------
def getCHDOType (ddid):
    CHDOType = {
       'C123' : 'C132',
       'C124' : 'C133',
       'C125' : 'C134',
       'C126' : 'C135',
       'C127' : 'C136',
       }
    return CHDOType[ddid]  

#-------------------------------------------------------------------------------------------------------------------
def getStat(fileName):
    # get outfile file's statistics
    byte = readFile(fileName)
    size = len(byte)
    obj, end = getHeader(byte, 0)
    ddid    = obj["data_description_id"]
    rcrdLen = obj["sfdu_length"]+20
    rcrdSize = round(size/rcrdLen)

    return getCHDOType (ddid.decode("utf-8")), rcrdLen, rcrdSize

#-------------------------------------------------------------------------------------------------------------------
def dataType(i):
    # data type format ids
    formatIds = {
         0  : "Uplink Carrier Phase",    
         1  : "Downlink Carrier Phase",    
         2  : "Uplink Sequential Ranging Phase",
         3  : "Downlink Sequential Ranging Phase",    
         4  : "Uplink PN Ranging Phase",    
         5  : "Downlink PN Ranging Phase",    
         6  : "Doppler",    
         7  : "Sequential Ranging",    
         8  : "Angles",    
         9  : "Ramps",    
         10 : "VLBI",    
         11 : "DRVID",    
         12 : "Smoothed Noise",    
         13 : "Allan Deviation",    
         14 : "PN Ranging",    
         15 : "Tone Ranging",    
         16 : "Carrier Observable",    
         17 : "Total Phase",    
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
def getHeader (data, start):
    byteLength = 20
    end = start + byteLength

    # Tracking SFDU Label (Table 3-1 in TRK2-34 document).
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

    return obj, end 
    
#-------------------------------------------------------------------------------------------------------------------
def getAggregate (data, start):
    byteLength = 4 
    end = start + byteLength

    # Aggregation CHDO Label (Table 3-2 in TRK2-34 document).
    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),                   
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )

    return obj, end 

#-------------------------------------------------------------------------------------------------------------------
def getPrimary (data, start):
    byteLength = 8 
    end = start + byteLength

    # Primary CHDO (Table 3-2 in TRK2-34 document). 
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

    return obj, end 

#-------------------------------------------------------------------------------------------------------------------
def getSecSetUP(chdo_type):
    
    byteLength = 0
    setup = [()]
    
    # Secondary CHDO 132 (Table 3-5 in TRK2-34 document).
    if chdo_type == 132:
       
       byteLength = 66
       setup = [
       ( "orig_id", byteFmt.uint1 ),
       ( "last_modifier_id", byteFmt.uint1 ),
       ( "reserve1", "1s" ),
       ( "scft_id", byteFmt.uint1 ),
       ( "upl_rec_seq_num", byteFmt.uint4 ),
       ( "rec_seq_num", byteFmt.uint4 ),
       ( "year", byteFmt.uint2 ),
       ( "doy", byteFmt.uint2 ),
       ( "sec", byteFmt.float8 ),
       ( "rct_day", byteFmt.uint2 ),
       ( "rct_msec", byteFmt.uint4 ),
       ( "ul_dss_id", byteFmt.uint1 ),
       ( "ul_band", byteFmt.uint1 ),
       ( "ul_assembly_num", byteFmt.uint1 ),
       ( "transmit_num", byteFmt.uint1 ),
       ( "transmit_stat", byteFmt.uint1 ),
       ( "transmit_mode", byteFmt.uint1 ),
       ( "cmd_modul_stat", byteFmt.uint1 ),
       ( "rng_modul_stat", byteFmt.uint1 ),
       ( "fts_vld_flag", byteFmt.uint1 ),
       ( "reserve1a", "1s" ),
       ( "transmit_time_tag_delay", byteFmt.float8 ),
       ( "ul_zheight_corr", byteFmt.float4 ),
       ( "mod_day", byteFmt.uint2 ),
       ( "mod_msec", byteFmt.uint4 ),
       ( "version_num", byteFmt.uint1 ),
       ( "sub_version_num", byteFmt.uint1 ),
       ( "sub_sub_version_num", byteFmt.uint1 ),
       ( "reserve1b", "1s" ),
       ( "reserve4", "4s" ),
          ]
      
    # Secondary CHDO 133 (Table 3-6 in TRK2-34 document).   
    elif chdo_type == 133:
       byteLength = 110 
       setup = [
       ( "orig_id", byteFmt.uint1 ),
       ( "last_modifier_id", byteFmt.uint1 ),
       ( "reserve1", "1s" ),
       ( "scft_id", byteFmt.uint1 ),
       ( "dtt_rec_seq_num", byteFmt.uint4 ),
       ( "rec_seq_num", byteFmt.uint4 ),
       ( "year", byteFmt.uint2 ),
       ( "doy", byteFmt.uint2 ),
       ( "sec", byteFmt.float8 ),
       ( "rct_day", byteFmt.uint2 ),
       ( "rct_msec", byteFmt.uint4 ),
       ( "dl_dss_id", byteFmt.uint1 ),
       ( "dl_band", byteFmt.uint1 ),
       ( "dl_chan_num", byteFmt.uint1 ),
       ( "prdx_mode", byteFmt.uint1 ),
       ( "ul_prdx_stn", byteFmt.uint1 ),
       ( "ul_band_dl", byteFmt.uint1 ),
       ( "array_delay", byteFmt.float8 ),
       ( "fts_vld_flag", byteFmt.uint1 ),
       ( "carr_lock_stat", byteFmt.uint1 ),
       ( "array_flag", byteFmt.uint1 ),
       ( "polarization", byteFmt.uint1 ),
       ( "diplxr_stat", byteFmt.uint1 ),
       ( "lna_num", byteFmt.uint1 ),
       ( "rf_if_chan_num", byteFmt.uint1 ),
       ( "if_num", byteFmt.uint1 ),
       ( "rcv_time_tag_delay", byteFmt.float8 ),
       ( "dl_zheight_corr", byteFmt.float4 ),
       ( "vld_ul_stn", byteFmt.uint1 ),
       ( "vld_dop_mode", byteFmt.uint1 ),
       ( "vld_scft_coh", byteFmt.uint1 ),
       ( "scft_transpd_lock", byteFmt.uint1 ),
       ( "scft_transpd_num", byteFmt.uint1 ),
       ( "reserve1a", "1s" ),
       ( "scft_osc_freq", byteFmt.float8 ),
       ( "scft_transpd_delay", byteFmt.float8 ),
       ( "scft_transpd_turn_num", byteFmt.uint4 ),
       ( "scft_transpd_turn_den", byteFmt.uint4 ),
       ( "scft_twnc_stat", byteFmt.uint1 ),
       ( "scft_osc_type", byteFmt.uint1 ),
       ( "mod_day", byteFmt.uint2 ),
       ( "mod_msec", byteFmt.uint4 ),
       ( "version_num", byteFmt.uint1 ),
       ( "sub_version_num", byteFmt.uint1 ),
       ( "sub_sub_version_num", byteFmt.uint1 ),
       ( "lna_corr_value", byteFmt.uint1 ),
       ( "reserve4", "4s" ),
          ]
    # Secondary CHDO 134 (Table 3-4 in TRK2-34 document).           
    elif chdo_type == 134:
       byteLength = 124
       setup = [
       ( "orig_id", byteFmt.uint1 ),
       ( "last_modifier_id", byteFmt.uint1 ),
       ( "reserve1", "1s" ),
       ( "scft_id", byteFmt.uint1 ),
       ( "rec_seq_num", byteFmt.uint4 ),
       ( "year", byteFmt.uint2 ),
       ( "doy", byteFmt.uint2 ),
       ( "sec", byteFmt.float8 ),
       ( "rct_day", byteFmt.uint2 ),
       ( "rct_msec", byteFmt.uint4 ),
       ( "stn_stream_src", byteFmt.uint1 ),
       ( "ul_band", byteFmt.uint1 ),
       ( "ul_assembly_num", byteFmt.uint1 ),
       ( "transmit_num", byteFmt.uint1 ),
       ( "transmit_stat", byteFmt.uint1 ),
       ( "transmit_mode", byteFmt.uint1 ),
       ( "cmd_modul_stat", byteFmt.uint1 ),
       ( "rng_modul_stat", byteFmt.uint1 ),
       ( "transmit_time_tag_delay", byteFmt.float8 ),
       ( "ul_zheight_corr", byteFmt.float4 ),
       ( "dl_dss_id", byteFmt.uint1 ),
       ( "reserve1a", "1s" ),
       ( "dl_chan_num", byteFmt.uint1 ),
       ( "prdx_mode", byteFmt.uint1 ),
       ( "ul_prdx_stn", byteFmt.uint1 ),
       ( "ul_band_dl", byteFmt.uint1 ),
       ( "array_delay", byteFmt.float8 ),
       ( "fts_vld_flag", byteFmt.uint1 ),
       ( "carr_lock_stat", byteFmt.uint1 ),
       ( "array_flag", byteFmt.uint1 ),
       ( "lna_num", byteFmt.uint1 ),
       ( "rcv_time_tag_delay", byteFmt.float8 ),
       ( "dl_zheight_corr", byteFmt.float4 ),
       ( "vld_ul_stn", byteFmt.uint1 ),
       ( "vld_dop_mode", byteFmt.uint1 ),
       ( "vld_scft_coh", byteFmt.uint1 ),
       ( "vld_dl_band", byteFmt.uint1 ),
       ( "scft_transpd_lock", byteFmt.uint1 ),
       ( "scft_transpd_num", byteFmt.uint1 ),
       ( "reserve2", "2s" ),
       ( "scft_osc_freq", byteFmt.float8 ),
       ( "scft_transpd_delay", byteFmt.float8 ),
       ( "scft_transpd_turn_num", byteFmt.uint4 ),
       ( "scft_transpd_turn_den", byteFmt.uint4 ),
       ( "scft_twnc_stat", byteFmt.uint1 ),
       ( "scft_osc_type", byteFmt.uint1 ),
       ( "mod_day", byteFmt.uint2 ),
       ( "mod_msec", byteFmt.uint4 ),
       ( "cnt_time", byteFmt.float4 ),
       ( "version_num", byteFmt.uint1 ),
       ( "sub_version_num", byteFmt.uint1 ),
       ( "sub_sub_version_num", byteFmt.uint1 ),
       ( "lna_corr_value", byteFmt.uint1 ),
       ]
    
    # Secondary CHDO 135 (Table 3-7 in TRK2-34 document). 
    elif chdo_type == 135:
       byteLength = 88
       setup = [
       ( "orig_id", byteFmt.uint1 ),
       ( "last_modifier_id", byteFmt.uint1 ),
       ( "reserve1a", "1s" ),
       ( "scft_id", byteFmt.uint1 ),
       ( "rec_seq_num", byteFmt.uint4 ),
       ( "year", byteFmt.uint2 ),
       ( "doy", byteFmt.uint2 ),
       ( "sec", byteFmt.float8 ),
       ( "rct_day", byteFmt.uint2 ),
       ( "rct_msec", byteFmt.uint4 ),
       ( "ul_dss_id", byteFmt.uint1 ),
       ( "dl_dss_id", byteFmt.uint1 ),
       ( "dl_dss_id_2", byteFmt.uint1 ),
       ( "dl_band", byteFmt.uint1 ),
       ( "prdx_mode", byteFmt.uint1 ),
       ( "ul_band", byteFmt.uint1 ),
       ( "rec_type", byteFmt.uint1 ),
       ( "source_type", byteFmt.uint1 ),
       ( "fts_vld_flag", byteFmt.uint1 ),
       ( "reserve1b", "1s" ),
       ( "array_flag", byteFmt.uint1 ),
       ( "array_flag_2", byteFmt.uint1 ),
       ( "array_delay", byteFmt.float8 ),
       ( "array_delay_2", byteFmt.float8 ),
       ( "rcv_time_tag_del", byteFmt.float8 ),
       ( "rcv_time_tag_delay_2", byteFmt.float8 ),
       ( "mod_day", byteFmt.uint2 ),
       ( "mod_msec", byteFmt.uint4 ),
       ( "version_num", byteFmt.uint1 ),
       ( "sub_version_num", byteFmt.uint1 ),
       ( "sub_sub_version_num", byteFmt.uint1 ),
       ( "reserve1c", "1s" ),
       ( "reserve8", "8s" ),
       ]

    # Secondary CHDO 136 (Table 3-8 in TRK2-34 document).        
    elif chdo_type == 136:
      byteLength =98 
      setup = [
      ( "orig_id", byteFmt.uint1 ),
      ( "last_modifier_id", byteFmt.uint1 ),
      ( "reserve1", "1s" ),
      ( "scft_id", byteFmt.uint1 ),
      ( "rec_seq_num", byteFmt.uint4 ),
      ( "year", byteFmt.uint2 ),
      ( "doy", byteFmt.uint2 ),
      ( "sec", byteFmt.float8 ),
      ( "rct_day   ", byteFmt.uint2 ),
      ( "rct_msec", byteFmt.uint4 ),
      ( "dl_dss_id", byteFmt.uint1 ),
      ( "dl_band", byteFmt.uint1 ),
      ( "dl_chan_num", byteFmt.uint1 ),
      ( "prdx_mode", byteFmt.uint1 ),
      ( "ul_prdx_stn", byteFmt.uint1 ),
      ( "ul_band_dl", byteFmt.uint1 ),
      ( "rcv_time_tag_delay", byteFmt.float8 ),
      ( "array_delay", byteFmt.float8 ),
      ( "fts_vld_flag", byteFmt.uint1 ),
      ( "carr_lock_stat", byteFmt.uint1 ),
      ( "array_flag", byteFmt.uint1 ),
      ( "lna_num", byteFmt.uint1 ),
      ( "vld_ul_stn", byteFmt.uint1 ),
      ( "vld_dop_mode", byteFmt.uint1 ),
      ( "vld_scft_coh", byteFmt.uint1 ),
      ( "scft_transpd_lock", byteFmt.uint1 ),
      ( "scft_transpd_num", byteFmt.uint1 ),
      ( "reserve1a", "1s" ),
      ( "scft_osc_freq", byteFmt.float8 ),
      ( "scft_transpd_delay", byteFmt.float8 ),
      ( "scft_transpd_turn_num", byteFmt.uint4 ),
      ( "scft_transpd_turn_den", byteFmt.uint4 ),
      ( "scft_twnc_stat", byteFmt.uint1 ),
      ( "scft_osc_type", byteFmt.uint1 ),
      ( "mod_day", byteFmt.uint2 ),
      ( "mod_msec", byteFmt.uint4 ),
      ( "version_num", byteFmt.uint1 ),
      ( "sub_version_num", byteFmt.uint1 ),
      ( "sub_sub_version_num", byteFmt.uint1 ),
      ( "reserve1b", "1s" ),
      ( "reserve4", "4s" ),
      ]

    return byteLength, setup   

#-------------------------------------------------------------------------------------------------------------------
def getSecondary (data, start ):
    byteLength = 4 
    end = start + byteLength
    success = False

    # # Secondary CHDO (Table 3-4 to 3-8 in TRK2-34 document). 
    setup = [
       ( "chdo_type", byteFmt.uint2 ),
       ( "chdo_length", byteFmt.uint2 ),
       ]
    
    (_fields, _format) = getFormat(setup, byteLength)
    obj = unpack( _fields, _format, data, start )


    byteLength, setup = getSecSetUP (obj["chdo_type"])
    start = end
    end = start + byteLength

    if byteLength == 0 or len(data)-start < byteLength: 
       return obj, end, success
    else: success = True   

    (_fields, _format) = getFormat(setup, byteLength)
    obj.update(unpack( _fields, _format, data, start ))
     
    return obj, end, success 

#-------------------------------------------------------------------------------------------------------------------
def readFile(data):
     # read input file
     with open(data, "rb") as f: byte = f.read()
     f.close()
     return byte

#-------------------------------------------------------------------------------------------------------------------
def findStartByte(data, start, inputFile):
    
    # find start byte
    data = data[start:]
    dataLen = len( data )
    while start < dataLen:
          i = data.find(b"N")
          sfdu_class_id = struct.unpack('s', data[i+5:i+6])[0]
          if sfdu_class_id == b"I":
             startByte = i
             return start+startByte
          elif sfdu_class_id == b"K":
             i = data.find(b"CCSD$$MARKER$T-2-34$", i)
             startByte = i + 20 + 20
             return start+startByte
          else:
             if start !=0:
                print ("Warning: The file '%s' doesn't contain the appropriate header everywhere" %(inputFile))     
                return dataLen
             else:                
                print("Error: The file '%s' doesn't contain the appropriate header" %(inputFile))
                sys.exit()

#----------------------------------------------------------------------------------------------------------------
def catFiles (inp, out):
    try:
      if os.path.exists(out): os.remove(out)
      error=runCommands("cat %s  > %s" %(inp, out))
      if error:
         print (error)
         raise Exception(error)
    except Exception as err:
      print (fn.colorTxt("\nError while concatenating the output files.\nConcatenating Files are:\n%s\n\n"%(inp), 'lightred'))
      sys.exit()


#----------------------------------------------------------------------------------------------------------------
def docLid (tt):
    tt = Time(tt)
    if tt < Time("2008-02-29"): return "dsn.trk-2-34:2002-12-15"
    if tt >= Time("2008-02-29") and tt < Time("2013-11-07"): return "dsn.trk-2-34:2008-02-29"
    if tt >= Time("2013-11-07") and tt < Time("2015-10-27"): return "dsn.trk-2-34:2013-11-07"
    if tt >= Time("2015-10-27") and tt < Time("2017-05-03"): return "dsn.trk-2-34:2015-10-27"
    if tt >= Time("2017-05-03"): return "dsn.trk-2-34:2017-05-03"

#----------------------------------------------------------------------------------------------------------------
def fillFileArea (nameNoExe, creTime, start, stop, size, md5, dsn):
    
    tempFile = "_temp_%s.xml"%(nameNoExe)
    readFile  = open("/tmp/"+tempFile, "r")
    fileDate = readFile.read()
    readFile.close()
    year = int(start[:4])

    _dsn =  list(map(int, dsn))
    _dsn = str(_dsn)[1:-1] 

    newArea = fileDate.replace("FILEWITHNOEXE", nameNoExe)
    newArea = newArea.replace("CREATION_DAY", creTime[:10])
    newArea = newArea.replace("START_TIME", start)
    newArea = newArea.replace("STOP_TIME", stop)
    newArea = newArea.replace("FILE_NAME", nameNoExe+".dat")
    newArea = newArea.replace("CREATION_TIME", creTime[:10])
    newArea = newArea.replace("FILE_SIZE", str(size))
    newArea = newArea.replace("MD5", md5)
    newArea = newArea.replace("DSSUSED", str(_dsn))
    newArea = newArea.replace("DSNOBSSYS", dsnContext(dsn, year, nameNoExe))
    newArea = newArea.replace("DOCLID", "urn:nasa:pds:radiosci.documentation:"+docLid(start))

    outFile = open("/tmp/"+tempFile, "w")
    outFile.write(newArea)
    outFile.close()

    return tempFile

#----------------------------------------------------------------------------------------------------------------
def fillTable(tempOut, velTempDir, nameNoExe, dataType, table, offSet, rcrdLen, rcrdSize):

    fileName = "%s/table_tnf_g1g2g3.xml"%(velTempDir)
    outFile  = '%s/_temp_'%(tempOut)+nameNoExe+'_TableFilled'+str(table)+'.xml'
    readTab = open(fileName, "r")
    outTab  = open(outFile, "w")
    tabData = readTab.read()
    readTab.close()

    newTab = tabData.replace("DATA_TYPE", "%-33s"%dataType)
    newTab = newTab.replace("TABLENO", str(table))
    newTab = newTab.replace("OFFSET", str(offSet))
    newTab = newTab.replace("RECORDS", str(rcrdSize))
    newTab = newTab.replace("RECORDLENGTH", str(rcrdLen))

    outTab.write(newTab)
    outTab.close()

    return outFile 
#----------------------------------------------------------------------------------------------------------------
def getXMLTemp (tempOut, velTemp, msg, c, ind, nameNoExe, dataType='NAN', table=0, offSet=0, rcrdLen=0, rcrdSize=0, year=0):
    tempFile = "/tmp/_temp_%s.xml"%(nameNoExe)

    velTempDir = os.path.dirname(velTemp)
    if ind is not None:
       files2cat = fillTable(tempOut, velTempDir, nameNoExe, dataType, table, offSet, rcrdLen, rcrdSize)
       msg = msg + '%s %s/tnf_g4-%s.xml %s/tnf_g5-DT%s.xml %s/tableEnd.xml ' %(files2cat, velTempDir, c, velTempDir, ind, velTempDir)
       return msg
    else:
       msg = velTemp + msg + " %s/endXML.xml"%(velTempDir)
       catFiles (msg, tempFile)

#-------------------------------------------------------------------------------------------------------------------
def getTimeFormat (time):
    time = str(time).split()
    time = time[0]+"T"+time[1]

    return time

#-------------------------------------------------------------------------------------------------------------------
def writeVelocityTemp(tempOut, secPast, inputFile, outFile, velTemp, nameNoExe, creTime, size, dsn, md5):
    
    # get output file name
    outDir = os.path.dirname(outFile)
    xmlFile = outDir+"/"+nameNoExe+".xml"

    # get start and stop time of the output file
    startEpoch = getTimeFormat(min(secPast))
    stopEpoch  = getTimeFormat(max(secPast))
    del secPast

    startByte = 0
    # get statistics for each data type
    msg = ' '
    for i in range(18):
        fileName = "%s/stnfPY%s.out"%(tempOut,"%02d"%i)
        if os.path.exists(fileName):
           ddid, rcrdLen, rcrdSize = getStat(fileName)
           os.remove(fileName)
           msg = (getXMLTemp (tempOut, velTemp, msg, str(ddid), "%02d"%i, nameNoExe, dataType=dataType(i), 
                  table="%02d"%i, offSet=int(startByte), rcrdLen=int(rcrdLen), rcrdSize=int(rcrdSize)))
           
           startByte = startByte + rcrdLen*rcrdSize

    getXMLTemp (tempOut, velTemp, msg, xmlFile, None, nameNoExe, year=int(startEpoch[:4]))
    tempOut=fillFileArea (nameNoExe, creTime, startEpoch, stopEpoch, startByte, md5, dsn)

    for tempFile in glob.glob('/tmp/_temp_'+nameNoExe+'_TableFilled*'):
        os.remove(tempFile)
    
    return tempOut
    
#-------------------------------------------------------------------------------------------------------------------
# main program
#-------------------------------------------------------------------------------------------------------------------
def getTemplet(velTemp, inputFile, outFile, nameNoExe, creTime, size, verbose):
    
    # setup update bar, if requested
    if verbose:
       prog = fn.progressBar()                                  
       prog.start(1, "Reading and shorting TNF file")

    # create temp output directory
    _tempOut = "/tmp/_temp_"+nameNoExe
    fn.rmDir(_tempOut, createDir=True)
    fn.rmFiles(outFile)

    # read input file
    byte=readFile(inputFile)
    size = len(byte)
    
    # initialize arrays and variable
    outfile = [0]*18 
    nsfdu = [0]*18  
    secPast = []
    orec = 0
    dss = []
    
    # get start byte
    start = findStartByte(byte, 0, inputFile)
    
    # start loop
    while start < size:
          # read tracking SFDU label 
          objH, end = getHeader(byte, start)
          if objH["control_auth_id"].decode("utf-8") != 'NJPL':
             start = findStartByte(byte, start, inputFile)
             objH, end = getHeader(byte, start)          
    
          # read aggregation CHDO label
          objA, end = getAggregate (byte, end)
    
          # read primary CHDO data
          objP, end = getPrimary(byte, end)
    
          # read secondary CHDO data
          objS, end, success = getSecondary (byte, end)
    
          # update chunkSize
          chunkSize = objH['sfdu_length'] + 20
    
          # get measurement time 
          if success: secPast.append(dt.datetime(objS["year"], 1, 1) + dt.timedelta(days=objS["doy"] - 1, seconds=objS["sec"]))

          # output array
          orec = orec + 1
          oarray = byte[start:start+chunkSize]
    
          # This section screens for rogue records, prevents their being written 
          if chunkSize <= len(byte) and size-start >= chunkSize:
             if objH["data_description_id"].decode("utf-8") != "C123" and \
                objH["data_description_id"].decode("utf-8") != "C124" and \
                objH["data_description_id"].decode("utf-8") != "C125" and \
                objH["data_description_id"].decode("utf-8") != "C126" and \
                objH["data_description_id"].decode("utf-8") != "C127" or  \
                objP["chdo_type"] != 2 or  objP['chdo_length'] != 4 or  \
                objP["mjr_data_class"] != 6 or  objP["mnr_data_class"] != 14 :

                if not success:
                   print (objH["data_description_id"].decode("utf-8"))
             
             else:
                # write output buffer to output array
                outfilePath = "%s/stnfPY%s.out"%(_tempOut, "%02d"%objP["format_code"])
                outfile[objP["format_code"]] = open(outfilePath, "ab")
                outfile[objP["format_code"]].write(oarray)
                
                if 'dl_dss_id' in objS.keys(): dss.append(str(objS['dl_dss_id']))
                if 'ul_dss_id' in objS.keys(): dss.append(str(objS['ul_dss_id']))
    
                nsfdu[objP["format_code"]] = nsfdu[objP["format_code"]] + 1
             
          # update start byte
          start = start + chunkSize 
          if verbose: prog.setStep(start*1.0/size)
    
    # get station
    dsn = np.unique(dss) 

    # empty memory
    del byte
    closeOutFile(outfile)
    
    # concatenating the output files
    catFiles (_tempOut+'/stnfPY*.out', outFile)
    md5 = fn.calculateMD5(outFile)
    if verbose: prog.stop()
    
    # updating velocity template
    if verbose:
       prog = fn.progressBar()
       prog.start(1, "Updating velocity template")
    tempOut=writeVelocityTemp(_tempOut, secPast, inputFile, outFile, velTemp, nameNoExe, creTime, size, dsn, md5)

    # remove tempOut file
    fn.rmDir(_tempOut)
    if verbose: prog.stop()

    return tempOut
#-------------------------------------------------------------------------------------------------------------------
