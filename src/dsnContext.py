""" 
Context product for the NASA Deep Space Networks
"""
import numpy as np
#-------------------------------------------------------------------------------------------------------------------
def raiseErr (err, inputFile):
     err = "\n\n\nFile Name: %s\n"%(inputFile) + err
     fn.colorTxt('\n' + "Error:\n" + err + '\n', 'red')
     with open("ErrorLog.txt", 'a+') as f: f.write(err)

     raise KeyError()

#-------------------------------------------------------------------------------------------------------------------
def facility(dsn, inputFile):
    dsnFacility = {
       'dss33' :   'canberra.dss33_34m',            
       'dss34' :   'canberra.dss34_34m',
       'dss35' :   'canberra.dss35_34m',
       'dss36' :   'canberra.dss36_34m',
       'dss42' :   'canberra.dss42_26m',
       'dss42' :   'canberra.dss42_34m',
       'dss43' :   'canberra.dss43_64m',
       'dss43' :   'canberra.dss43_70m',
       'dss45' :   'canberra.dss45_34m',
       'dss46' :   'canberra.dss46_26m',
       'dss11' :   'goldstone.dss11_26m', 
       'dss12' :   'goldstone.dss12_26m',
       'dss12' :   'goldstone.dss12_34m',
       'dss13' :   'goldstone.dss13_26m',
       'dss13' :   'goldstone.dss13_34m',
       'dss14' :   'goldstone.dss14_64m',
       'dss14' :   'goldstone.dss14_70m',
       'dss15' :   'goldstone.dss15_34m',
       'dss16' :   'goldstone.dss16_26m',
       'dss23' :   'goldstone.dss23_34m',
       'dss24' :   'goldstone.dss24_34m',
       'dss25' :   'goldstone.dss25_34m',
       'dss26' :   'goldstone.dss26_34m',
       'dss27' :   'goldstone.dss27_34m',
       'dss28' :   'goldstone.dss28_34m',
       'dss53' :   'madrid.dss53_34m', 
       'dss54' :   'madrid.dss54_34m',
       'dss55' :   'madrid.dss55_34m',
       'dss56' :   'madrid.dss56_34m',
       'dss61' :   'madrid.dss61_26m',
       'dss61' :   'madrid.dss61_34m',
       'dss63' :   'madrid.dss63_64m',
       'dss63' :   'madrid.dss63_70m',
       'dss65P' :  'madrid.dss65_34m-post2005',
       'dss65' :   'madrid.dss65_34m',
       'dss66' :   'madrid.dss66_26m',
       'dss74' :   'new_norcia.dss74_35m',
       'dss84' :   'malarguee.dss84_35m',
       'dss47' :   'paul_wild.dss47_2008',
    }
    try:
       return dsnFacility[dsn]
    except KeyError:  
       msg = "\n====>The context product is not found for %s.<===="
       raiseErr (msg%(dsn), inputFile)

def commnCntx (dsn):
    msg = """
                  <Observing_System_Component>
                         <name>DSN Instrumentation</name>
                         <type>Instrument</type>
                         <Internal_Reference>
                               <lid_reference>urn:nasa:pds:context:instrument:dsn.rss</lid_reference>
                               <reference_type>is_instrument</reference_type>
                         </Internal_Reference>
                   </Observing_System_Component>

                   <Observing_System_Component>
                         <name>NASA Deep Space Network</name>
                         <type>Observatory</type>
                         <Internal_Reference>
                               <lid_reference>urn:nasa:pds:context:facility:observatory.dsn</lid_reference>
                               <reference_type>is_facility</reference_type>
                         </Internal_Reference>
                   </Observing_System_Component>
                   """
    if np.size and np.any(dsn == '47'):
        msg = msg + """ 
                   <Observing_System_Component>
                         <name>Paul Wild Observatory</name>
                         <type>Observatory</type>
                         <Internal_Reference>
                               <lid_reference>urn:nasa:pds:context:facility:observatory.paul_wild</lid_reference>
                               <reference_type>is_facility</reference_type>
                         </Internal_Reference>
                   </Observing_System_Component>
                   """

    return msg
        
def dsnContext(dsn, year, inputFile):
    msg = """ 
                   <Observing_System_Component>
                         <name>DSN Antenna DSS %s</name>
                         <type>Telescope</type>
                         <Internal_Reference>
                               <lid_reference>urn:nasa:pds:context:telescope:%s</lid_reference>
                               <reference_type>is_telescope</reference_type>
                         </Internal_Reference>
                   </Observing_System_Component>"""

    m=commnCntx(dsn)
    for d in dsn:
        dss = 'dss'+d
        if int(d) == 65 and year > 2005: dss = dss+"P"

        m = m + msg%(d, facility(dss, inputFile)) 
    
    return m

