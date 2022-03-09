# rs_migration_tools
Tools to migrate radio science data from PDS3 standards to PDS4 standards.


Currently, 5 ancillary data types and 4 raw data types are supported:

########################################################################
The following ancillary data types are supported:
########################################################################
1. SPK: Spacecraft ephemeris (transfer file formate)
2. CKF: Spacecraft attitude (transfer file formate)
3. TRK-2-21: Earth Orientation Parameter (EOP)
4. TRK-2-23: Media Calibration Interface (Ionosphere)
5. TRK-2-23: Media Calibration Interface (Troposphere)

########################################################################
The following raw data types are supported:
########################################################################
1. TRK-2-25: Closed-loop Archival Data (also known as TDF or ATDF)
2. TRK-2-18: Closed-loop Orbital Data File (ODF)
3. TRK-2-34: Closed-loop DSN Tracking System Data Archival Format (TNF)
4. 0159-SCIENCE: Open-loop Radio Science Receiver (RSR)

########################################################################
Requirements
########################################################################
Language: Python 3.0 or above
Docgen: The igpp-docgen document generator tool based on Apache Velocity
Java: With Java Runtime Environment (JRE)
OS: Will run on any platform with a supported Java Runtime Environment (JRE)
Input data: PDS3 RS data with PDS3 label 


########################################################################
Usuage
########################################################################

Python createPDS4Label.py -i inputDataFile -v velocityTemepate

where, the inputDataFile can be any supported input file, and the velocityTemepate 
can be any predefined velocity template (example templates are provided in the 
velocityTemplates folder).

Running example:

python createPDS4Label.py -i ../sampleData/odf/s15digs2005_283_0900x25mv1.odf -v ../velocityTemplates/odf/velocityTemplate.vm

Output should looks like this:

          Reading And Unpacking Odf File: [████████████████████] 100% Time taken: 0 Hr: 0 Min: 10 Sec 

              Updating Velocity Template: [████████████████████] 100% Time taken: 0 Hr: 0 Min: 0 Sec 

                      Writing Pds4 Label: [████████████████████] 100% Time taken: 0 Hr: 0 Min: 1 Sec 


########################################################################
Support
########################################################################
If you have questions, comments, or would like to report a bug, please contact:

Ashok Kumar Verma, Ph.D
University of California, Los Angeles
Email: ashokverma@ucla.edu
