mv context-pds4 context-pds4.old 
wget -r -nH --cut-dirs=2 --no-parent --no-check-certificate --reject="index.html*" -e robots=off https://pds.nasa.gov/data/pds4/context-pds4/
python read_context-pds4.py
