import os
import shutil
import sys
import ebltable.ebl_from_model as ebl

# 1. Locate where ebltable was installed
package_dir = os.path.dirname(ebl.__file__)
data_dir_source = os.path.join(package_dir, 'data')

# Target filename
filename = 'tau_dominguez11.fits'
source_file = os.path.join(data_dir_source, filename)

# Destination in project
dest_dir = "data/ebl_models"
dest_file = os.path.join(dest_dir, filename)

print(f"--- EBL File Extractor ---")
print(f"Searching in: {source_file}")

if os.path.exists(source_file):
    # Create destination folder if it doesn't exist
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
        
    # Copy file
    shutil.copy2(source_file, dest_file)
    print(f"\n[SUCCESS] File copied to: {dest_file}")
    print("Now you can run utils/ebl.py normally!")
else:
    print(f"\n[ERROR] File not found in the installed package.")
    print(f"Package data folder content: {os.listdir(data_dir_source)}")
    
    # Plan B: Try to find any .fits file that looks like Dominguez
    print("\n Trying smart search..")
    for f in os.listdir(data_dir_source):
        if 'dominguez' in f.lower() and 'tau' in f.lower():
            print(f"Alternative find: {f}")
            shutil.copy2(os.path.join(data_dir_source, f), os.path.join(dest_dir, f))
            print(f"Copied as substitute. Update utils/ebl.py with the name: {f}")
