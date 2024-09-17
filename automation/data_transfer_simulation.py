"""
Script to simulate file transfer to test directory watcher
Simply transfers files from one location to another (making a copy).
"""

import os
import shutil
import time

source_path = 'G:\\My Drive\\Job\\IGP\\Scanner data\\data dia 16\\'
destination_path = 'G:\\My Drive\\Job\\IGP\\Scanner data\\watch directory\\'

scan_dirs = os.listdir(source_path)

for scan_dir in scan_dirs:
    scan_path = os.path.join(source_path, scan_dir)
    destination_scan_path = os.path.join(destination_path, scan_dir)
    if not os.path.exists(destination_scan_path):
        os.mkdir(destination_scan_path)

    files = os.listdir(scan_path)

    for file in files:
        destination_file_path = os.path.join(destination_scan_path, file)
        source_file_path = os.path.join(scan_path, file)

        shutil.copyfile(source_file_path, destination_file_path)

        print('Transferring file: {}'.format(file))

        time.sleep(1)

    complete_file = os.path.join(destination_scan_path, 'complete.txt')
    open(complete_file, 'w').close()


