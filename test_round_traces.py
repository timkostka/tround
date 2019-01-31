"""
Perform unit tests on round_traces.py

"""

import os
import sys
import glob

# delete backup files

# execute as a script
if __name__ == "__main__":
    # output Python version
    os.system('python --version')
    # output test directory
    test_directory = os.path.dirname(os.path.realpath(sys.argv[0]))
    # delete backup files
    for filename in glob.glob(os.path.join(test_directory, 'backup_*.brd.*')):
        os.remove(filename)
