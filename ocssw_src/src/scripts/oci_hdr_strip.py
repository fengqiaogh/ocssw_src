# Util script to strip OCI headers off of .L0 files

from sys import argv
from os import getcwd

cwd = getcwd()
argc = len(argv)

tmp_filename = cwd + "tmp.L0"

with open(tmp_filename, "ab") as new_file:
    for file in argv[0:]:
        out = open(file, "rb")
        new_file.write(out.read()[64:])
        out.close
    new_file.close()
