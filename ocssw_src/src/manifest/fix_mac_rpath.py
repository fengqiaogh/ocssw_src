#! /usr/bin/env python3

# This program looks for shared library references that need to be changed to RPATH.
# It also copies libs from homebrew and X11 into OCSSW/opt/lib for seadas distribution

import subprocess
import os
import shutil

rerunNeeded = True
while rerunNeeded:
    rerunNeeded = False

    # set the rpath in opt/libs
    os.chdir(os.path.join(os.environ['LIB3_DIR'], "lib"))
    for fileName in os.listdir('.'):
        if os.path.isfile(fileName):
            if ".dylib" in fileName:
                #print (fileName)
                p = subprocess.Popen(["otool", "-D", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                parts = out.decode("utf-8").split()
                parts = parts[1].split("/")
                rpath = parts[0]
                if not "@rpath" in parts[0]:
                    name = parts[-1]
                    id = "@rpath/" + name
                    # print(fileName, id)
                    subprocess.call(["install_name_tool", "-id", id, fileName])

                p = subprocess.Popen(["otool", "-L", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                lines = out.decode("utf-8").split('\n')
                for line in lines:
                    if '/opt/homebrew/' in line or '/usr/local/' in line or '/opt/X11/lib/' in line or os.environ['LIB3_DIR'] in line or ('/' not in line and 'compatibility' in line):
                        #print('  ' + line)
                        libPath = line.split()[0]
                        parts = libPath.split('/')
                        libName = parts[-1]
                        newName = '@rpath/' + libName
                        #print('  ' + libPath + ' -> ' + newName)
                        subprocess.call(["install_name_tool", "-change", libPath, newName, fileName])
                        if os.environ['LIB3_DIR'] not in libPath:
                            if not os.path.isfile(os.environ['LIB3_DIR'] + "/lib/" + libName):
                                rerunNeeded = True
                                print("copying", libPath, os.environ['LIB3_DIR'] + "/lib" )
                                shutil.copy(libPath, os.environ['LIB3_DIR'] + "/lib")


    # set the rpath in opt/bin
    os.chdir(os.path.join(os.environ['LIB3_DIR'], "bin"))
    for fileName in os.listdir('.'):
        if os.path.isfile(fileName):
            line = subprocess.check_output(['file', fileName]).decode("utf-8")
            if "Mach-O 64-bit executable" in line:
                #print ('------' + fileName)
                p = subprocess.Popen(["otool", "-L", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                lines = out.decode("utf-8").split('\n')
                for line in lines:
                    if '/opt/homebrew/' in line or '/usr/local/' in line or '/opt/X11/lib/' in line or os.environ['LIB3_DIR'] in line:
                        print('  ' + line)
                        libPath = line.split()[0]
                        parts = libPath.split('/')
                        libName = parts[-1]
                        newName = '@rpath/' + libName
                        #print('  ' + libPath + ' -> ' + newName)
                        subprocess.call(["install_name_tool", "-change", libPath, newName, fileName])
                        if os.environ['LIB3_DIR'] not in libPath:
                            if not os.path.isfile(os.environ['LIB3_DIR'] + "/lib/" + libName):
                                rerunNeeded = True
                                print("copying", libPath, os.environ['LIB3_DIR'] + "/lib" )
                                shutil.copy(libPath, os.environ['LIB3_DIR'] + "/lib")


    # set the rpath in OCSSW/libs
    os.chdir(os.path.join(os.environ['OCSSWROOT'], "lib"))
    for fileName in os.listdir('.'):
        if os.path.isfile(fileName):
            if ".dylib" in fileName:
                #print (fileName)
                p = subprocess.Popen(["otool", "-D", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                parts = out.decode("utf-8").split()
                parts = parts[1].split("/")
                rpath = parts[0]
                if not "@rpath" in parts[0]:
                    name = parts[-1]
                    id = "@rpath/" + name
                    # print(fileName, id)
                    subprocess.call(["install_name_tool", "-id", id, fileName])

                p = subprocess.Popen(["otool", "-L", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                lines = out.decode("utf-8").split('\n')
                for line in lines:
                    if '/opt/homebrew/' in line or '/usr/local/' in line or '/opt/X11/lib/' in line or os.environ['LIB3_DIR'] in line:
                        #print('  ' + line)
                        libPath = line.split()[0]
                        parts = libPath.split('/')
                        libName = parts[-1]
                        newName = '@rpath/' + libName
                        #print('  ' + libPath + ' -> ' + newName)
                        subprocess.call(["install_name_tool", "-change", libPath, newName, fileName])
                        if os.environ['LIB3_DIR'] not in libPath:
                            if not os.path.isfile(os.environ['LIB3_DIR'] + "/lib/" + libName):
                                rerunNeeded = True
                                print("copying", libPath, os.environ['LIB3_DIR'] + "/lib" )
                                shutil.copy(libPath, os.environ['LIB3_DIR'] + "/lib")

                        
    # set the rpath in OCSSW/bin
    os.chdir(os.path.join(os.environ['OCSSWROOT'], "bin"))
    for fileName in os.listdir('.'):
        if os.path.isfile(fileName):
            line = subprocess.check_output(['file', fileName]).decode("utf-8")
            if "Mach-O 64-bit executable" in line:
                #print ('------' + fileName)
                p = subprocess.Popen(["otool", "-L", fileName], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                #print("out=",out)
                lines = out.decode("utf-8").split('\n')
                for line in lines:
                    if '/opt/homebrew/' in line or '/usr/local/' in line or '/opt/X11/lib/' in line or os.environ['LIB3_DIR'] in line or '%s/lib' % (os.environ['OCSSWROOT']) in line:
                        #print('  ' + line)
                        libPath = line.split()[0]
                        parts = libPath.split('/')
                        libName = parts[-1]
                        newName = '@rpath/' + libName
                        #print('  ' + libPath + ' -> ' + newName)
                        subprocess.call(["install_name_tool", "-change", libPath, newName, fileName])
                        if os.environ['LIB3_DIR'] not in libPath:
                            if not os.path.isfile(os.environ['LIB3_DIR'] + "/lib/" + libName):
                                rerunNeeded = True
                                print("copying", libPath, os.environ['LIB3_DIR'] + "/lib" )
                                shutil.copy(libPath, os.environ['LIB3_DIR'] + "/lib")

