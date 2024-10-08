# PEX for ODPS

In order to distribute python code to ODPS **without** requiring the production
system to carry an over-bloated conda environment, we build reasonably portable
PEX (Python EXecutable) binaries.  The following describes how this is done.

## The Python environment

Currently, ODPS runs Ubuntu 20.04, which includes python 3.8.  The build
environment needs to match that as closely as possible.  

### virtual environment
Set up a python virtual environment with the depednecies for all scripts included.

```
$ python3 -m venv <path to>/odps
$ source <path to>/odps/bin/activate
```

Once in the appropriate environment, time to make sure the required packages are
included.  The current list is relatively small (some of these require others):

 * affine
 * astropy
 * GDAL
 * netCDF4
 * numpy
 * pandas
 * pex
 * requests
 * scipy
 * MetPy
 * setuptools

-----------OLD----------
...and until we can get rid of merge_met_reanalysis:
 * pyhdf (install with conda install -c conda-forge pyhdf)
-----------OLD----------

Most of these can be installed individually, e.g.:
```
$ pip install netCDF4
```
but...you could just use pip and the odps_requirements.txt file:
```
pip install -r $OCSSWROOT/src/scripts/pex/odps_requirements.txt
```

### CMake to run PEX
For PEX to work, you need:
 * a valid setup.py file
 * a requirements.txt file
 * a CMakeLists.txt file that defines the "build"

Best to look at the ones already created for examples ....but here's some basics:

#### setup.py
The setup.py file should be generated by CMake using a setup.py.in file:

```
from distutils.core import setup

setup(name='seabass2L1B',
      version='${PACKAGE_VERSION}',
      package_dir={ '': '${CMAKE_CURRENT_BINARY_DIR}' },
      packages=['seabass'],
      scripts=['$ENV{OCSSWROOT}/src/matchup_tools/seabass2L1B.py']
)
```
#### requirements.txt
This file contains the list of python packages required.  For example, if you need netCDF4:

```
$ cat requirements.txt
netCDF4
numpy

(We're not using 'versions' in the requirements.txt, as we want to build multi-python compatible PEX binaries)
```

You could use pipreqs to generate this file:
```
pip install pipreqs
pipreqs /path/to/project
```

#### The wheel(house)
In order to enusre the PEX binaries get the versions of the python packages we need, we will use our own wheel files.
We'll tell pex to not go off an pull down the dependancies (with --no-index -f <path to wheel directory>)

The wheels are generated ala:
```
python -m pip wheel  --wheel-dir=<path to wheel directory> <packagename[==version]>
```

It *should* work to use the same odps_requiremnts.txt file:
```
python -m pip wheel  --wheel-dir=<path to wheel directory> -r $OCSSWROOT/src/scripts/pex/odps_requirements.txt
```

This is done for each python venv to build against.


#### CMakeLists.txt
The CMakeLists.txt file in the root directory ($OCSSWROOT/src/scripts/pex) defines the binaries to create.
It is independent from the core OCSSW build process.  As an added level of "make-sure-this-is-what-you-want",
cmake needs to be called with -DBUILD_PEX=1,  e.g.:

```
$ mkdir $OCSSWROOT/src/scripts/pex/build
$ cd $OCSSWROOT/src/scripts/pex/build
$ cmake ..
$ make install
```

The CMakeLists.txt file in the individual program subdirectories uses `configure_file(<IN> <OUT>)` to expand the variables in the setup.py.in file.
/accounts/swbaile1/ocssw/bin/l1bgen_hico
If your code requires a module, you can copy the required files into a build version:
```
    set(MODULES     "${CMAKE_CURRENT_BINARY_DIR}/seabass")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    file(MAKE_DIRECTORY ${MODULES})
    configure_file("${SCRIPTSDIR}/SB_support.py" ${MODULES} COPYONLY)
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E touch "${MODULES}/__init.py__")
```

Running PEX requires the `execute_process` command:
```
install(CODE "execute_process(COMMAND pex ${CMAKE_CURRENT_BINARY_DIR} -v
-r ${PYREQUIRE}  -o ${OUTPUT} -c seabass2L1B.py --no-compile --disable-cache)")
```

Here `-r` points to the requirements.txt files; `-o` is the output PEX binary; `-c` is the "console" script name for the "entry point" - the bit to run

