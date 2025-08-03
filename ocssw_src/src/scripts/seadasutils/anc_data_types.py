import re

# Holds information about valid ancillary data types and their bitwise
# values that indicate if the file is missing.
# It also holds regular expression strings to be used in anc_utils.py
# to parse for valid ancillary data types.
#
# You do not need to initialize this class to get the info.
# Instead, just import and then call the variables to get the information 
class AncillaryDataTypes:

    """
    Check the database return status and see if any files are missing.
    Also used to check if ancillary file types are valid based in the dict's keys.

    findweb() can ignore any new ancillary file types that are not in the list. 
    This is so getanc wont crash when the ancillary API gets updated. 
    Ancillary files are checked to see if any are missing using bitwise values.
    When new types are added, the bitwise dict will not have the new type's key. 

    DB return status bitwise values:
        all bits off means all is well in the world
        value of -1 means have not checked for ancfiles yet
        Ancillary:
            bit 0 - missing one or more MET
            bit 1 - missing one or more OZONE
            bit 2 - missing SST
            bit 3 - missing NO2
            bit 4 - missing ICE
            bit 5 - missing GEO
            bit 6 - missing one or more AER
            ...
            bit n - add new type
        Attitude-Ephemeris
            bit 0 - predicted attitude selected
            bit 1 - predicted ephemeris selected
            bit 2 - no attitude found
            bit 3 - no ephemeris found
            bit 4 - invalid mission
            ...
            bit n - add new type
    """
    bitwiseErrorValues = {
        "atm": 1,
        "met": 1,       # bit 0 - 00 0001
        "ozone": 2,
        "sstfile": 4,
        "no2file": 8,
        "icefile": 16,  
        "geo": 32,
        "aer": 64,      # bit 6 - 10 00000
        # atteph
        "att": 1,
        "eph": 2,
        # aquarius
        "sssfile": 32,
        "xrayfile": 64,
        "scat": 128,
        "tecfile": 256,
        "swhfile": 512,
        "frozenfile": 1024,
        "geosfile": 2048,
        "argosfile": 4096,
        "sif": 8192,  # sif_file
        "pert": 16384,  # l2_uncertainties_file
        "sssmatchup": 32768,  # sss_matchup_file
        "rim_file": 65536,
    }


    # Add additional types to grab if new ones gets added
    # Current anc types the api returns and parses
    validApiAncTypesRegEx = "met|ozone|sst|ice|aer|geo"

    validApiAttEphTypesRegEx = "att|eph"

    
    # Class method so you do not need an instance of this class to use
    # Checks if the passed in ancillary type is in the lookup table for
    # missing files
    @staticmethod
    def isValidType(ancType):

        # lookup table only has the base name (ie. mat for mat1, mat2, etc)
        # extract the base name for names that end in numbers and check the table
        if (re.search("\d$", ancType)):
            ancType = ancType[0: len(ancType)-1]
        
        if ancType not in AncillaryDataTypes.bitwiseErrorValues:
            return False
        return True
    

    # performs bitwise "or" on the current database status with the current
    # missing ancillary type file. Flip the bit value to indicate it's missing
    # ie: No missing files          ==  00 0000     ==  db_status = 0
    # ie: Missing met               ==  00 0001     ==  db_status = 1
    # ie: Missing met, ozone        ==  00 0011     ==  db_status = 3
    # ie: Missing met, ozone, aer   ==  10 0011     ==  db_status = 35
    @staticmethod
    def updateMissingFileStatus(currDatabaseStatus, ancType, isAttEph):
        
        # NOT attitude ephemeris
        # attitude ephemeris has its own bitwise values so read correct ancillary files only
        if not isAttEph:
            if re.search('met', ancType):
                return  currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["met"]       # 0000001
            if re.search('ozone', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["ozone"]      # 0000010
            if re.search('sst', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["sstfile"]    # 0000100
            if re.search('no2', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["no2file"]
            if re.search('ice', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["icefile"]
            if re.search('geo', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["geo"]
            if re.search('aer', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["aer"]
        
        # IS attitude ephemeris
        else:
            if re.search('att', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["att"]
            if re.search('eph', ancType):
                return currDatabaseStatus | AncillaryDataTypes.bitwiseErrorValues["eph"]