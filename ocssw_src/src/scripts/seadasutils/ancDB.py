import sqlite3
import re
from seadasutils.anc_data_types import AncillaryDataTypes

class ancDB:



    def __init__(self, dbfile=None, local=False):
        """A small set of functions to generate, update, and read from a local SQLite database of ancillary
        file information"""
        self.dbfile = dbfile
        self.local = local
        self.conn = None
        self.cursor = None

    def openDB(self):
        """
        Open connection to the ancillary DB and initiate a cursor
        """
        conn = sqlite3.connect(self.dbfile, timeout=30)
        self.conn = conn
        c = conn.cursor()
        c.execute("""PRAGMA foreign_keys = ON""")
        self.cursor = c
        return

    def closeDB(self):
        """
        Close the DB connection, committing changes.
        """
        conn = self.conn
        cursor = self.cursor
        conn.commit()
        cursor.close()

    def create_db(self):
        """
        Create the ancillary DB
        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor
        # Create  satfiles table
        c.execute(
            """CREATE TABLE IF NOT EXISTS satfiles
            (satid INTEGER PRIMARY KEY,
            filename TEXT ,
            starttime TEXT,
            stoptime TEXT,
            status INTEGER,
            attephstat INTEGER)"""
        )

        # Create  ancfiles table
        c.execute(
            """CREATE  TABLE IF NOT EXISTS ancfiles
            (ancid INTEGER PRIMARY KEY,
            filename TEXT ,
            path TEXT  ,
            type TEXT)"""
        )

        # Create  satancinfo table
        c.execute(
            """CREATE  TABLE IF NOT EXISTS satancinfo
            (satid INTEGER  ,
            ancid INTEGER  ,
            optimal INTEGER,
            FOREIGN KEY(satID) REFERENCES satfiles(satid),
            FOREIGN KEY(ancID) REFERENCES ancfiles(ancid) ON DELETE RESTRICT)"""
        )

    def insert_record(
        self,
        satfile=None,
        starttime=None,
        stoptime=None,
        dbstat=0,
        ancfile=None,
        ancpath=None,
        anctype=None,
        atteph=False,
    ):
        """
        Insert record into ancillary DB
        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor
        satid = self.check_file(satfile, starttime=starttime)
        ancid = self.check_file(ancfile, anctype=anctype)

        if satid is None:
            inputdbstat = dbstat
            attephstat = -1
            if atteph:
                attephstat = dbstat
                inputdbstat = -1

            c.execute(
                "INSERT INTO satfiles VALUES (NULL,?,?,?,?,?)",
                [satfile, starttime, stoptime, inputdbstat, attephstat],
            )
            self.conn.commit()
            satid = ancDB.check_file(self, satfile, starttime=starttime)

        else:
            if atteph:
                c.execute(
                    """UPDATE satfiles SET attephstat = ?
                                 WHERE satid = ?""",
                    [dbstat, satid],
                )
            else:
                c.execute(
                    """UPDATE satfiles SET status = ?
                                 WHERE satid = ?""",
                    [dbstat, satid],
                )

            self.conn.commit()

        if ancid is None:
            c.execute(
                "INSERT INTO ancfiles VALUES (NULL,?,?,?)", [ancfile, ancpath, anctype]
            )
            self.conn.commit()
            ancid = ancDB.check_file(self, ancfile, anctype=anctype)

        opt = self.check_dbrtn_status(dbstat, anctype)

        result = c.execute(
            "SELECT * from satancinfo where satid = ? and ancid = ?", [satid, ancid]
        )
        r = result.fetchone()

        if r is None:
            c.execute("INSERT INTO satancinfo VALUES (?,?,?)", [satid, ancid, opt])

    def delete_record(self, filename, anctype=None, starttime=None):
        """
        Deletes records from ancillary DB

        Given an anctype, an ancillary file and its dependent satfiles will be deleted.
        Given a starttime and no anctype, or only a filename, a satfile will be deleted, and its ancillary files will be deleted only if no other satfiles depend on them.

        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor
        conn = self.conn

        if anctype:  # Deleting an ancillary file
            ancid = self.check_file(filename, anctype=anctype)
            c.execute("DELETE from satancinfo WHERE ancid = ?", [ancid])

            satids = conn.execute(
                "SELECT satid FROM satancinfo WHERE ancid = ?", [ancid]
            ).fetchall()
            for satid in satids:  # Purge satfiles that depended on this ancfile
                c.execute("DELETE from satfiles WHERE satid = ?", [satid[0]])

            c.execute("DELETE from ancfiles WHERE ancid = ?", [ancid])

        else:  # Deleting a satfile
            satid = self.check_file(filename, starttime=starttime)

            ancids = conn.execute(
                "SELECT ancid FROM satancinfo WHERE satid = ?", [satid]
            ).fetchall()

            c.execute("DELETE from satancinfo WHERE satid = ?", [satid])
            c.execute("DELETE from satfiles WHERE satid = ?", [satid])

            # Purge associated ancfiles if and only if they're not needed by other satfiles
            for ancid in ancids:
                remaining_deps = conn.execute(
                    "SELECT COUNT(*) FROM satancinfo WHERE ancid = ?", [ancid[0]]
                ).fetchone()[0]

                if remaining_deps == 0:
                    c.execute("DELETE from ancfiles WHERE ancid = ?", [ancid[0]])

        conn.commit()
    
    def isAncillaryTypeValid(self, ancillaryType):
        if re.search("\d$", ancillaryType):
            ancillaryType = ancillaryType[0 : len(ancillaryType) - 1]
        if ancillaryType not in self.ancillaryTypeCheck:
            return False
        return True

    def check_dbrtn_status(self, dbstat, anctype):
        
        # grab ancillary types that can have multiple files and strip the number
        # ie. MET1, MET2, etc. to MET
        if re.search("\d$", anctype):
            anctype = anctype[0 : len(anctype) - 1]
        # anctype exists in the lookup table. Check if it is missing and report it
        if dbstat & AncillaryDataTypes.bitwiseErrorValues[anctype]:
            return 0
        else:
            return 1

    def check_file(self, filename, anctype=None, starttime=None):
        """
        Check database for existing file, return ID if exists
        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor

        table = "satfiles"
        id = "satid"
        if anctype is None:
            if filename:
                query = " ".join(
                    [
                        "select",
                        id,
                        "from",
                        table,
                        "where filename =",
                        '"' + filename + '"',
                    ]
                )
            else:
                query = " ".join(
                    [
                        "select",
                        id,
                        "from",
                        table,
                        "where starttime =",
                        '"' + starttime + '"',
                    ]
                )

        else:
            table = "ancfiles"
            id = "ancid"
            if filename:
                query = " ".join(
                    [
                        "select",
                        id,
                        "from",
                        table,
                        "where filename =",
                        '"' + filename + '"',
                        " and type = ",
                        '"' + anctype + '"',
                    ]
                )
            else:
                return None

        result = c.execute(query)
        r = result.fetchone()

        if r is None:
            return None
        else:
            if len(r) > 1:
                print(
                    "more than one entry for this starttime - this may be a problem.?"
                )
            return r[0]

    def get_status(self, filename, atteph=False, starttime=None):
        """
        Check the stored database return status
        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor
        if atteph:
            if filename:
                query = " ".join(
                    [
                        "select attephstat from satfiles where filename =",
                        '"' + filename + '"',
                    ]
                )
            else:
                query = " ".join(
                    [
                        "select attephstat from satfiles where starttime =",
                        '"' + starttime + '"',
                    ]
                )
        else:
            if filename:
                query = " ".join(
                    [
                        "select status from satfiles where filename =",
                        '"' + filename + '"',
                    ]
                )
            else:
                query = " ".join(
                    [
                        "select status from satfiles where starttime =",
                        '"' + starttime + '"',
                    ]
                )

        result = c.execute(query)
        r = result.fetchone()

        if r is None:
            return None
        else:
            return r[0]

    def get_filetime(self, filename, starttime=None):
        """
        return the stored file start and stop times
        """
        if self.conn is None:
            print("No connection to database!")
            return 110

        c = self.cursor
        if filename:
            query = " ".join(
                [
                    "select starttime,stoptime from satfiles where filename =",
                    '"' + filename + '"',
                ]
            )
        else:
            query = " ".join(
                [
                    "select starttime,stoptime from satfiles where starttime =",
                    '"' + starttime + '"',
                ]
            )

        result = c.execute(query)
        r = result.fetchone()
        return [r[0], r[1]]

    def get_ancfiles(self, filename, atteph=False, starttime=None):
        """
        Return the ancillary files associated with a given input file
        """
        import os

        if self.conn is None:
            print("No connection to database!")
            return None

        c = self.cursor

        satID = self.check_file(filename, starttime=starttime)
        if satID is None:
            return None

        filehash = {}
        result = c.execute(
            "SELECT a.type, a.path, a.filename from ancfiles a, satancinfo s where a.ancid = s.ancid and s.satid = ?",
            [satID],
        )
        for row in result:
            anctype = row[0]
            if atteph and not re.search("(att|eph)", anctype, re.IGNORECASE):
                continue
            elif not atteph and re.search("(att|eph)", anctype, re.IGNORECASE):
                continue

            filehash[row[0]] = os.path.join(row[1], row[2])

        return filehash


if __name__ == "__main__":
    import os
    import tempfile

    # Create a temporary database file
    with tempfile.NamedTemporaryFile(suffix=".sqlite.db", delete=False) as temp_db:
        db_path = temp_db.name

    try:
        # Initialize the database
        db = ancDB(dbfile=db_path)
        db.openDB()
        db.create_db()

        print("Testing ancDB functions:")

        # Test insert_record
        print("\nInserting records...")
        db.insert_record(
            satfile="A2002365234500.L1A_LAC",
            starttime="2002365234500",
            stoptime="2002365235000",
            ancfile="N200236518_MET_NCEPN_6h.hdf",
            ancpath="/Users/Shared/python/OCSSW_Scripts",
            anctype="met1",
        )
        db.insert_record(
            satfile="A2002365234500.L1A_LAC",
            starttime="2002365234500",
            stoptime="2002365235000",
            ancfile="N200300100_MET_NCEPN_6h.hdf",
            ancpath="/Users/Shared/python/OCSSW_Scripts",
            anctype="att1",
            atteph=True,
        )

        # Insert another satellite file that uses the same ancillary file
        db.insert_record(
            satfile="B2002365234500.L1A_LAC",
            starttime="2002365234500",
            stoptime="2002365235000",
            ancfile="N200236518_MET_NCEPN_6h.hdf",
            ancpath="/Users/Shared/python/OCSSW_Scripts",
            anctype="met1",
        )

        # Test check_file
        print("\nChecking files...")
        print("Satellite file A ID:", db.check_file("A2002365234500.L1A_LAC"))
        print("Satellite file B ID:", db.check_file("B2002365234500.L1A_LAC"))
        print(
            "Ancillary file ID:",
            db.check_file("N200236518_MET_NCEPN_6h.hdf", anctype="met1"),
        )

        # Test case: Ensure deletion of a satfile only deletes the ancfile if no other satfiles depend on it
        print("\nTesting deletion of satfile...")
        db.delete_record(filename="A2002365234500.L1A_LAC")
        print(
            "Checking if satellite file A still exists:",
            db.check_file("A2002365234500.L1A_LAC"),
        )
        print(
            "Checking if satellite file B still exists:",
            db.check_file("B2002365234500.L1A_LAC"),
        )
        print(
            "Checking if shared ancillary file still exists:",
            db.check_file("N200236518_MET_NCEPN_6h.hdf", anctype="met1"),
        )
        print(
            "Checking if unshared ancillary file still exists:",
            db.check_file("N200300100_MET_NCEPN_6h.hdf", anctype="att1"),
        )

        # Test case: Ensure deletion of an ancillary file deletes all satfiles that reference it
        print("\nTesting deletion of ancillary file...")
        db.delete_record(filename="N200236518_MET_NCEPN_6h.hdf", anctype="met1")
        print(
            "Checking if ancillary file still exists:",
            db.check_file("N200236518_MET_NCEPN_6h.hdf", anctype="met1"),
        )
        print(
            "Checking if satellite file B still exists:",
            db.check_file("B2002365234500.L1A_LAC"),
        )

        # Close the database connection
        db.closeDB()

    finally:
        # Clean up the temporary database file
        os.unlink(db_path)

    print("\nTests completed.")
