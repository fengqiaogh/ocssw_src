/*******************************************************************************
 *
 * NAME:  VcstTime.h
 *
 * DESCRIPTION: This is the header file for the VcstTime object class.
 *
 * Adapted directly from various IDPS files, InfUtil_Tim.h, InfUtil_TimApi.h,
 * etc., published by Raytheon Company.
 *
 *******************************************************************************/

#ifndef _VcstTime_h
#define _VcstTime_h

/**
 @addtogroup TimeCpp
 @{
 */

#include <map>
#include <set>
#include <string>
#include <exception>
#include <time.h>

#define MAX_LEAP 156
#define MAXBUFFER 256

/*
 *  The following are all valid error numbers for the time utility
 */
#define TIMERR_ROOT         -1600

#define TIMERR_BADDAY       TIMERR_ROOT -  3
#define TIMERR_BADHOUR      TIMERR_ROOT -  4
#define TIMERR_BADMIN       TIMERR_ROOT -  5
#define TIMERR_BADMONTH     TIMERR_ROOT -  6
#define TIMERR_BADMSEC      TIMERR_ROOT -  7
#define TIMERR_BADSEC       TIMERR_ROOT -  8
#define TIMERR_BADYEAR      TIMERR_ROOT -  9
#define TIMERR_LEAPACCESS   TIMERR_ROOT - 10
#define TIMERR_LEAPBADYEAR  TIMERR_ROOT - 11
#define TIMERR_LEAPBADMONTH TIMERR_ROOT - 12
#define TIMERR_LEAPBADDAY   TIMERR_ROOT - 13
#define TIMERR_LEAPBADLEAP  TIMERR_ROOT - 14
#define TIMERR_LEAPOPEN     TIMERR_ROOT - 15
#define TIMERR_MONTHCALC    TIMERR_ROOT - 16
#define TIMERR_BADUTC       TIMERR_ROOT - 17
#define TIMERR_UTCCOLON2    TIMERR_ROOT - 18
#define TIMERR_UTCDASH1     TIMERR_ROOT - 19
#define TIMERR_UTCDOT1      TIMERR_ROOT - 20
#define TIMERR_UTCMONTH     TIMERR_ROOT - 21
#define TIMERR_UTCT         TIMERR_ROOT - 22
#define TIMERR_NOINIT       TIMERR_ROOT - 23
#define TIMERR_BADTAI       TIMERR_ROOT - 24

/**
 * The GDT format is Gregorian Date / Time format. It is individual
 * numbers for each portion of a calendar day. To a prevent large list
 * of parameters, this structure is used.
 */
struct TimGdtStruct {
    int year;
    int month;
    int day;
    int hour;
    int min;
    int sec;
    int mill;
    int msec;
    int wkday;
    int numday;

};

/**
 * The format for the leapsecond file was taken from the Naval observatory.
 * This means that new releases can simply use a file replace. Since all
 * leap seconds start with the first second of a day, only year, month
 * and day are necessary for the definition. "Leap" is the number of leap
 * seconds for that given moment in time.
 */

struct LeapStruct {
    int year;
    int month;
    int day;
    double tai;
    int leap;
};

/**
 * The list of error messages are kept in this structural format.
 */
struct TimError {
    int errNo;
    char errmsg[MAXBUFFER];
};

/**
 * This struct ltGDT is used for less than comparisons on TimGdtStruct*
 * by STL.
 */

struct ltGDT {

    bool operator()(const TimGdtStruct* gdt1, const TimGdtStruct* gdt2) const {
        if (gdt1 == NULL) {
            if (gdt2 == NULL) {
                return false;
            }

            return true; // gdt1 == null && gdt2 !=null
        } else if (gdt2 == NULL) //gdt1 != null && gdt2 == null
        {
            return false;
        } else {
            if (gdt1->year < gdt2->year) {
                return true;
            } else if ((gdt1->year == gdt2->year)
                    && (gdt1->month < gdt2->month)) {
                return true;
            } else if ((gdt1->year == gdt2->year)
                    && (gdt1->month == gdt2->month)
                    && (gdt1->day < gdt2->day)) {
                return true;
            } else if ((gdt1->year == gdt2->year)
                    && (gdt1->month == gdt2->month) && (gdt1->day == gdt2->day)
                    && (gdt1->hour < gdt2->hour)) {
                return true;
            } else if ((gdt1->year == gdt2->year)
                    && (gdt1->month == gdt2->month) && (gdt1->day == gdt2->day)
                    && (gdt1->hour == gdt2->hour) && (gdt1->min < gdt2->min)) {
                return true;
            } else if ((gdt1->year == gdt2->year)
                    && (gdt1->month == gdt2->month) && (gdt1->day == gdt2->day)
                    && (gdt1->hour == gdt2->hour) && (gdt1->min == gdt2->min)
                    && (gdt1->sec < gdt2->sec)) {
                return true;
            } else {
                // We stop at second because we want to see if this GDT is
                // during the same second.  NOT that it is EXACTLY the same.
                return false;
            }
        }
    }
};

/**
 * This struct ltLeap is used for less than comparisons on LeapStruct*
 * by STL.
 */
struct ltLeap {

    bool operator()(const LeapStruct* leap1, const LeapStruct* leap2) const {
        if (leap1 == NULL) {
            if (leap2 == NULL) {
                return false;
            }
            return true; // leap1 == null && leap2 !=null
        } else if (leap2 == NULL) //leap1 != null && leap2 == null
        {
            return false;
        } else {
            if (leap1->tai < leap2->tai) {
                // If the TAI is the same ALL other values should be the same
                // in a consistent CDS Struct.
                return true;
            }
            return false;
        }
    }
};

/**
 * This struct ltCDS is used for less than comparisons on Int64's (CDS
 * representations by STL.
 */
struct lsCDS {

    bool operator()(const long long l1, const long long l2) const {
        return l1 < l2;
    }
};

typedef std::set<const TimGdtStruct*, ltGDT> GDTSet;
typedef std::set<unsigned long long, lsCDS> CDSSet;

/**
 * This class manages the get time, and time conversion functions.
 */

class VcstTime {
public:

    /**
     * This method instantiates a singleton object. If it is called
     * multiple times, it will return the same object pointer.
     *
     * @return Singleton instance of VcstTime.
     */
    static VcstTime* getInstance(void);

    /**
     * This method finds the Naval Observatory
     * leap second file to initialize the time utility.  NOTE: This
     * method is intended to be called once within a process space
     * and should not be called multiple time within a multi-threaded
     * environment.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int initialize(void);

    /**
     * This method contains the steps necessary to initialize
     * and set up for use the time conversions. NOTE: This
     * method is intended to be called once within a process space
     * and should not be called multiple time within a multi-threaded
     * environment.
     *
     * @param inLeapFile path name for leap second file.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int initialize(const char* inLeapFile);

    /**
     * Requests the system time in International Atomic Time (TAI) format.
     *
     * @return Current time in TAI format.
     *
     * @throw TimException time exception, if error.
     */
    double taiNow(void);

    /**
     * Converts Time (TAI) format time into Unix SystemTime
     *
     * @param Current time in TAI format.
     *
     * @return tspec pointer to timespec structure for return
     */
    void tai2Unix(double inTai, timespec*);

    /**
     * Requests the system time in International Atomic Time (TAI) format.
     *
     * @param tspec pointer to timespec structure for return
     *
     * @return Current time in TAI format.
     */
    double unix2Tai(timespec *tspec);

    /**
     * Add microseconds to TAI format, maintain leap seconds.
     *
     * @param ioIet  TAI format time to be added to
     * @param inDelta Number of Microseconds being added
     *
     * @return Current time in TAI format.
     */
    int taiAddDelta(double *ioTai, double inDelta);

    /**
     * Subtract microseconds to TAI format, maintain leap seconds.
     *
     * @param ioIet  TAI format time to be added to
     * @param inDelta Number of Microseconds being subtracted.
     *
     * @return Current time in TAI format.
     */
    int taiSubtractDelta(double *ioTai, double inDelta);

    /**
     * This method contains the algorithm to convert the date/time
     * in CCSDS Day Segmentation (CDS) format to International Atomic Time (TAI) format.
     *
     * @param inCds  Time in CDS format for conversion
     * @param outTai Place to put time in TAI Format
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int cds2Tai(unsigned long long inCds, double *outTai);

    /**
     * This method contains the algorithm to convert the International Atomic
     * Time (TAI) format to date/time in CCSDS Day Segmentation (CDS) format.
     *
     * @param inTai  Time in TAI format for conversion.
     * @param outCds Place to put time in CDS format.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int tai2Cds(double inTai, unsigned long long *outCds);

    /**
     * This method Displays a CDS structure
     *
     * @param inCds Place to put time in CDS format.
     *
     * @return Status. 0 if OK, error code if not.
     */
    int dispCds(unsigned long long inCds);

    /**
     * This method contains the algorithm to convert the date/time
     * in International Atomic Time (TAI) format to Gregorian Date Time (GDT) format.
     *
     * @param inTai  Time in TAI format for conversion
     * @param outGdt place to put GDT version.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int tai2Gdt(double inTai, TimGdtStruct *outGdt);

    /**
     * This method contains the algorithm to convert Gregorian Date
     * Time (GDT) to International Atomic Time (TAI).
     *
     * @param inGdt time in GDT format for conversion
     * @param outTai place to put TAI format time.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int gdt2Tai(TimGdtStruct *inGdt, double *outTai);

    /**
     * This method contains algorithm to convert International Atomic Time (TAI)
     * formatted time to Universal Time Coordinated (UTC) formatted time.
     *
     * @param inTai Time in TAI format for conversion
     * @param outString Place to put UTC String output.
     * @param inFormat format of UTC string:
     *           format1:  With Dashes, with Microseconds
     *           format2:  With Dashes, NO   Microseconds
     *           format3:  NO   Dashes, with Microseconds
     *           format4:  NO   Dashes, NO   Microseconds
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int tai2String(double inTai, char* outString, int inFormat = 1);

    /**
     * This method contains algorithms to convert Coordinated Universal Time
     * (UTC) formatted time to International Atomic Time (TAI)
     * formatted time.
     *
     * @param inString Time in UTC format for conversion
     * @param outTai Time converted to TAI format.
     *
     * @return Status. 0 if OK, error code if not.
     *
     * @throw TimException time exception, if error.
     */
    int string2Tai(std::string inString, double *outTai);

    /**
     * Calling this method prints the current TAI and GDT time to standard out.
     */
    void printInfo();

    /**
     * Return the system time in the INF standard format
     *
     * @return system time in the UNIX standard format
     */
    timespec getSystemTime();

    /**
     * Return a string describing the error
     *
     * @param errNo error number for search
     *
     * @return string describing error.
     */
    static std::string errorMessage(int errNo = 1);

private:

    /// The possible type of tokens that may be returned.

    enum tokenType {
        NONE, INT, FLOAT, STRING, EOS
    };

    // Constant memebers
    static const char *TOKEN_SEPARATORS;

    // Constant memebers
    static const char *TOKEN_SEPARATORS_UTC;

    /// The structure used to contain the current token.
    /// Contains a type and a value.

    struct TokenStruct {
        tokenType type;

        union tokenVal {
            char *s;
            int i;
            float f;
            long long l;
            double d;
        } val;
    };

    /**
     * Constructor()
     *
     */
    VcstTime(void);

    /**
     * Calculate the TAI without leap seconds for a GDT value.
     * ( When reading in the leap second table, leap calcs aren't practical)
     *
     * @param outTai place to put TAI format time.
     *
     * @return Status. 0 if OK, error code if not.
     */
    void calcTai(TimGdtStruct *inGdt, double *outTai);

    /**
     * Set Feb 29th, if the year give is a leap year
     *
     * @param inYear Year for leap year check
     *
     * @return status (zero for OK), < 0 for error code
     */
    int adjustLeapDay(int inYear);

    /**
     * Return number of leap days since 1958 for a given year
     *
     * @param inYear Year for calculation
     *
     * @return status (zero for OK), < 0 for error code
     */
    int calcLeapDays(int inYear);

    /**
     * Return number of leap seconds for the given TAI Time
     *
     * @param inTai time in TAI format for leap second retrieval
     * @param withLeap TRUE, tai is leapsecond adjusted, FALSE, not.
     *
     * @return status (zero for OK), < 0 for error code
     */
    int getTaiLeapSecond(double inTai, int withLeap);

    /**
     * Return whether it is an exact leap second moment
     *
     * @param inTai time in TAI format for leap second check
     *
     * @return TRUE if an exact leap second moment, FALSE if not
     */
    int isExactLeapSecond(double inTai);

    /**
     * Return the number of leap seconds for the given year/month/day
     *
     * @param inYear  year for leapsecond lookup
     * @param inMonth month for leapsecond lookup
     * @param inDay   day for leapsecond lookup
     *
     * @return status (zero for OK), < 0 for error code
     */
    int getLeapSecond(int inYear, int inMonth, int inDay);

    /**
     * Read leap second file (File Format the same as the Naval Observatory.
     *
     * @param inFilename mame of leap second file
     *
     * @return status 0 if OK
     *
     * @throw TimException time exception, if error.
     */
    int readLeapSecondFile(const char * inFilename);

    /**
     * Return the number of leap days since 1958 for a given year
     *
     * @param *inBuffer   String from file
     * @param *outYear    year, parsed from file
     * @param *outMonth   month, parsed from file
     * @param *outDay     day, parsed from file
     * @param *outLeap    leap seconds for that date
     *
     * @return status (zero for OK), < 0 for error code
     *
     * @throw TimException time exception, if error.
     */
    int parseLeapFileLine(char *inBuffer, int *outYear, int *outMonth,
            int *outDay, int *outLeap);

    /**
     * Do a normal parse of a UTC String. The format is as follows:
     * YEAR<dash>MONTH<dash>DAY<space>
     * HOUR<colon>MIN<colon>SEC<period>MILLI+MICRO<Z>
     *    Example long : 1958-01-01 00:00:00.000000000Z
     *    Example short: 1958-01-01 00:00:00Z
     *
     * @param inString is the UTC string to be parsed
     *
     * @return a GDT object with parsed fields within
     *
     * @throw TimException time exception, if error.
     */
    TimGdtStruct utcParse(std::string inString);

    /**
     * Outputs a display of the leapsecond table. (for debug)
     */
    void displayLeap(void);

    /**
     * Outputs a display of the GDT Structure
     *
     * @param *inPrompt Prompt to display
     * @param *inGdt structure to display
     *
     */
    void displayGdt(char *inPrompt, TimGdtStruct *inGdt);

    /**
     * parse a string, return substring ending before delin character
     *
     * @param *inBuffer string for search
     * @param  delin    deliniation character
     * @param *outSub   Substring returned.
     *
     */
    char *nextField(char *inBuffer, char delin, char *outSub);

    /**
     * Parses for the next token in the string myBuffer.  When the
     * String has been finished, the token returned will be EOS (End Of String)
     *
     * @param tokens The tokens to use to parse the string.
     *
     * @return TokenStruct - the next token parsed.
     */
    TokenStruct getNextToken(const char* tokens);

    /**
     * Parses for the next token in the string myBuffer.  myBuffer is
     * initialized to myBuffer before obtaining the first token.  If the
     * String has been finished, the token returned will be EOS (End Of String)
     *
     * @param parseString The string to be parsed.
     * @param tokens The tokens to use to parse the string.
     *
     * @return TokenStruct - the next token parsed.
     */
    TokenStruct getNextToken(char * parseString, const char* tokens);

    /**
     * Returns true if the given GDT is on or during a Negative Leap Second
     *
     * @param inGdt time in GDT format for conversion
     *
     * @return true if the given GDT is on or during a Negative Leap Second,
     *         false otherwise.
     */
    bool isNegativeLeapSecond(TimGdtStruct *inGdt);

    /**
     * Returns true if the given CDS is on or during a Negative Leap Second
     *
     * @param inCDS time in CDS format for conversion
     *
     * @return true if the given CDS is on or during a Negative Leap Second
     *         false otherwise.
     */
    bool isNegativeLeapSecondCDS(long long inCDS);

    /**
     * Initializes the Negative Leaps.
     *
     * @throw TimException
     */
    void initNegativeLeaps(void);

private:
    // Members

    /// Structure array of leap second information
    struct LeapStruct myLeapList[MAX_LEAP + 1];

    /// Number of leap seconds in the leap second structure array.
    int myNleap;

    // The pointer to the string currently being parsed.
    char *myBuffer;

    /// List of month strings for UTC conversion
    static char myMonthStr[13][13];

    /// Number of days in each month (Feb changes on leap years)
    static int myDayInMonth[13];

    /// Whether init has been called flag.
    bool myInit;

    // A Set of TimGdtStructs representing negative Leaps
    GDTSet* myGDTNegativeLeaps;

    // A Set of long long (CDS representations) Leaps representing negative Leaps
    CDSSet* myCDSNegativeLeaps;

    // The hash map of month to integer conversion.
    std::map<std::string, int> myMonthMap;
};

/**
 * Generic error string for this exception
 */
extern const std::string EXCEPTION_STRING;

/**
 * Error message for Invalid Day
 */
extern const std::string ERR_BADDAY_MSG;

/**
 * Error message for Invalid Hour
 */
extern const std::string ERR_BADHOUR_MSG;

/**
 * Error message for Invalid Minute
 */
extern const std::string ERR_BADMIN_MSG;

/**
 * Error message for Invalid Month
 */
extern const std::string ERR_BADMONTH_MSG;

/**
 * Error message for Invalid Microsecond
 */
extern const std::string ERR_BADMSEC_MSG;

/**
 * Error message for Invalid Second
 */
extern const std::string ERR_BADSEC_MSG;

/**
 * Error message for Invalid Year
 */
extern const std::string ERR_BADYEAR_MSG;

/**
 * Error message for Leap Second file access
 */
extern const std::string ERR_LEAPACCESS_MSG;

/**
 * Error message for Leap Second read failed
 */
extern const std::string ERR_LEAPREAD_MSG;

/**
 * Error message for Leap Second file invalid day
 */
extern const std::string ERR_LEAPBADDAY_MSG;

/**
 * Error message for Leap Second file invalid leap second
 */
extern const std::string ERR_LEAPBADLEAP_MSG;

/**
 * Error message for Leap Second file invalid month
 */
extern const std::string ERR_LEAPBADMONTH_MSG;

/**
 * Error message for Leap Second file invalid year
 */
extern const std::string ERR_LEAPBADYEAR_MSG;

/**
 * Error message for Leap Second file open failed
 */
extern const std::string ERR_LEAPOPEN_MSG;

/**
 * Error message for MONTH CALC LOOP FAILED
 */
extern const std::string ERR_MONTHCALC_MSG;

/**
 * Error message for Invalid UTC string
 */
extern const std::string ERR_BADUTC_MSG;

/**
 * Error message for Invalid Second Colon in UTC string
 */
extern const std::string ERR_UTCCOLON2_MSG;

/**
 * Error message for No First Dash in UTC string
 */
extern const std::string ERR_UTCDASH1_MSG;

/**
 * Error message for Invalid microseconds dot in UTC string
 */
extern const std::string ERR_UTCDOT1_MSG;

/**
 * Error message for Invalid Month in UTC string
 */
extern const std::string ERR_UTCMONTH_MSG;

/**
 * Error message for Invalid 'T' delimiter in UTC string
 */
extern const std::string ERR_UTCT_MSG;

/**
 * Error message for Init must be called before this function
 */
extern const std::string ERR_NOINIT_MSG;

/**
 * Error message for negative TAI input
 */
extern const std::string ERR_BADTAI_MSG;

/**
 * Error message for an TAI value that is greater than the upper limit
 * of supported values
 */
extern const std::string MAX_TAI_VALUE_EXCEEDED_MSG;

/**
 * Error message for put / output exception
 */
extern const std::string IOEXCEPTION_Msg;

/**
 * TimException class is the main class for all
 * exceptions defined in the Time Utility.
 */

class TimException : public std::exception {
public:

    /**
     * Generic error code for this exception
     */
    static const int EXCEPTION_CODE;

    /**
     * Error code if no name was found.
     */
    static const int NAME_NOT_FOUND;

    /**
     * Error code for Invalid Day
     */
    static const int ERR_BADDAY;

    /**
     * Error code for Invalid Hour
     */
    static const int ERR_BADHOUR;

    /**
     * Error code for Invalid Minute
     */
    static const int ERR_BADMIN;

    /**
     * Error code for Invalid Month
     */
    static const int ERR_BADMONTH;

    /**
     * Error code for Invalid Microsecond
     */
    static const int ERR_BADMSEC;

    /**
     * Error code for Invalid Second
     */
    static const int ERR_BADSEC;

    /**
     * Error code for Invalid Year
     */
    static const int ERR_BADYEAR;

    /**
     * Error code for Leap Second file access
     */
    static const int ERR_LEAPACCESS;

    /**
     * Error code for Leap Second read failed
     */
    static const int ERR_LEAPREAD;

    /**
     * Error code for Leap Second file invalid day
     */
    static const int ERR_LEAPBADYEAR;

    /**
     * Error code for Leap Second file invalid leap second
     */
    static const int ERR_LEAPBADMONTH;

    /**
     * Error code for Leap Second file invalid month
     */
    static const int ERR_LEAPBADDAY;

    /**
     * Error code for Leap Second file invalid year
     */
    static const int ERR_LEAPBADLEAP;

    /**
     * Error code for Leap Second file open failed
     */
    static const int ERR_LEAPOPEN;

    /**
     * Error code for MONTH CALC LOOP FAILED
     */
    static const int ERR_MONTHCALC;

    /**
     * Error code for Invalid UTC string
     */
    static const int ERR_BADUTC;

    /**
     * Error code for Invalid Second Colon in UTC string
     */
    static const int ERR_UTCCOLON2;

    /**
     * Error code for No First Dash in UTC string
     */
    static const int ERR_UTCDASH1;

    /**
     * Error code for Invalid microseconds dot in UTC string
     */
    static const int ERR_UTCDOT1;

    /**
     * Error code for Invalid Month in UTC string
     */
    static const int ERR_UTCMONTH;

    /**
     * Error code for Invalid 'T' delimiter in UTC string
     */
    static const int ERR_UTCT;

    /**
     * Error code for Init must be called before this function
     */
    static const int ERR_NOINIT;

    /**
     * Error code for negative TAI input
     */
    static const int ERR_BADTAI;

    /**
     * Error code for an TAI value that is greater than the upper limit
     * of supported values
     */
    static const int MAX_TAI_VALUE_EXCEEDED;

    /**
     * Error code for Input / output exception
     */
    static const int IOEXCEPTION;

    /**
     *    Default Constructor.
     *
     */
    TimException(void);

    /**
     *    String and code Constructor. This is the main one used
     *    by Inf_Exception.
     *
     * @param errorString error string for this exception
     * @param errorCode error code for this exception
     *
     */
    TimException(const std::string& errorString, int errorCode);

    /**
     *    Destructor
     *
     */

    ~TimException(void) throw () {
        if (myCause) {
            delete myCause;
            myCause = 0;
        }
    }
    ;

    inline int getErrorCode() const {
        return myErrorCode;
    }

    inline std::string getErrorString() const {
        return myErrorString;
    }

private:

    /**
     * Specific error code for this exception.  Appropriate error codes
     * are documented in the INF Error Handling charts.
     */
    int myErrorCode;

    /**
     * Specific description for this exception
     */
    std::string myErrorString;

    /**
     * Optional cause of this exception.
     */
    std::exception* myCause;
};

/**
 * TimCDS union.
 */

union TimCDS {
public:
    unsigned long long ll;
    unsigned char uc[8];

    //getters
    unsigned short getDays();
    unsigned int getMillis();
    unsigned short getMicros();

    // Setters
    void setDays(unsigned short);
    void setMillis(unsigned int);
    void setMicros(unsigned short);

    /**
     *    Stuff day, millisecond and microsecond into a byte buffer
     *
     * @param inDay Number of days since 1958
     * @param inMsec Number of Microseconds
     * @param *outCds where to put cds format output
     *
     * @return Status. 0 if OK, error code if not.*/

    static int makeCds(unsigned int inDay, unsigned long long inMsec,
            unsigned long long *outCds);

    /**
     *    Parse day and microseconds out of byte buffer known as CDS.
     *
     * @param inCds byte buffer for display.
     * @param outDay Number of days since 1958
     * @param outMsec Number of Microseconds
     *
     * @return Status. 0 if OK, error code if not.
     */

    static int parseCds(unsigned long long inCds, unsigned int *outDay,
            unsigned long long *outMsec);

    /**
     *    Display the internals of the byte buffer known as CDS.
     *
     * @param inCds - byte buffer for display.
     *
     * @return Status. 0 if OK, error code if not.
     */

    static int displayCds(unsigned long long inCds);

};

#endif

