#!/usr/bin/env python3

'''
# check PACE telemetry limits after they are ingested into pace_tlm DB
# program returns 100 if there are telemetries alarms, i.e. over limits
# Liang Hong, 10/10/2023
# Corrine Rojas, 10/30/2023: Changed verbosity to be the default for output txt file
# Corrine Rojas, 07/31/2025: Replaced 'd.populates' to 'd.id' column in sql query for NoneType format handling
'''

import argparse
import os
import sys
from datetime import datetime

import pyodbc

__version__ = '0.1.3 (2025-07-31)'
alert_types = {'green': 0, 'yellow': 1, 'red': 2, '0': 'green', '1': 'yellow', '2': 'red'}


def find_out_of_limits(hsk_file, conn, c, str_alert):
    if str_alert == 'red':
        # beyond red limits
        str_limit_check = (
            " (cast(v.value as float) > l.redhigh or cast(v.value as float) < l.redlow) "
        )
    elif str_alert == 'yellow':
        # between yellow and red limits
        str_limit_check = " ((cast(v.value as float) <= l.redhigh and cast(v.value as float) > l.yellowhigh) or \
        (cast(v.value as float) >= l.redlow and cast(v.value as float) < l.yellowlow)) "

    sql = (
        "select d.id,b.value, a.time_val, a.time_end from (select  v.var as var, min(v.time_val) \
    as time_val, max(v.time_val) as time_end from hsk_values v join hsk_limit l on l.id=v.var where v.filename='%s' and \
    %s group by v.var) a join hsk_mnemonics d on a.var = d.id \
    join hsk_values b on b.var=a.var and b.time_val=a.time_val;"
        % (hsk_file, str_limit_check)
    )
    c.execute(sql)
    d = c.fetchall()

    # set the alert flag if a telemetry has yellow/red alert
    sql = (
        "update hsk_values set v.alert_type = %d from hsk_values v join hsk_limit l on l.id=v.var where \
    v.filename = '%s' and %s "
        % (alert_types[str_alert], hsk_file, str_limit_check)
    )
    c.execute(sql)
    conn.commit()

    return d


def print_alerts(str_alert, hkt_alerts, bVerbose, output, conn, c):
    # str_alert = 'red'/'yellow'
    if len(hkt_alerts[str_alert]) > 0:
        tmpstr = '====== %s alerts found in %s ======' % (
            str_alert.upper(),
            ','.join(list(hkt_alerts[str_alert].keys())),
        )
        print(tmpstr)
        if output:
            output.write(tmpstr + "\n")
        ### if bVerbose: ###
        for ihkt in hkt_alerts[str_alert].keys():
            tmpstr = (
                '{:50}'.format('var')
                + '{:20}'.format('value')
                + '{:30}'.format('time [first alert]')
                + '{:30}'.format('time [last alert]')
            )
            print(tmpstr)
            if output:
                output.write(tmpstr + "\n")

            tmpstr = '-' * 140
            print(tmpstr)
            if output:
                output.write(tmpstr + "\n")

            for ivar in hkt_alerts[str_alert][ihkt]:
                tmpstr = (
                    '{:50}'.format(ivar[0])
                    + '{:20}'.format(ivar[1])
                    + ('%s' % ivar[2]).ljust(30)
                    + '%s' % ivar[3]
                )
                print(tmpstr)
                if output:
                    output.write(tmpstr + "\n")

                # log red/yellow alert record into hsk_alerts table, start date and end date, assuming HSK file doesn't span over 2 days
                try:
                    sql = (
                        "insert into hsk_alerts (date_val,alert_type,id) values ('%s',%d,'%s')"
                        % (ivar[2], alert_types[str_alert], ivar[0])
                    )
                    c.execute(sql)
                    conn.commit()
                except:
                    pass
                try:
                    sql = (
                        "insert into hsk_alerts (date_val,alert_type,id) values ('%s',%d,'%s')"
                        % (ivar[3], alert_types[str_alert], ivar[0])
                    )
                    c.execute(sql)
                    conn.commit()
                except:
                    pass

            print("\n")
            if output:
                output.write("[Check the telemetry plots at:]\n")
                output.write(
                    "https://oceandata.sci.gsfc.nasa.gov/tlm/plots?strDate=%s&hkt=%s\n"
                    % (ivar[2].strftime('%Y%m%d'), ivar[0])
                )
                output.write("\n")
    else:
        print('No %s alerts found.' % str_alert.upper())


def pace_hkt_check_limits(args):
    print('\nChecking...')
    status = 0
    sbcon = ''

    try:
        # parse DB access config file
        with open(args.db_config, 'r') as f:
            lines = f.read().splitlines()
            sbcon = ";".join(lines)
        if sbcon.lower().find('driver') < 0:
            sbcon += (
                ";DRIVER="
                + os.environ["SYBASE"]
                + "/DataAccess64/ODBC/lib/libsybdrvodb-sqllen4.so.fbo"
            )
        if sbcon.lower().find("port") < 0:
            sbcon += ";PORT=5000"
        conn = pyodbc.connect(sbcon, autocommit=True)
        c = conn.cursor()
    except:
        print("Not able to set up DB connection. Check --db_config settings.")
        return 1

    sql = "use pace_tlm"
    c.execute(sql)

    hkt_alerts = {'red': {}, 'yellow': {}}

    # loop through HKT files and check each telemetry in that file
    # against limits
    for hsk_file in args.hsk_files:
        sql = "select * from hsk_files where filename='%s'" % hsk_file
        c.execute(sql)
        try:
            d = c.fetchall()[0]
            print("%s: %s through %s" % (d[0], d[1], d[2]))

            # initiate hsk_alerts table with 0=green status for that day HSK file covers
            try:
                sql = (
                    "insert into hsk_alerts (date_val,alert_type,id) values ('%s',0,'checked')"
                    % d[1]
                )
                c.execute(sql)
                conn.commit()
            except:
                pass
            try:
                sql = (
                    "insert into hsk_alerts (date_val,alert_type,id) values ('%s',0,'checked')"
                    % d[1]
                )
                c.execute(sql)
                conn.commit()
            except:
                pass

        except Exception as e:
            print(e)
            print(hsk_file)

        # check against HKT limits, retrieve all mnemonics that have red/yellow alarms
        red_alerts = find_out_of_limits(hsk_file, conn, c, 'red')
        # print(red_alerts)
        if len(red_alerts) > 0:
            hkt_alerts['red'][hsk_file] = red_alerts

        yellow_alerts = find_out_of_limits(hsk_file, conn, c, 'yellow')
        # print(yellow_alerts)
        if len(yellow_alerts) > 0:
            hkt_alerts['yellow'][hsk_file] = yellow_alerts

    print('\nPACE HKT limit check done.')

    output = None
    if len(hkt_alerts['red']) > 0 or len(hkt_alerts['yellow']) > 0:
        status = 100
        if args.output:
            # if args.verbose:
            print("Writing output file: %s\n" % args.output)
            output = open(args.output, 'w')

        if args.sms:
            sms = open(args.sms, 'w')
            tmpstr = datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '  '
            if len(hkt_alerts['red'].keys()) > 0:
                tmpstr += "Red alerts in %s. " % (','.join(hkt_alerts['red'].keys()))
            if len(hkt_alerts['yellow'].keys()) > 0:
                tmpstr += "Yellow alerts in %s. " % (
                    ','.join(hkt_alerts['yellow'].keys())
                )
            sms.write(tmpstr)
            sms.close()

    # output alert info
    print_alerts('red', hkt_alerts, args.verbose, output, conn, c)
    print_alerts('yellow', hkt_alerts, args.verbose, output, conn, c)

    if output:
        output.close()

    c.close()
    conn.close()

    return status


if __name__ == '__main__':
    print("pace_hkt_check_limits", __version__)
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Checks telemetry in ingested HKT files, sends out alerts if values are over limits',
        epilog="""
    EXIT Status:
    0   : All green
    1   : Dunno, something horrible occurred
    100 : There are telemetry alarms    
    """,
    )
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--hsk_files', nargs='+', help='list of hkt file names')
    parser.add_argument('--db_config', type=str, help='path to the db config file')
    parser.add_argument(
        '--output', type=str, help='path to the optional output text file'
    )
    parser.add_argument('--sms', type=str, help='path to the optional SMS text')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    if not args.hsk_files:
        print("Please specify HKT files to be checked.")
        sys.exit(1)

    status = pace_hkt_check_limits(args)
    if status:
        sys.exit(status)