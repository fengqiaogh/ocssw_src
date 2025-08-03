#!/usr/bin/env python3

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
from datetime import datetime
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter

DEFAULT_BASE_URL = "https://oceandata.sci.gsfc.nasa.gov/manifest/tags"
MANIFEST_BASENAME = "manifest.json"


#  ------------------ DANGER -------------------
#
# The next 5 functions:
#    getSession
#    isRequestAuthFailure
#    httpdl
#    uncompressFile
#    get_file_time
#
# exist in two places:
#    OCSSWROOT/src/manifest/manifest.py
#    OCSSWROOT/src/scripts/seadasutils/ProcUtils.py
#
# Make sure changes get into both files.
#

DEFAULT_CHUNK_SIZE = 131072

# requests session object used to keep connections around
obpgSession = None

def getSession(verbose=0, ntries=5):
    global obpgSession

    if not obpgSession:
        # turn on debug statements for requests
        if verbose > 1:
            logging.basicConfig(level=logging.DEBUG)

        obpgSession = requests.Session()
        obpgSession.mount('https://', HTTPAdapter(max_retries=ntries))

        if verbose:
            print("OBPG session started")
    else:
        if verbose > 1:
            print("reusing existing OBPG session")

    return obpgSession

#  ------------------ DANGER -------------------
# See comment above
def isRequestAuthFailure(req) :
    ctype = req.headers.get('Content-Type')
    if ctype and ctype.startswith('text/html'):
        if "<title>Earthdata Login</title>" in req.text:
            return True
    return False

#  ------------------ DANGER -------------------
# See comment above
def httpdl(server, request, localpath='.', outputfilename=None, ntries=5,
           uncompress=False, timeout=30., verbose=0, force_download=False,
           chunk_size=DEFAULT_CHUNK_SIZE, timestamp=False):

    status = 0
    urlStr = 'https://' + server + request

    global obpgSession
    localpath = Path(localpath)
    getSession(verbose=verbose, ntries=ntries)

    modified_since = None
    headers = {}

    if not force_download:
        if outputfilename:
            ofile = localpath / outputfilename
            modified_since = get_file_time(ofile)
        else:
            rpath = Path(request.rstrip())
            if 'requested_files' in request:
                rpath = Path(request.rstrip().split('?')[0])
            ofile = localpath / rpath.name
            if re.search(r'(?<=\?)(\w+)', ofile.name):
                ofile = Path(ofile.name.split('?')[0])

            modified_since = get_file_time(ofile)

        if modified_since:
            headers = {"If-Modified-Since":modified_since.strftime("%a, %d %b %Y %H:%M:%S GMT")}

    with obpgSession.get(urlStr, stream=True, timeout=timeout, headers=headers) as req:

        if req.status_code != 200:
            status = req.status_code
        elif isRequestAuthFailure(req):
            status = 401
        else:
            if not Path.exists(localpath):
                os.umask(0o02)
                Path.mkdir(localpath, mode=0o2775, parents=True)

            if not outputfilename:
                cd = req.headers.get('Content-Disposition')
                if cd:
                    outputfilename = re.findall("filename=(.+)", cd)[0]
                else:
                    outputfilename = urlStr.split('/')[-1]

            ofile = localpath / outputfilename

            # This is here just in case we didn't get a 304 when we should have...
            download = True
            if 'last-modified' in req.headers:
                remote_lmt = req.headers['last-modified']
                remote_ftime = datetime.strptime(remote_lmt, "%a, %d %b %Y %H:%M:%S GMT").replace(tzinfo=None)
                if modified_since and not force_download:
                    if (remote_ftime - modified_since).total_seconds() < 0:
                        download = False
                        if verbose:
                            print("Skipping download of %s" % outputfilename)

            if download:
                total_length = req.headers.get('Content-Length')
                if total_length is None:
                    total_length = len(req.content)
                length_downloaded = 0
                total_length = int(total_length)
                if verbose >0:
                    print("Downloading %s (%8.2f MBs)" % (outputfilename,total_length /1024/1024))

                with open(ofile, 'wb') as fd:

                    for chunk in req.iter_content(chunk_size=chunk_size):
                        if chunk: # filter out keep-alive new chunks
                            length_downloaded += len(chunk)
                            fd.write(chunk)
                            if verbose > 0:
                                percent_done = int(50 * length_downloaded / total_length)
                                sys.stdout.write("\r[%s%s]" % ('=' * percent_done, ' ' * (50-percent_done)))
                                sys.stdout.flush()

                if uncompress:
                    if ofile.suffix in {'.Z', '.gz', '.bz2'}:
                        if verbose:
                            print("\nUncompressing {}".format(ofile))
                        compressStatus = uncompressFile(ofile)
                        if compressStatus:
                            status = compressStatus
                else:
                    status = 0

                if timestamp:
                    if 'Last-Modified' in req.headers:
                        last_modified = req.headers["Last-Modified"]
                        last_modified_time = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
                        last_modified_timestamp = time.mktime(last_modified_time.timetuple())
                        file_path = os.path.join(localpath, os.path.basename(request))
                        if os.path.exists(file_path):
                            os.utime(file_path, (last_modified_timestamp, last_modified_timestamp))

    return status


#  ------------------ DANGER -------------------
# See comment above
def uncompressFile(compressed_file):
    """
    uncompress file
    compression methods:
        bzip2
        gzip
        UNIX compress
    """

    compProg = {".gz": "gunzip", ".Z": "gunzip", ".bz2": "bunzip2"}
    exten = Path(compressed_file).suffix
    unzip = compProg[exten]
    cmd = [unzip, "-f", str(compressed_file.resolve())]
    p = subprocess.Popen(cmd, shell=False)
    status = os.waitpid(p.pid, 0)[1]
    if status:
        print("Warning! Unable to decompress %s" % compressed_file)
        return status
    else:
        return 0

#  ------------------ DANGER -------------------
# See comment above
def get_file_time(localFile):
    ftime = None
    localFile = Path(localFile)
    if not Path.is_file(localFile):
        while localFile.suffix in {'.Z', '.gz', '.bz2'}:
            localFile = localFile.with_suffix('')

    if Path.is_file(localFile):
        ftime = datetime.fromtimestamp(localFile.stat().st_mtime)

    return ftime

#  ------------------ Start of non-duplicated code -------------------

def run():
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=download)
    subparsers = parser.add_subparsers()

    _add_subparser_reprint(subparsers)
    _add_subparser_update_file(subparsers)
    _add_subparser_add_tag(subparsers)
    _add_subparser_get_value(subparsers)
    _add_subparser_get_first_tag(subparsers)
    _add_subparser_list(subparsers)
    _add_subparser_clean(subparsers)
    _add_subparser_download(subparsers)
    _add_subparser_generate(subparsers)
    _add_subparser_list_tags(subparsers)

    options, args = parser.parse_known_args()
    return options.func(options, args)

def _add_subparser_reprint(subparsers):
    parser_reprint = subparsers.add_parser('reprint')
    parser_reprint.add_argument("manifest", help="manifest to reprint")
    parser_reprint.set_defaults(func=reprint)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_reprint.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_update_file(subparsers):
    parser_update_file = subparsers.add_parser('update-file')
    parser_update_file.add_argument("manifest", help="manifest to update")
    parser_update_file.add_argument("path", help="file to update")
    parser_update_file.set_defaults(func=update_file)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_update_file.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_add_tag(subparsers):
    parser_add_tag = subparsers.add_parser('add-tag')
    parser_add_tag.add_argument("-m", "--manifest", help="manifest to update")
    parser_add_tag.add_argument("tag", help="tag to add to tags attribute")
    parser_add_tag.set_defaults(func=add_tag)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_add_tag.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_get_value(subparsers):
    parser_get_value = subparsers.add_parser('get-value')
    parser_get_value.add_argument("-m", "--manifest", help="manifest from which to retrieve the value")
    parser_get_value.add_argument("xpath", help="key to print, colon separated for nested values")
    parser_get_value.set_defaults(func=get_value)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_get_value.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_get_first_tag(subparsers):
    parser_get_first_tag = subparsers.add_parser('get-first-tag')
    parser_get_first_tag.add_argument("-m", "--manifest", help="manifest from which to retrieve the first tag")
    parser_get_first_tag.set_defaults(func=get_first_tag)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_get_first_tag.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_list(subparsers):
    parser_list = subparsers.add_parser('list')
    parser_list.add_argument("manifest", help="manifest to list")
    parser_list.add_argument("-i", "--info", action="store_const", const=1, help="include extra info")
    parser_list.add_argument("-t", "--tag", help="tag to list files for")
    parser_list.set_defaults(func=list)
    if os.path.isfile(MANIFEST_BASENAME):
        parser_list.set_defaults(manifest=MANIFEST_BASENAME)

def _add_subparser_clean(subparsers):
    parser_clean = subparsers.add_parser('clean')
    parser_clean.add_argument("-d", "--dry-run", action="store_const", const=1, help="don't actually delete files")
    parser_clean.add_argument("directory", default=".", nargs='?', help="directory to clean (must contain %s)" % MANIFEST_BASENAME)
    parser_clean.add_argument("-e", "--exclude", nargs="+", action='append', help="relative paths to ignore")
    parser_clean.add_argument("-i", "--include", nargs="+", action='append', help="relative paths to include (ignore *)")
    parser_clean.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser_clean.set_defaults(func=clean)

def _add_subparser_download(subparsers):
    parser_download = subparsers.add_parser('download')
    parser_download.add_argument("-d", "--dest-dir", help="destination directory")
    parser_download.add_argument("-t", "--tag", help="tag to download")
    parser_download.add_argument("-b", "--base-url", default=DEFAULT_BASE_URL, help="base URL")
    parser_download.add_argument("-n", "--name", help="bundle name")
    parser_download.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE, help="download chunk size")
    parser_download.add_argument("-s", "--save-dir", help="save a copy of the manifest files to this directory")
    parser_download.add_argument("-l", "--local-dir", help="directory containing local manifest files")
    parser_download.add_argument("--timestamp", default=False, action="store_true", help="preserve timestamps of downloaded files")
    parser_download.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser_download.add_argument("files", action="append", nargs="*", default=None, type=str, help="files to download if needed")

    parser_download.set_defaults(func=download)
    parser_download.set_defaults(dest_dir=".")

def _add_subparser_generate(subparsers):
    parser_gen = subparsers.add_parser('generate')
    parser_gen.add_argument("-b", "--base-manifest", help="base manifest file")
    parser_gen.add_argument("-c", "--checksum-bytes", default=1000000, help="how many bytes to checksum per file")
    parser_gen.add_argument("-t", "--tag", required=True, help="new tag for manifest")
    parser_gen.add_argument("-f", "--force", action="store_const", const=1, help="generate manifest despite warnings")
    parser_gen.add_argument("-e", "--exclude", nargs="+", action='append', help="relative paths to ignore")
    parser_gen.add_argument("-i", "--include", nargs="+", action='append', help="relative paths to include (ignore *)")
    parser_gen.add_argument("-n", "--name", help="bundle name")
    parser_gen.add_argument("directory", help="directory to generate a manifest for")
    parser_gen.set_defaults(func=generate)

def _add_subparser_list_tags(subparsers):
    parser_list_tags = subparsers.add_parser('list_tags')
    parser_list_tags.add_argument("-b", "--base-url", default=DEFAULT_BASE_URL, help="base URL")
    parser_list_tags.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE, help="download chunk size")
    parser_list_tags.add_argument("--timestamp", default=False, action="store_true", help="preserve timestamps of downloaded files")
    parser_list_tags.add_argument("-v", "--verbose", action="count", default=0, help="increase output verbosity")
    parser_list_tags.set_defaults(func=list_tags)

def create_default_options():
    options = argparse.Namespace(
                verbose=0,
                dest_dir=None,
                tag=None,
                base_url=DEFAULT_BASE_URL,
                name=None,
                chunk_size=DEFAULT_CHUNK_SIZE,
                save_dir=None,
                local_dir=None,
                timestamp=False,
                files=None,
                func=list_tags)
    return options

def run_command(command):
    proc = subprocess.run(command, shell=True)
    if proc.returncode != 0:
        print("Error: return =", proc.returncode, ": trying to run command =", command)
        sys.exit(1)

def reprint(options, args):
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        print(json.dumps(manifest, indent=4, sort_keys=True))

def update_file(options, args):
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        current_entry = manifest['files'].get(options.path)
        if os.path.islink(options.path):
            linkValue = os.readlink(options.path)
            if not current_entry or current_entry.get("symlink") != linkValue:
                info = {"symlink": linkValue, "tag": options.tag}
                manifest['files'][options.path] = info
        else:
            checksum = _get_checksum(manifest, options.path)
            if not current_entry or current_entry.get('checksum') != checksum:
                info = {
                    "checksum": checksum,
                    "size": os.stat(options.path).st_size,
                    "mode": os.stat(options.path).st_mode,
                    "tag": manifest['tag']
                }
                manifest['files'][options.path] = info

        print(json.dumps(manifest, indent=4, sort_keys=True))

def add_tag(options, args):
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        if options.tag not in manifest["tags"]:
            manifest["tags"].append(options.tag)
        else:
            print("%s is already in the tags attribute" % (options.tag), file=sys.stderr)

        print(json.dumps(manifest, indent=4, sort_keys=True))

def get_value(options, args):
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        for part in options.xpath.split(":"):
            if part in manifest:
                manifest = manifest[part]
            else:
                print("Path not found, invalid part: %s" % part)
                return
        print(manifest)

def get_first_tag(options, args):
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        print(manifest['tags'][0])

def getFileList(excludeList=None, includeList=None):
    allFiles = []

    for root, _, files in os.walk(".", followlinks=True):
        for f in files:
            if '/' in root:
                name = root[2:]+'/'+f
            else:
                name = f

            # exclude files if not in include list
            addIt = True
            if excludeList:
                for exclude in excludeList:
                    if exclude[0] == "." or name.startswith(exclude[0]):
                        addIt = False
                        if includeList:
                            for include in includeList:
                                if name.startswith(include[0]):
                                    addIt = True
                                    break
                        if not addIt:
                            break
            if addIt:
                if "__pycache__" not in name:
                    allFiles.append(name)

    return allFiles

def clean(options, args):
    os.chdir(options.directory)

    # check for exclude wild card
    if options.exclude:
        for exclude in options.exclude:
            if exclude[0] == ".":
                return

    if not os.path.isfile(MANIFEST_BASENAME):
        print("directory needs to contain a", MANIFEST_BASENAME)
        return 1

    with open(MANIFEST_BASENAME, 'rb') as manifest:
        manifest = json.load(manifest)
        files = manifest["files"]
        for f in getFileList(options.exclude, options.include):
            if f == MANIFEST_BASENAME:
                continue
            if not files.get(f):
                if options.verbose or options.dry_run:
                    print("cleaning %s" % (f))
                if not options.dry_run:
                    try:
                        os.remove(f)
                    except FileNotFoundError:
                        pass # just ignore the file if it does not exist

def list(options, args):
    if os.path.isdir(options.manifest):
        options.manifest = "%s/%s" % (options.manifest, MANIFEST_BASENAME)
    with open(options.manifest, 'rb') as manifest:
        manifest = json.load(manifest)
        if options.info:
            for f, info in manifest["files"].items():
                if not options.tag or info["tag"] == options.tag:
                    if info.get('symlink'):
                        print("%s %s, -> %s" % (f, info["tag"], info["symlink"]))
                    else:
                        print("%s %s, %s bytes, %s" % (f, info["tag"], info["size"], info["checksum"]))
        elif options.tag:
            for f, info in manifest["files"].items():
                if info["tag"] == options.tag:
                    print(f)
        else:
            for f in manifest["files"]:
                print(f)

def generate(options, args):
    if not options.base_manifest and os.path.isfile("%s/%s" % (options.directory, MANIFEST_BASENAME)):
        options.base_manifest = "%s/%s" % (options.directory, MANIFEST_BASENAME)

    manifest = None
    if options.base_manifest and os.path.isfile(options.base_manifest) and os.path.getsize(options.base_manifest):
        with open(options.base_manifest, 'rb') as base_manifest:
            manifest = json.load(base_manifest)
    else:
        manifest = {"checksum_bytes": options.checksum_bytes, "tags": []}

    manifest["tags"] = [options.tag]

    os.chdir(options.directory)

    all_files = getFileList(options.exclude, options.include)

    if options.name:
        manifest['name'] = options.name

    files_entries = manifest.get("files", {})

    # delete entries not in the directory
    files_to_delete = []
    if "files" in manifest:
        for path, info in manifest["files"].items():
            if path not in all_files:
                files_to_delete.append(path)
    for path in files_to_delete:
        del files_entries[path]

    for f in all_files:
        if os.path.basename(f) == MANIFEST_BASENAME:
            continue

        current_entry = files_entries.get(f)
        if os.path.islink(f):
            linkValue = os.readlink(f)
            if not current_entry or current_entry.get("symlink") != linkValue:
                info = {"symlink": linkValue, "tag": options.tag}
                files_entries[f] = info
        else:
            fileSize = os.path.getsize(f)
            checksum = _get_checksum(manifest, f)
            if not current_entry or current_entry.get('size') != fileSize or current_entry.get('checksum') != checksum:
                info = {
                    "checksum": checksum,
                    "size": fileSize,
                    "mode": os.stat(f).st_mode,
                    "tag": options.tag
                }
                files_entries[f] = info
    manifest["files"] = files_entries
    print(json.dumps(manifest, indent=4, sort_keys=True))

def download(options, args):
    manifest = None
    manifest_filename = "%s/%s" % (options.dest_dir, MANIFEST_BASENAME)

    if not os.path.isdir(options.dest_dir):
        os.makedirs(options.dest_dir)

    if options.local_dir:
        if options.save_dir:
            print("Error: Can not have --local_dir and --save_dir")
            return 1

    if not options.tag or not options.name:
        if not os.path.isfile(manifest_filename):
            print("must have -t and -n or %s" % (manifest_filename))
            return 1
        with open(manifest_filename, 'rb') as manifest:
            manifest = json.load(manifest)
        if not options.tag:
            options.tag = manifest['tags'][-1]
        if not options.name:
            options.name = manifest['name']

    if not _download_file(options, MANIFEST_BASENAME):
        return 1

    with open(manifest_filename, 'rb') as manifest:
        manifest = json.load(manifest)

    modified_files = _check_directory_against_manifest(options, options.dest_dir, manifest)

    # if files on command line only look at those
    if options.files and options.files[0]:
        newList = {}
        for f in options.files[0]:
            try:
                newList[f] = modified_files[f]
            except:
                pass
        modified_files = newList

    if not modified_files:
        if options.verbose:
            print("No files require downloading")
    else:
        _download_files(options, modified_files)

    if options.save_dir:
        for path, info in manifest['files'].items():
            if info.get('checksum'):
                src = "%s/%s" % (options.dest_dir, path)
                dest = "%s/%s/%s/%s" % (options.save_dir, info["tag"], options.name, path)
                destDir = os.path.dirname(dest)
                if not os.path.isdir(destDir):
                    os.makedirs(destDir)
                shutil.copy2(src, dest)
                os.chmod(dest, info["mode"])

def get_tags(options, args):
    tag_list = []
    tempDir = tempfile.TemporaryDirectory(prefix="manifest-")
    status = 0
    url = options.base_url + "/"
    parts = urllib.parse.urlparse(url)
    host = parts.netloc
    request = parts.path
    status = httpdl(host, request, localpath=tempDir.name,
                    outputfilename="index.html",
                    verbose=options.verbose,
                    force_download=True,
                    chunk_size=options.chunk_size,
                    timestamp=options.timestamp)
    if status == 0:
        with open("%s/index.html" % (tempDir.name)) as f:
            inBody = False
            for line in f:
                if "<body>" in line:
                    inBody = True
                if "</body>" in line:
                    break
                if inBody:
                    if line.startswith("<a href="):
                        parts = line.split('"')
                        s = parts[1].split("/")[0]
                        if s != "..":
                            tag_list.append(s)
        return tag_list
    else:
        print("Error downloading list of tags : return code =", status)
        return tag_list

def list_tags(options, args):
    for tag in get_tags(options, args):
        print(tag)

def check_tag(options, args):
    for tag in get_tags(options, args):
        if tag == options.tag:
            return True
    return False

def _get_checksum(manifest, path):
    checksum = hashlib.sha256()
    with open(path, 'rb') as current_file:
        checksum.update(current_file.read(manifest['checksum_bytes']))
    return checksum.hexdigest()

def _check_directory_against_manifest(options, directory, manifest):
    modified_files = {}
    for path, info in manifest['files'].items():
        dest = os.path.join(directory, path)
        if os.path.islink(dest):
            if info.get('symlink') != os.readlink(dest):
                modified_files[path] = info
        elif os.path.isfile(dest):
            if info.get('size') != os.path.getsize(dest) or info.get('checksum') != _get_checksum(manifest, dest) or info.get('mode') != os.stat(dest).st_mode:
                modified_files[path] = info
        else:
            modified_files[path] = info
    return modified_files

def _download_file(options, fileName):
    dest = "%s/%s" % (options.dest_dir, fileName)
    dest_dir = os.path.dirname(dest)
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)

    if options.local_dir:
        src = "%s/%s/%s/%s" % (options.local_dir, options.tag, options.name, fileName)
        if options.verbose:
            print("Copying %s from %s" % (fileName, src))
        if options.timestamp:
            shutil2.copy(src, dest)
        else:
            shutil.copy(src, dest)
        return True

    url = "%s/%s/%s/%s" % (options.base_url, options.tag, options.name, fileName)
    if options.verbose:
        print("Downloading %s from %s" % (fileName, url))
    parts = urllib.parse.urlparse(url)
    #host = "%s://%s" % (parts.scheme, parts.netloc)
    host = parts.netloc
    request = parts.path
    status = httpdl(host, request, localpath=dest_dir,
                    outputfilename=os.path.basename(dest),
                    verbose=options.verbose,
                    force_download=True,
                    chunk_size=options.chunk_size,
                    timestamp=options.timestamp)
    if status != 0:
        print("Error downloading", dest, ": return code =", status)
        return False

    if options.save_dir:
        src = "%s/%s" % (options.dest_dir, fileName)
        dest = "%s/%s/%s/%s" % (options.save_dir, options.tag, options.name, fileName)
        destDir = os.path.dirname(dest)
        if not os.path.isdir(destDir):
            os.makedirs(destDir)
        shutil.copy2(src, dest)
    return True

def _download_files(options, file_list):
    if options.local_dir:
        for path, info in file_list.items():
            dest = "%s/%s" % (options.dest_dir, path)
            dest_dir = os.path.dirname(dest)
            if not os.path.isdir(dest_dir):
                os.makedirs(dest_dir)
            if info.get('checksum'):
                src = "%s/%s/%s/%s" % (options.local_dir, info["tag"], options.name, path)
                shutil.copy2(src, dest)
                os.chmod(dest, info["mode"])
            else:
                src = info['symlink']
                os.symlink(src, dest)
        return

    for path, info in file_list.items():
        dest = "%s/%s" % (options.dest_dir, path)
        dest_dir = os.path.dirname(dest)
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)

        if info.get('checksum'):
            if os.path.islink(dest) or os.path.exists(dest):
                os.remove(dest)
            url = "%s/%s/%s/%s" % (options.base_url, info["tag"], options.name, path)
            if options.verbose:
                print("Downloading %s from %s" % (path, url))
            parts = urllib.parse.urlparse(url)
            host = parts.netloc
            request = parts.path
            status = httpdl(host, request, localpath=dest_dir,
                        outputfilename=os.path.basename(dest),
                        verbose=options.verbose,
                        force_download=True,
                        chunk_size=options.chunk_size,
                        timestamp=options.timestamp)
            if status == 0:
                os.chmod(dest, info["mode"])
            else:
                print("Error downloading", dest, ": return code =", status)
        else:
            src = info['symlink']
            if options.verbose:
                print("Making symlink %s -> %s" % (dest, src))
            if os.path.islink(dest) or os.path.exists(dest):
                os.remove(dest)
            os.symlink(src, dest)


if __name__ == "__main__":
    sys.exit(run())
