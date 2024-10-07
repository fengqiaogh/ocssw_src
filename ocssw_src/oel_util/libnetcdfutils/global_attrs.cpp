#include "global_attrs.h"
#include <ctime>
#include <netcdf>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include <cstdio>
#include <iostream>
#include <libgen.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace rapidjson;

string call_sequence(int argc, char* argv[]) {

  // get processing timestamp, to the second
  time_t tnow = time(nullptr);
  char prodtime[80];
  strftime(prodtime, 80, "%Y-%m-%dT%XZ", gmtime(&tnow));

  // append calling sequence
  string callseq = prodtime;
  callseq.append(":");
  for (int i=0; i<argc; i++) {
    callseq.append(" ");
    callseq.append(basename(argv[i]));
  }
  return callseq;
}

string get_history(NcFile *ncfile) {
  string history;
  NcGroupAtt att = ncfile->getAtt("history");
  if (!att.isNull()) att.getValues(history);
  if (history.length() > 0) { history.append("; "); }
  return history;
}

void set_global_attrs(NcFile *outfile,
                      string history,
                      string doi,
                      string pversion) {

  // read standard global attributes from JSON file
  char *dataRoot;
  if ((dataRoot = getenv("OCDATAROOT")) == NULL) {
    printf("-E- OCDATAROOT environment variable is not defined.\n");
    exit(EXIT_FAILURE);
  }
  string filename = (string) dataRoot + "/common/global_metadata.json";
  char readBuffer[65536];
  Document document;

  try {
    FILE* fp = fopen((char*) filename.c_str(), "r");
    if (fp == nullptr) {
      cerr << "Error: reading " << filename << endl;
      exit(EXIT_FAILURE);
    }
    FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    fclose(fp);
    document.ParseStream(is);
  }
  catch (std::exception const & e) {
    cerr << "Exception: " << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  // set standard global attributes
  for (Value::ConstMemberIterator itr = document.MemberBegin(); itr != document.MemberEnd(); ++itr)
    outfile->putAtt(itr->name.GetString(), itr->value.GetString());  // fix for non-string attrs

  // set pversion info
  if (!pversion.empty()) {
    outfile->putAtt("processing_version", pversion);
  }
  // set doi info
  if (!doi.empty()) {
    outfile->putAtt("identifier_product_doi_authority", "https://dx.doi.org");
    outfile->putAtt("identifier_product_doi", doi);
  }
  // set processing history
  if (!history.empty()) {
    outfile->putAtt("history", history);
  }

  // set date_created
  outfile->putAtt("date_created", unix2isodate(now(), 'G'));

}

void set_global_attrs(string filename,
                      string history,
                      string doi,
                      string pversion) {

  // open file for append
  NcFile *outfile = new NcFile(filename, NcFile::write);

  // set attributes
  set_global_attrs(outfile, history, doi, pversion);

  // close file
  outfile->close();
  delete(outfile);
}
