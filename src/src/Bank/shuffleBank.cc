// $Id$

// This program outputs reads/contigs from a bank 
// into specified directories

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <getopt.h>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <string>
#include <math.h>
#include <functional>
#include "foundation_AMOS.hh"
#include <fstream>
#include <ctype.h>
#include <sys/resource.h>

using namespace std;
using namespace AMOS;
using namespace HASHMAP;

// endl is a sham: die endl, die!
#define endl "\n"

const Size_t NUM_BANK_FILES = 10;
const string UNKNOWN_PREFIX = "UNKNOWN";

map<string, string> globals; // global variables

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

void printHelpText()
{
  cerr << "\n.USAGE.\n"
       << "  shuffleFasta -b <bank_name>\n"
       << "\n.DESCRIPTION.\n"
       << "  shuffleFasta - generates a .fasta (and .qual) or fastq file from the bank\n"
       << "\n.OPTIONS.\n"
       << "  -h, -help     print out help message\n"
       << "  -b <bank_name>, -bank     bank where assembly is stored\n"
       << "  -p <prefix>   Common file prefix to add to the output\n"
       << "  -c            Dump contigs from the bank (default)\n"
       << "  -r            Dump reads from the bank\n"
       << "  -eid          report eids\n"
       << "  -iid          report iids (default)\n"
       << "  -f            Dump in fastq format\n"
       << "  -Q <int>      Use this as the min base quality (default: 33 / Sanger FASTQ)\n"
       << "  -a            Ignore clear range and dump entire sequence\n"
       << "  -d            Display details on header line\n"
       << "  -L <num>      Set the maximum number of bases per line (Default: 70)\n"
       << "  -E	<fofn>      List of files specifying by EID where to write\n"
       << "  -I <fofn>     List of files specifying by EID where to write\n"
       << "\n.KEYWORDS.\n"
       << "  AMOS bank, Converters\n"
       << endl;
}

bool GetOptions(int argc, char ** argv)
{  
  int option_index = 0;
  static struct option long_options[] = {
    {"help",      0, 0, 'h'},
    {"h",         0, 0, 'h'},
    {"d",         0, 0, 'd'},
    {"f",         0, 0, 'f'},
    {"Q",         1, 0, 'Q'},
    {"a",         0, 0, 'a'},
    {"L",         1, 0, 'L'},
    {"c",         0, 0, 'c'},
    {"r",         0, 0, 'r'},
    {"bank",      1, 0, 'b'},
    {"eid",       0, 0, 'e'},
    {"iid",       0, 0, 'i'},
    {"E",         1, 0, 'E'},
    {"I",         1, 0, 'I'},
    {"p",         1, 0, 'p'},
    {0, 0, 0, 0}
  };
  
  int c;
  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index))!= -1){
    switch (c){
    case 'h':
      printHelpText();
      exit(0);
      break;
    case 'b':
      globals["bank"] = string(optarg);
      break;
    case 'p':
      globals["prefix"] = "." + string(optarg);
      break;
    case 'd':
      globals["details"] = "true";
      break;
    case 'f':
      globals["fastq"] = "true";
      break;
    case 'c':
      globals["type"] = "CTG";
      break;
    case 'r':
      globals["type"] = "RED";
      break;
    case 'a':
      globals["noclear"] = "true";
      break;
    case 'Q':
      globals["qualOffset"] = string(optarg);
      break;
    case 'L':
      globals["lineLength"] = string(optarg);
      break;
    case 'e':
      globals["eid"] = "true";
      globals["iid"] = "false";
      break;
    case 'i':
      globals["iid"] = "true";
      globals["eid"] = "false";
      break;
    case 'E':
      globals["eidfile"] = string(optarg);
      break;
    case 'I':
      globals["iidfile"] = string(optarg);
      break;
    case '?':
      return false;
    }
  }

  return true;
} // GetOptions

string getQualValue(char qualVal) {
   if (globals["fastq"] == "false") {
      return convertInt((int) (toascii(qualVal) - toascii(AMOS::MIN_QUALITY))) + " ";
   } else {
      char converted = (char) (toascii(qualVal) - toascii(AMOS::MIN_QUALITY) + atoi(globals["qualOffset"].c_str()));
      stringstream ss;
      ss << converted;
      return ss.str();
   }
}

void outputFastX(const string &seqHeader, const string &seq, const string &qualHeader, const string &quals, ofstream &fasta, ofstream &qual, int lineLen = 0) {
  int nout = 0;

  if (globals["fastq"] == "true") {
    fasta << "@";
  } else {
    fasta << ">";
  }
  fasta << seqHeader;
  for (int i = 0; i < seq.length(); i++)
  {
    if (seq[i] == '-') { continue; }
    if (lineLen > 0 && nout % lineLen == 0)
      fasta << endl;

    nout++;
    fasta << seq[i];
  }
  fasta << endl;
  
  nout = 0;
  if (qual.is_open()) {
     if (globals["fastq"] == "true") {
        qual << "+";
     } else {
        qual << ">";
     }
     qual << qualHeader;
     for (int i = 0; i < quals.length(); i++) {
        string qv = getQualValue(quals[i]);
        if (lineLen > 0 && nout % (lineLen/qv.length()) == 0 && qual.is_open())
           qual << endl;
        nout++;
        qual << qv;
     }
     qual << endl;
  }
}

void printRead(Read_t & red, ofstream & fasta, ofstream & qual)
{
  if ( red . getLength( ) <= 0 ) {
      cerr << "WARNING: read with IID " << red . getIID( )
           << " has no sequence, skipped\n";
      return;
  }
  if ( globals.find("noclear") != globals.end() && red . getClearRange( ) . getLength( ) <= 0 ) {
      cerr << "WARNING: read with IID " << red . getIID( )
           << " has no clear range sequence, skipped\n";
      return;
  }

  stringstream fastaHeader;
  stringstream qualHeader;
  string seq = (globals.find("noclear") == globals.end() ? red.getSeqString() : red.getSeqString(red.getClearRange()));
  string quals = (globals.find("noclear") == globals.end() ? red.getQualString() : red.getQualString(red.getClearRange()));

  if (globals["eid"] == "true") {
    fastaHeader << red.getEID();
    qualHeader << red.getEID();
  } else {
    fastaHeader << red.getIID();
    qualHeader << red.getIID();
  }

  if (globals.find("details") != globals.end()) {
     fastaHeader << "len=" << seq.length()
         << " clear_start=" << red.getClearRange().getLo() + 1
         << " clear_end=" << red.getClearRange().getHi();
  }
  outputFastX(fastaHeader.str(), seq, qualHeader.str(), quals, fasta, qual, atoi(globals["lineLength"].c_str()));
}

void printContig(const Contig_t & ctg, ofstream & fasta, ofstream & qual)
{
  string seq = ctg.getSeqString();
  string quals = ctg.getQualString();
  stringstream fastaHeader;
  stringstream qualHeader;

  if (globals["eid"] == "true") { 
    fastaHeader << ctg.getEID(); 
    qualHeader << ctg.getEID();
  } else { 
    fastaHeader << ctg.getIID(); 
    qualHeader << ctg.getIID();
  }

  if (globals.find("details") != globals.end())
  {
    int len = ctg.getUngappedLength();

    std::vector<Tile_t>::const_iterator ti;
    const std::vector<Tile_t> & tiling = ctg.getReadTiling();

    int nreads = ctg.getReadTiling().size();
    double covsum = 0.0;

    for (ti = tiling.begin();
         ti != tiling.end();
         ti++)
    {
      covsum += ti->range.getLength();
    }

    double cov = len ? covsum / len : 0.0;

    fastaHeader << " len=" << len
         << " nreads=" << nreads
         << " cov=" << cov;
  }

  outputFastX(fastaHeader.str(), seq, qualHeader.str(), quals, fasta, qual, atoi(globals["lineLength"].c_str()));
} // printFasta

//----------------------------------------------
int main(int argc, char **argv)
{
  ofstream outqual;
  globals["iid"] = "true";
  globals["eid"] = "false";
  globals["fastq"] = "false";
  globals["type"] = "CTG";
  globals["qualOffset"] = "33";
  globals["lineLength"] = "60";
  globals["prefix"] = "";

  if (! GetOptions(argc, argv)){
    cerr << "Command line parsing failed" << endl;
    printHelpText();
    exit(1);
  }

  // open necessary files
  if (globals.find("bank") == globals.end()){ // no bank was specified
    cerr << "A bank must be specified" << endl;
    exit(1);
  }

  Bank_t bank (globals["type"]);
  BankStream_t stream (globals["type"]);
  if (! bank.exists(globals["bank"])){
    cerr << "No account found in bank " << globals["bank"] << endl;
    exit(1);
  }

  try 
  {
    stream.open(globals["bank"], B_READ);
  } 
  catch (Exception_t & e)
  {
    cerr << "Failed to open account in bank " << globals["bank"] 
         << ": " << endl << e << endl;
    exit(1);
  }
  

  try
  {
    ifstream file;
    string id;
    map <string, string> idToFasta;
    map <string, ofstream *> prefixToFastaFile;
    map <string, ofstream *> prefixToQualFile;
    string unknownPrefix = "";

    if (!globals["eidfile"].empty())
    {
      file.open(globals["eidfile"].c_str());
      
      if (!file)
      {
        throw Exception_t("Couldn't open EID File", __LINE__, __FILE__);
      }

      while (file >> id)
      {
        ifstream currFofn;
        currFofn.open(id.c_str());
        if (!currFofn) {
           throw Exception_t("Couldn't open EID fofn File", __LINE__, __FILE__);
        }
        string prefix = id.substr(0, id.find_last_of(".")) + globals["prefix"];
        while (currFofn >> id) {
           if (!stream.existsEID(id)) {
              cerr << "WARNING: skipping " << id << " it does not exist in the bank!" << endl;
              continue;
           }
           idToFasta[id] = prefix;
           prefixToFastaFile[prefix] = new ofstream();
           prefixToQualFile[prefix] = new ofstream();
           if (prefix.find(UNKNOWN_PREFIX) != std::string::npos) {
              unknownPrefix = prefix;
           }
        }
        currFofn.close();
      }
      file.close();
    }
    else if (!globals["iidfile"].empty())
    {
      file.open(globals["iidfile"].c_str());

      if (!file)
      {
        throw Exception_t("Couldn't open IID File", __LINE__, __FILE__);
      }

      while (file >> id)
      {
        ifstream currFofn;
        currFofn.open(id.c_str());
        if (!currFofn) {
           throw Exception_t("Couldn't open EID fofn File", __LINE__, __FILE__);
        }
        string prefix = id.substr(0, id.find_last_of(".")) + globals["prefix"];

        while (currFofn >> id)
        {
           if (!stream.existsIID(atoi(id.c_str()))) {
              cerr << "WARNING: skipping " << id << " it does not exist in the bank!" << endl;
              continue;
           }
           idToFasta[id] = prefix;
           prefixToFastaFile[prefix] = new ofstream();
           prefixToQualFile[prefix] = new ofstream();
           if (prefix.find(UNKNOWN_PREFIX) != std::string::npos) {
              unknownPrefix = prefix;
           }
        }
        currFofn.close();
      }
      file.close();
    }

    IBankable_t *seq = NULL;
    Contig_t ctg;
    Read_t read;
    if (globals["type"] == "CTG") {
       seq = &ctg;
    } else if (globals["type"] == "RED") {
       seq = &read;
    }
    rlimit limits;
    getrlimit(RLIMIT_NOFILE, &limits);
    uint32_t maxFiles = limits.rlim_cur - NUM_BANK_FILES;
    if (globals["fastq"] == "false") {
       maxFiles /= 2;
    }
    queue<string> openFiles;
    map<string, bool> append;
    while (stream >> (*seq)) {
       string prefix;
       if (!globals["eidfile"].empty()) {
          if (idToFasta.find(seq->getEID()) != idToFasta.end()) {
             prefix = idToFasta[seq->getEID()];
          } else {
             cerr << "WARNING: " << seq->getEID() << " it has not been assigned, marking as UNKNOWN!" << endl;
             prefix = unknownPrefix;
          }
       } else if (!globals["iidfile"].empty()) {
          string id = convertInt(seq->getIID());
          if (idToFasta.find(id) != idToFasta.end()) {
             prefix = idToFasta[id];
          } else {
             cerr << "WARNING: " << seq->getEID() << " it has not been assigned, marking as UNKNOWN!" << endl;
             prefix = unknownPrefix;
          }
       }
       if (prefixToFastaFile.find(prefix) == prefixToFastaFile.end() || prefixToQualFile.find(prefix) == prefixToQualFile.end()) {
          cerr << "WARNING: no file defined for prefix " << prefix << endl;
          continue;
        }
        ofstream *fastaFile = prefixToFastaFile[prefix];
        ofstream *qualFile = prefixToQualFile[prefix];
        if (!fastaFile->is_open() || (globals["fastq"] == "false" && !qualFile->is_open())) {
           if (openFiles.size()+1 >= maxFiles) {
              string toClose = openFiles.front();
              openFiles.pop();
              if (!prefixToFastaFile[toClose]->is_open() || !prefixToQualFile[toClose]->is_open()) {
                 cerr << "ERROR: Expected that file " << toClose << " is already open but it is closed!";
                 exit(1);
              }
              prefixToFastaFile[toClose]->close();
              prefixToQualFile[toClose]->close();
              append[toClose] = true;
           }
           fastaFile->open(string(prefix + (globals["fastq"] == "true" ? ".fastq" : ".fasta")).c_str(), ofstream::out | (append.find(prefix) == append.end() ? ofstream::trunc : ofstream::app));
           if (!fastaFile->is_open())
           {
              throw Exception_t("Couldn't open EID File", __LINE__, __FILE__);
           }
           if (globals["fastq"] == "true") {
              prefixToQualFile[prefix] = fastaFile;
              qualFile = prefixToQualFile[prefix];
           } else {
              qualFile->open(string(prefix + ".qual").c_str(), ofstream::out | (append.find(prefix) == append.end() ? ofstream::trunc : ofstream::app));
              if (!qualFile->is_open())
              {
                 throw Exception_t("Couldn't open EID File", __LINE__, __FILE__);
              }
          }
          openFiles.push(prefix);
       }
       if (globals["type"] == "CTG") {
          printContig(*((Contig_t *)seq), (*fastaFile), (*qualFile));
       } else {
          printRead(*((Read_t *)seq), (*fastaFile), (*qualFile));
       }
    }

    for(std::map<string, ofstream*>::iterator iter = prefixToFastaFile.begin(); iter != prefixToFastaFile.end(); ++iter) {
       iter->second->close();
       delete iter->second;
    }
    if (globals["fastq"] == "false") {
       for(std::map<string, ofstream*>::iterator iter = prefixToQualFile.begin(); iter != prefixToQualFile.end(); ++iter) {
          iter->second->close();
          delete iter->second;
       }
    }

    stream.close();
  }
  catch (Exception_t & e)
  {
    cerr << "ERROR: -- Fatal AMOS Exception --\n" << e;
    return EXIT_FAILURE;
  }


  return(0);
} // main
