//
//  mapparser.cpp
//  express
//
//  Created by Adam Roberts on 8/19/11.
//  Copyright 2011 Adam Roberts. All rights reserved.
//

#include "mapparser.h"
#include "main.h"
#include "fragments.h"
#include "targets.h"
#include "threadsafety.h"
#include "library.h"
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

const size_t BUFF_SIZE = 9999;

/**
 * A helper functon that calculates the length of the reference spanned by the
 * read and populates the indel vectors (for SAM input).
 * @param cigar_str a pointer to the char array containing the cigar string.
 * @param inserts an empty Indel vector into which to add inserts.
 * @param deletes an empty Indel vector into which to add deletions.
 */
size_t cigar_length(const char* cigar_str, vector<Indel>& inserts,
                    vector<Indel>& deletes) {
  inserts.clear();
  deletes.clear();
  const char* p_cig = cigar_str;
  size_t i = 0; // read index
  size_t j = 0; // genomic index
  while (*p_cig) {
    char* t;
    size_t op_len = (size_t)strtol(p_cig, &t, 10);
    char op_char = toupper(*t);
    switch(op_char) {
      case 'I':
      case 'S':
        inserts.push_back(Indel(i, op_len));
        i += op_len;
        break;
      case 'D':
        deletes.push_back(Indel(i, op_len));
        j += op_len;
        break;
      case 'M':
      case 'N':
        i += op_len;
        j += op_len;
        break;
    }
    p_cig = t + 1;
  }
  return j;
}

/**
 * A helper functon that calculates the length of the reference spanned by the
 * read and populates the indel vectors (for BAM input).
 * @param cigar_vec a vector containing the split cigar string.
 * @param inserts an empty Indel vector into which to add inserts.
 * @param deletes an empty Indel vector into which to add deletions.
 */
size_t cigar_length(vector<BamTools::CigarOp>& cigar_vec,
                    vector<Indel>& inserts, vector<Indel>& deletes) {
  inserts.clear();
  deletes.clear();
  size_t i = 0; // read index
  size_t j = 0; // genomic index
  for (size_t k = 0; k < cigar_vec.size(); ++k) {
    char op_char = cigar_vec[k].Type;
    size_t op_len = cigar_vec[k].Length;
    switch(op_char) {
      case 'I':
      case 'S':
        inserts.push_back(Indel(i, op_len));
        i += op_len;
        break;
      case 'D':
        deletes.push_back(Indel(i, op_len));
        j += op_len;
        break;
      case 'M':
      case 'N':
        i += op_len;
        j += op_len;
        break;
    }
  }
  return j;
}

MapParser::MapParser(Library* lib, bool write_active)
    : _lib(lib), _write_active(write_active) {

  string in_file = lib->in_file_name;
  string out_file = lib->out_file_name;
  bool is_sam = false;
  if (in_file.size() == 0) {
    logger.info("No alignment file specified. Expecting streaming input on "
                "stdin...\n");
    _parser.reset(new SAMParser(&cin));
    is_sam = true;
  } else {
    logger.info_stderr("Attempting to read '%s' in BAM format...", in_file.c_str());
    BamTools::BamReader* reader = new BamTools::BamReader();
    if (reader->Open(in_file)) {
      logger.info_stderr("Parsing BAM header...");
      _parser.reset(new BAMParser(reader));
      if (out_file.size()) {
        out_file += ".bam";
        BamTools::BamWriter* writer = new BamTools::BamWriter();
        if (writer->Open(out_file, reader->GetHeader(),
                         reader->GetReferenceData())) {
          bool sample = out_file.substr(out_file.length()-8,4) == "samp";
          _writer.reset(new BAMWriter(writer, sample));
        } else {
          logger.severe("Unable to open output BAM file '%s'.",
                        out_file.c_str());
        }
      }
    } else {
      delete reader;
      logger.info_stderr("Input is not in BAM format. Trying SAM...");
      ifstream* ifs = new ifstream(in_file.c_str());
      if (!ifs->is_open()) {
        logger.severe("Unable to open input SAM file '%s'.", in_file.c_str());
      }
      _parser.reset(new SAMParser(ifs));
      is_sam = true;
    }
  }

  if (is_sam && out_file.size()) {
    out_file += ".sam";
    ofstream* ofs = new ofstream(out_file.c_str());
    if (!ofs->is_open()) {
      logger.severe("Unable to open output SAM file '%s'.", out_file.c_str());
    }
    *ofs << _parser->header();
    bool sample = out_file.substr(out_file.length()-8,4) == "samp";
    _writer.reset(new SAMWriter(ofs, sample));
  }
}

void MapParser::threaded_parse(ParseThreadSafety* thread_safety_p,
                               size_t stop_at,
                               size_t num_neighbors) {
  ParseThreadSafety& pts = *thread_safety_p;
  bool fragments_remain = true;
  size_t n = 0;
  size_t still_out = 0;

  TargetTable& targ_table = *(_lib->targ_table);

  while (!stop_at || n < stop_at) {
    Fragment* frag = NULL;
    while (fragments_remain) {
      frag = new Fragment(_lib);
      fragments_remain = _parser->next_fragment(*frag);
      if (frag->num_hits()) {
        break;
      }
      delete frag;
      frag = NULL;
    }
    for (size_t i = 0; frag && i < frag->hits().size(); ++i) {
      FragHit& m = *(frag->hits()[i]);
      
      if (m.first_read() && m.first_read()->seq.length() > max_read_len) {
        logger.severe("Length of first read for fragment '%s' is longer than "
                      "maximum allowed read length (%d vs. %d). Increase the "
                      "limit using the '--max-read-len,L' option.",
                      m.frag_name().c_str(),  m.first_read()->seq.length(),
                      max_read_len);
      }
      if (m.second_read() && m.second_read()->seq.length() > max_read_len) {
        logger.severe("Length of second read for fragment '%s' is longer than "
                      "maximum allowed read length (%d vs. %d). Increase the "
                      "limit using the '--max-read-len,L' option.",
                      m.frag_name().c_str(),  m.second_read()->seq.length(),
                      max_read_len);
      }
      
      Target* t = targ_table.get_targ(m.target_id());
      if (!t) {
        logger.severe("Target sequence at index '%s' not found. Verify that it "
                      "is in the SAM/BAM header and FASTA file.",
                      m.target_id());
      }
      m.target(t);
      assert(t->id() == m.target_id());
      
      // Add num_neighbors targets on either side to the neighbors list.
      // Used for experimental feature.
      vector<const Target*> neighbors;
      for (TargID j = 1; j <= num_neighbors;  j++) {
        if (j <= m.target_id()) {
          neighbors.push_back(targ_table.get_targ(m.target_id() - j));
        }
        if (j + m.target_id() < targ_table.size()) {
          neighbors.push_back(targ_table.get_targ(m.target_id() + j));
        }
      }
      m.neighbors(neighbors);
    }
    boost::scoped_ptr<Fragment> done_frag(pts.proc_out.pop(false));
    while (done_frag) {
      if (_writer && _write_active) {
        _writer->write_fragment(*done_frag);
      }
      still_out--;
      done_frag.reset(pts.proc_out.pop(false));
    }

    if (!frag) {
      break;
    }

    pts.proc_in.push(frag);
    n++;
    still_out++;
  }

  pts.proc_in.push(NULL);

  while (still_out) {
    boost::scoped_ptr<Fragment> done_frag(pts.proc_out.pop(true));
    if (_writer && _write_active) {
      _writer->write_fragment(*done_frag);
    }
    still_out--;
  }
}

BAMParser::BAMParser(BamTools::BamReader* reader) : _reader(reader) {
  BamTools::BamAlignment a;

  size_t index = 0;
  foreach(const BamTools::RefData& ref, _reader->GetReferenceData()) {
    if (_targ_index.count(ref.RefName)) {
      logger.severe("Target '%s' appears multiple times in BAM header.",
                    ref.RefName.c_str());
    }
    _targ_index[ref.RefName] = index++;
    _targ_lengths[ref.RefName] = ref.RefLength;
  }

  // Get first valid ReadHit
  _read_buff = new ReadHit();
  do {
    if (!_reader->GetNextAlignment(a)) {
      logger.severe("Input BAM file contains no valid alignments.");
    }
  } while(!map_end_from_alignment(a));
}

bool BAMParser::next_fragment(Fragment& nf) {
  nf.add_map_end(_read_buff);

  BamTools::BamAlignment a;
  _read_buff = new ReadHit();

  while(true) {
    if (!_reader->GetNextAlignment(a)) {
      return false;
    } else if (!map_end_from_alignment(a)) {
      continue;
    } else if (!nf.add_map_end(_read_buff)) {
      return true;
    }
    _read_buff = new ReadHit();
  }
}

bool BAMParser::map_end_from_alignment(BamTools::BamAlignment& a) {
  ReadHit& r = *_read_buff;

  if (!a.IsMapped()) {
    return false;
  }

  bool is_paired = a.IsPaired();

  if (is_paired && !a.IsProperPair()) {
    return false;
  }

  if (is_paired && (direction == F || direction == R)) {
    return false;
  }
  
  bool is_reversed = a.IsReverseStrand();
  bool is_mate_reversed = a.IsMateReverseStrand();
  if (is_paired && (!a.IsMateMapped() || a.RefID != a.MateRefID ||
                    is_reversed == is_mate_reversed ||
                    (is_reversed && a.MatePosition > a.Position) ||
                    (is_mate_reversed && a.MatePosition < a.Position))) {
    return false;
  }
  
  bool left_first = (!is_paired && !is_reversed) ||
                    (a.IsFirstMate() && !is_reversed)||
                    (a.IsSecondMate() && is_reversed);

  if (((direction == RF || direction == R) && left_first) ||
      ((direction == FR || direction == F) && !left_first)) {
    return false;
  }
    
  r.name = a.Name;
  if (boost::algorithm::ends_with(r.name, "\1") ||
      boost::algorithm::ends_with(r.name, "\2")) {
    r.name = r.name.substr(r.name.size()-2);
  }

  r.reversed = is_reversed;
  r.first = !is_paired || a.IsFirstMate();
  r.targ_id = a.RefID;
  r.left = a.Position;
  r.mate_l = a.MatePosition;
  r.seq.set(a.QueryBases, is_reversed);
  r.bam = a;
  r.right = r.left + cigar_length(a.CigarData, r.inserts, r.deletes);
  
  foreach (Indel& indel, r.inserts) {
    if (indel.len > max_indel_size) {
      return false;
    }
  }
  foreach (Indel& indel, r.deletes) {
    if (indel.len > max_indel_size) {
      return false;
    }
  }
  
  return true;
}

void BAMParser::reset() {
  _reader->Rewind();

  // Get first valid FragHit
  BamTools::BamAlignment a;
  delete _read_buff;
  _read_buff = new ReadHit();
  do {
    _reader->GetNextAlignment(a);
  } while(!map_end_from_alignment(a));
}

SAMParser::SAMParser(istream* in) {
  _in = in;

  char line_buff[BUFF_SIZE];
  _read_buff = new ReadHit();
  _header = "";

  // Parse header
  size_t index = 0;
  while(_in->good()) {
    _in->getline(line_buff, BUFF_SIZE-1, '\n');
    if (line_buff[0] != '@') {
      break;
    }
    _header += line_buff;
    _header += "\n";

    string str(line_buff);
    size_t idx = str.find("SN:");
    if (idx!=string::npos) {
      string name = str.substr(idx+3);
      name = name.substr(0,name.find_first_of("\n\t "));
      if (_targ_index.count(name)) {
        logger.severe("Target '%s' appears multiple times in SAM header.",
                      str.c_str());
      }
      _targ_index[name] = index++;
      idx = str.find("LN:");
      if (idx != string::npos) {
        string len = str.substr(idx+3);
        len = len.substr(0,len.find_first_of("\n\t "));
        _targ_lengths[name] = atoi(len.c_str());
      }
    }
  }

  // Load first aligned read
  while(!map_end_from_line(line_buff)) {
    if (!_in->good()) {
      logger.severe("Input SAM file contains no valid alignments.");
    }
    _in->getline(line_buff, BUFF_SIZE-1, '\n');
  }
}

bool SAMParser::next_fragment(Fragment& nf) {
  nf.add_map_end(_read_buff);

  _read_buff = new ReadHit();
  char line_buff[BUFF_SIZE];

  while(_in->good()) {
    _in->getline (line_buff, BUFF_SIZE-1, '\n');
    if (!map_end_from_line(line_buff)) {
      continue;
    }
    if (!nf.add_map_end(_read_buff)) {
      break;
    }
    _read_buff = new ReadHit();
  }

  return _in->good();
}

bool SAMParser::map_end_from_line(char* line) {
  ReadHit& r = *_read_buff;
  string sam_line(line);
  char *p = strtok(line, "\t");
  int sam_flag = 0;
  bool paired = 0;
  bool left_first = 0;
  bool other_reversed = 0;

  int i = 0;
  while (p && i <= 9) {
    switch(i++) {
      case 0: {
        r.name = p;
        if (boost::algorithm::ends_with(r.name, "\1") ||
            boost::algorithm::ends_with(r.name, "\2")) {
          r.name = r.name.substr(r.name.size()-2);
        }
        break;
      }
      case 1: {
        sam_flag = atoi(p);
        if (sam_flag & 0x4) {
          goto stop;
        }
        paired = sam_flag & 0x1;
        if (paired && (direction == F || direction == R)) {
          goto stop;
        }
        if (paired && (!(sam_flag & 0x2) || sam_flag & 0x8)) {
          goto stop;
        }
        r.reversed = sam_flag & 0x10;
        other_reversed = sam_flag & 0x20;
        if (paired && r.reversed == other_reversed) {
          goto stop;
        }
        r.first = !paired || sam_flag & 0x40;
        left_first = ((!paired && !r.reversed) || (r.first && !r.reversed) ||
                        (!r.first && r.reversed));
        if (((direction == RF || direction == R) && left_first) ||
            ((direction == FR || direction == F) && !left_first)) {
  		    goto stop;
        }
        break;
      }
      case 2: {
        if(p[0] == '*') {
          goto stop;
        }
        if (!_targ_index.count(p)) {
          logger.severe("Target sequence '%s' not found. Verify that it is in "
                        "the SAM/BAM header and FASTA file.", p);
        }
        r.targ_id = _targ_index[p];
        break;
      }
      case 3: {
        r.left = (size_t)(atoi(p)-1);
        break;
      }
      case 4: {
        break;
      }
      case 5: {
        r.right = r.left + cigar_length(p, r.inserts, r.deletes);
        foreach (Indel& indel, r.inserts) {
          if (indel.len > max_indel_size) {
            goto stop;
          }
        }
        foreach (Indel& indel, r.deletes) {
          if (indel.len > max_indel_size) {
            goto stop;
          }
        }
        break;
      }
      case 6: {
        // skip if only one end mapped of paired-end read
        if(paired && p[0] != '=') {
          goto stop;
        }
        break;
      }
      case 7: {
        r.mate_l = atoi(p)-1;
        if (paired && ((r.reversed && r.left < (size_t)r.mate_l) ||
                       (other_reversed && r.left > (size_t)r.mate_l))) {
          goto stop;
        }
        break;
      }
      case 8: {
        break;
      }
      case 9: {
        r.seq.set(p, r.reversed);
        r.sam = sam_line;
        goto stop;
      }
    }
    p = strtok(NULL, "\t");
  }
 stop:
  return i == 10;
}

void SAMParser::reset() {
  // Rewind input file
  _in->clear();
  _in->seekg(0, ios::beg);

  // Load first alignment
  char line_buff[BUFF_SIZE];
  delete _read_buff;
  _read_buff = new ReadHit();

  while(_in->good()) {
    _in->getline(line_buff, BUFF_SIZE-1, '\n');
    if (line_buff[0] != '@') {
      break;
    }
  }

  while(!map_end_from_line(line_buff)) {
    _in->getline(line_buff, BUFF_SIZE-1, '\n');
  }
}

BAMWriter::BAMWriter(BamTools::BamWriter* writer, bool sample)
   : _writer(writer) {
  _sample = sample;
}

BAMWriter::~BAMWriter() {
  _writer->Close();
}

void BAMWriter::write_fragment(Fragment& f) {
  if (_sample) {
    const FragHit* hit = f.sample_hit();
    PairStatus ps = hit->pair_status();
    if (ps != RIGHT_ONLY) {
      _writer->SaveAlignment(hit->left_read()->bam);
    }
    if (ps != LEFT_ONLY) {
      _writer->SaveAlignment(hit->right_read()->bam);
    }
  } else {
    double total = 0;
    foreach(FragHit* hit, f.hits()) {
      total += sexp(hit->params()->posterior);
      PairStatus ps = hit->pair_status();
      if (ps != RIGHT_ONLY) {
        hit->left_read()->bam.AddTag("XP","f",(float)sexp(hit->params()->posterior));
        _writer->SaveAlignment(hit->left_read()->bam);
      }
      if (ps != LEFT_ONLY) {
        hit->right_read()->bam.AddTag("XP","f",(float)sexp(hit->params()->posterior));
        _writer->SaveAlignment(hit->right_read()->bam);
      }
    }
    assert(approx_eq(total, 1.0));
  }
}

SAMWriter::SAMWriter(ostream* out, bool sample) : _out(out) {
  _sample = sample;
}

SAMWriter::~SAMWriter() {
  _out->flush();
}

void SAMWriter::write_fragment(Fragment& f) {
  if (_sample) {
    const FragHit* hit = f.sample_hit();
    PairStatus ps = hit->pair_status();

    if (ps != RIGHT_ONLY) {
      *_out << hit->left_read()->sam << endl;
    }
    if (ps != LEFT_ONLY) {
      *_out << hit->right_read()->sam << endl;
    }
  } else {
    double total = 0;
    foreach(const FragHit* hit, f.hits()) {
      total += sexp(hit->params()->posterior);
      PairStatus ps = hit->pair_status();
      if (ps != RIGHT_ONLY) {
        *_out << hit->left_read()->sam << " XP:f:"
              << (float)sexp(hit->params()->posterior) << endl;
      }
      if (ps != LEFT_ONLY) {
        *_out << hit->right_read()->sam << " XP:f:"
              << (float)sexp(hit->params()->posterior) << endl;
      }
    }
    assert(approx_eq(total, 1.0));
  }
}
