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
    cout << "No alignment file specified. Expecting streaming input on stdin...\n\n";
    _parser.reset(new SAMParser(&cin));
    is_sam = true;
  } else {
    cout << "Attempting to read '" << in_file << "' in BAM format...\n";
    BamTools::BamReader* reader = new BamTools::BamReader();
    if (reader->Open(in_file)) {
      cout << "Parsing BAM header...\n";
      _parser.reset(new BAMParser(reader));
      if (out_file.size()) {
        out_file += ".bam";
        BamTools::BamWriter* writer = new BamTools::BamWriter();
        if (writer->Open(out_file, reader->GetHeader(),
                         reader->GetReferenceData())) {
          bool sample = out_file.substr(out_file.length()-8,4) == "samp";
          _writer.reset(new BAMWriter(writer, sample));
        } else {
          cerr << "ERROR: Unable to open output BAM file '"
               << out_file << "'.\n";
          exit(1);
        }
      }
    } else {
      delete reader;
      cout << "Input is not in BAM format. Trying SAM...\n";
      ifstream* ifs = new ifstream(in_file.c_str());
      if (!ifs->is_open()) {
        cerr << "ERROR: Unable to open input SAM file '" << in_file << "'.\n";
        exit(1);
      }
      _parser.reset(new SAMParser(ifs));
      is_sam = true;
    }
  }

  if (is_sam && out_file.size()) {
    out_file += ".sam";
    ofstream* ofs = new ofstream(out_file.c_str());
    if (!ofs->is_open()) {
      cerr << "ERROR: Unable to open output SAM file '" << out_file << "'.\n" ;
      exit(1);
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
      Target* t = targ_table.get_targ(m.targ_id);
      if (!t) {
        cerr << "ERROR: Target sequence at index '" << m.targ_id
             << "' not found. Verify that it is in the SAM/BAM header and "
             << "FASTA file.\n";
        exit(1);
      }
      m.targ = t;
      assert(t->id() == m.targ_id);
      
      // Add num_neighbors targets on either side to the neighbors list.
      // Used for experimental feature.
      for (TargID j = 1; j <= num_neighbors;  j++) {
        if (j <= m.targ_id) {
          m.neighbors.push_back(targ_table.get_targ(m.targ_id - j));
        }
        if (j + m.targ_id < targ_table.size()) {
          m.neighbors.push_back(targ_table.get_targ(m.targ_id + j));
        }
      }
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
    _targ_index[ref.RefName] = index++;
    _targ_lengths[ref.RefName] = ref.RefLength;
  }
    
  // Get first valid FragHit
  _frag_buff = new FragHit();
  do {
    if (!_reader->GetNextAlignment(a)) {
      cerr << "ERROR: Input BAM file contains no valid alignments.\n";
      exit(1);
    }
  } while(!map_end_from_alignment(a));
}

bool BAMParser::next_fragment(Fragment& nf) {
  nf.add_map_end(_frag_buff);
    
  BamTools::BamAlignment a;
  _frag_buff = new FragHit();

  while(true) {
    if (!_reader->GetNextAlignment(a)) {
      return false;
    } else if (!map_end_from_alignment(a)) {
      continue;
    } else if (!nf.add_map_end(_frag_buff)) {
      return true;
    }
    _frag_buff = new FragHit();
  }
}

bool BAMParser::map_end_from_alignment(BamTools::BamAlignment& a) {
  FragHit& f = *_frag_buff;
    
  if (!a.IsMapped()) {
    return false;
  }
    
  if (a.IsPaired() && (!a.IsMateMapped() || a.RefID != a.MateRefID ||
                       a.IsReverseStrand() == a.IsMateReverseStrand())) {
    return false;
  }
  
  f.left_first = (!a.IsPaired() && !a.IsReverseStrand()) ||
                 (a.IsFirstMate() && !a.IsReverseStrand())||
                 (a.IsSecondMate() && a.IsReverseStrand());
    
  if ((direction == RF && f.left_first) || (direction == FR && !f.left_first)) {
    return false;
  }
 
  f.name = a.Name;
  f.targ_id = a.RefID;
  f.left = a.Position;
  f.mate_l = a.MatePosition;

  if (a.IsReverseStrand()) {
    f.seq_r.set(a.QueryBases, 1);
    f.bam_r = a;
    f.right = f.left + cigar_length(a.CigarData, f.inserts_r, f.deletes_r);
  } else {
    f.seq_l.set(a.QueryBases, 0);
    f.bam_l = a;
    f.right = f.left + cigar_length(a.CigarData, f.inserts_l, f.deletes_l);
  }
    
  return true;
}

void BAMParser::reset() {
  _reader->Rewind();
    
  // Get first valid FragHit
  BamTools::BamAlignment a;
  _frag_buff = new FragHit();
  do {
    _reader->GetNextAlignment(a);
  } while(!map_end_from_alignment(a));
}


SAMParser::SAMParser(istream* in) {
  _in = in;
  
  char line_buff[BUFF_SIZE];
  _frag_buff = new FragHit();
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
      if (_targ_index.find(name) == _targ_index.end()) {
        _targ_index[name] = index++;
      } else {
        cerr << "Warning: Target '" << str
             << "' appears twice in the SAM index.\n";
      }
      
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
      cerr << "ERROR: Input SAM file contains no valid alignments.\n";
      exit(1);
    }
    _in->getline(line_buff, BUFF_SIZE-1, '\n');
  }
}

bool SAMParser::next_fragment(Fragment& nf) {
  nf.add_map_end(_frag_buff);
  
  _frag_buff = new FragHit();
  char line_buff[BUFF_SIZE];
  
  while(_in->good()) {
    _in->getline (line_buff, BUFF_SIZE-1, '\n');
    if (!map_end_from_line(line_buff)) {
      continue;
    }
    if (!nf.add_map_end(_frag_buff)) {
      break;
    }
    _frag_buff = new FragHit();
  }
  
  return _in->good();
}

bool SAMParser::map_end_from_line(char* line) {
  FragHit& f = *_frag_buff;
  string sam_line(line);
  char *p = strtok(line, "\t");
  int sam_flag=0;
  bool paired=0;
  
  int i = 0;
  while (p && i <= 9) {
    switch(i++) {
      case 0: {
        f.name = p;
        break;
      }
      case 1: {
        sam_flag = atoi(p);
        if (sam_flag & 0x4) {
          goto stop;
        }
        paired = sam_flag & 0x1;
        if (paired && (sam_flag & 0x8)) {
          goto stop;
        }
        bool reversed = sam_flag & 0x10;
        bool other_reversed = sam_flag & 0x20;
        if (paired && reversed == other_reversed) {
          goto stop;
        }
        bool first = sam_flag & 0x40;
        f.left_first = ((!paired && !reversed) || (first && !reversed) ||
                        (!first && reversed));
        if ((direction == RF && f.left_first) ||
            (direction == FR && !f.left_first)) {
  		    goto stop;
        }
        break;
      }
      case 2: {
        if(p[0] == '*') {
          goto stop;
        }
        if (!_targ_index.count(p)) {
          cerr << "ERROR: Target sequence '" << p << "' not found. Verify that "
              << "it is in the SAM/BAM header and FASTA file.\n";
          exit(1);
        }
        f.targ_id = _targ_index[p];
        break;
      }
      case 3: {
        f.left = (size_t)(atoi(p)-1);
        break;
      }
      case 4: {
        break;
      }
      case 5: {
        if (sam_flag & 0x10) {
          f.right = f.left + cigar_length(p, f.inserts_r, f.deletes_r);
        } else {
          f.right = f.left + cigar_length(p, f.inserts_l, f.deletes_l);
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
        f.mate_l = atoi(p)-1;
        break;
      }
      case 8: {
        break;
      }
      case 9: {
        if (sam_flag & 0x10) {
          f.seq_r.set(p, 1);
          f.sam_r = sam_line;
        } else {
          f.seq_l.set(p, 0);
          f.sam_l = sam_line;
        }
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
  _frag_buff = new FragHit();
  
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
      _writer->SaveAlignment(hit->bam_l);
    }
    if (ps != LEFT_ONLY) {
      _writer->SaveAlignment(hit->bam_r);
    }
  } else {
    double total = 0;
    foreach(FragHit* hit, f.hits()) {
      total += hit->probability;
      PairStatus ps = hit->pair_status();
      if (ps != RIGHT_ONLY) {
        hit->bam_l.AddTag("XP","f",(float)hit->probability);
        _writer->SaveAlignment(hit->bam_l);
      }
      if (ps != LEFT_ONLY) {
        hit->bam_r.AddTag("XP","f",(float)hit->probability);
        _writer->SaveAlignment(hit->bam_r);
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
      *_out << hit->sam_l << endl;
    }
    if (ps != LEFT_ONLY) {
      *_out << hit->sam_r << endl;
    }
  } else {
    double total = 0;
    foreach(FragHit* hit, f.hits()) {
      total += hit->probability;
      PairStatus ps = hit->pair_status();
      if (ps != RIGHT_ONLY) {
        *_out << hit->sam_l << " XP:f:" << (float)hit->probability << endl;
      }
      if (ps != LEFT_ONLY) {
        *_out << hit->sam_r << " XP:f:" << (float)hit->probability << endl;
      }
    }
    assert(approx_eq(total, 1.0));
  }
}
