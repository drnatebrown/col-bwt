#include <zlib.h>
#include <kseq.h>
KSEQ_INIT(gzFile, gzread);

//*********************** Fasta reader ***************************************
class PatternProcessor {
public:
    PatternProcessor(const std::string &filename) : l(0) {
        gzFile fp = gzopen(filename.data(), "r");
        seq = kseq_init(fp);
    }

    bool read() {
        return (l = kseq_read(seq)) >= 0;
    }

    std::pair<const char*, size_t> get_seq() {
        if (l >= 0) {
            return std::make_pair(seq->seq.s, seq->seq.l);
        }
        return std::make_pair(nullptr, 0);
    }

    std::string get_id() {
        return std::string(seq->name.s, seq->name.l);
    }

    std::string get_string() {
        return std::string(get_seq().first, get_seq().second);
    }

private:
    kseq_t* seq;
    int l;
};
//********** end fasta reader ********************