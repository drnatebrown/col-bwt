
struct LCS_Info {
    ulint lcs;
    ulint k;
    ulint d;

    LCS_Info(ulint k, ulint d, ulint curr_lcs) : k(k), d(d), curr_lcs(curr_lcs) {}
};

class LCS_Iter {
    ulint k, d, curr_lcs;
    std::vector<FL_run> FL_runs; // Assuming FL_run is a struct or class you have defined
    bool hasNext;

public:
    LCS_Iter(std::vector<FL_run> runs) : FL_runs(runs), k(0), d(0), curr_lcs(0), hasNext(true) {}

    void begin() {
        k = 0;
        d = 0;
        curr_lcs = 0;
        hasNext = true;
    }

    std::tuple<ulint, ulint, ulint> next() {
        if (!hasNext) {
            throw std::out_of_range("End of iterator reached.");
        }

        if (d == 0) {
            curr_lcs = 0;
        }
        else {
            ++curr_lcs;
        }

        std::pair<ulint, ulint> new_pos = FL(k, d); // Assuming FL is a function you have defined
        k = new_pos.first;
        d = new_pos.second;

        if (k == 0 && d == 0) {
            hasNext = false;
        }

        return LCS_Info(k, d, curr_lcs);
    }

    bool has_next() {
        return hasNext;
    }
};