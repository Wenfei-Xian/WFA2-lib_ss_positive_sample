#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <zlib.h>
#include <cctype>
#include "bindings/cpp/WFAligner.hpp"

using namespace std;
using namespace wfa;

// 用于保存 FASTQ 记录的结构体
struct FastqRecord {
    string header;
    string seq;
    string plus;
    string qual;
};

// 读取一条 FASTQ 记录（4行）
bool read_record(gzFile file, FastqRecord &record) {
    char buffer[4096];
    if (gzgets(file, buffer, sizeof(buffer)) == nullptr) return false;
    record.header = buffer;
    if (!record.header.empty() && record.header.back() == '\n')
        record.header.pop_back();
    
    if (gzgets(file, buffer, sizeof(buffer)) == nullptr) return false;
    record.seq = buffer;
    if (!record.seq.empty() && record.seq.back() == '\n')
        record.seq.pop_back();
    
    if (gzgets(file, buffer, sizeof(buffer)) == nullptr) return false;
    record.plus = buffer;
    if (!record.plus.empty() && record.plus.back() == '\n')
        record.plus.pop_back();
    
    if (gzgets(file, buffer, sizeof(buffer)) == nullptr) return false;
    record.qual = buffer;
    if (!record.qual.empty() && record.qual.back() == '\n')
        record.qual.pop_back();
    
    return true;
}

// 从 FASTQ 头部（去掉开头的 '@'）提取 ID
string get_id(const string &header) {
    return (header.empty() || header[0] != '@') ? header : header.substr(1);
}

// 检查两个 ID 是否构成合法的正反链对：
// 去掉最后 3 个字符后应相同，且第一个记录后缀为 "fwd"，第二个为 "rev"
bool valid_pair(const string &id1, const string &id2) {
    if (id1.size() < 3 || id2.size() < 3) return false;
    string base1 = id1.substr(0, id1.size()-3);
    string base2 = id2.substr(0, id2.size()-3);
    string suf1 = id1.substr(id1.size()-3);
    string suf2 = id2.substr(id2.size()-3);
    return (base1 == base2 && suf1 == "fwd" && suf2 == "rev");
}

// 返回 DNA 碱基的互补
char complement(char base) {
    switch(base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default: return base;
    }
}

// 计算 DNA 序列的反向互补
string reverse_complement(const string &seq) {
    string rc(seq.rbegin(), seq.rend());
    for (char &base : rc)
        base = complement(base);
    return rc;
}

// 直接基于 CIGAR 字符串处理：
// 对于每个错配 X，取其左侧最多 window 个字符和右侧最多 window 个字符，
// 统计其中 M 的个数；如果左侧 M 数等于 window 且右侧 M 数也等于 window，输出 "#Good <pattern_id> <mismatch_index>"
void report_x_positions_window_simple(const string &cigar, const string &pattern_id, int window = 50) {
    for (size_t i = 0; i < cigar.size(); i++) {
        if (cigar[i] == 'X') {
            size_t left_start = (i < static_cast<size_t>(window) ? 0 : i - static_cast<size_t>(window));
            string left_window = cigar.substr(left_start, i - left_start);
            int leftM = 0;
            for (char c : left_window) {
                if (c == 'M')
                    leftM++;
            }
            size_t right_len = ((i + 1 + static_cast<size_t>(window)) > cigar.size() ? cigar.size() - (i + 1) : static_cast<size_t>(window));
            string right_window = cigar.substr(i + 1, right_len);
            int rightM = 0;
            for (char c : right_window) {
                if (c == 'M')
                    rightM++;
            }
            cout << "Mismatch (X) at cigar index " << i 
                 << ": left window M count = " << leftM 
                 << ", right window M count = " << rightM << endl;
            if (leftM == window && rightM == window) {
                cout << "#Good " << pattern_id << " " << i << endl;
            }
        }
    }
}

// 根据 CIGAR 字符串重构对齐结果（简单版本），以类似 C 示例格式显示
void pretty_print_alignment_display(const string &pattern, const string &query, const string &cigar) {
    string aligned_pattern, aligned_query;
    size_t i = 0, j = 0;
    for (char op : cigar) {
        if (op == 'M' || op == 'X') {
            if (i < pattern.size() && j < query.size()) {
                aligned_pattern.push_back(pattern[i]);
                aligned_query.push_back(query[j]);
                i++; j++;
            }
        } else if (op == 'I') {
            if (j < query.size()) {
                aligned_pattern.push_back('-');
                aligned_query.push_back(query[j]);
                j++;
            }
        } else if (op == 'D') {
            if (i < pattern.size()) {
                aligned_pattern.push_back(pattern[i]);
                aligned_query.push_back('-');
                i++;
            }
        }
    }
    cout << "  PATTERN  " << aligned_pattern << "\n";
    cout << "  TEXT     " << aligned_query << "\n";
}

// 使用 WFA2 C++ 绑定进行比对：
// 使用 local alignment (extension alignment)，将 fwd 作为 reference，
// 将 rev（经过反向互补）作为 query；pattern_id 为 fwd 的 ID
void align_wfa(const string &fwd, const string &rev, const string &pattern_id) {
    // 创建可修改的局部变量，用于 alignExtension
    string pattern = fwd;  // reference
    string query = rev;    // query（已为反向互补）
    
    // 创建 Gap-Affine 模型比对对象：mismatch=4, gap-opening=6, gap-extension=2
    WFAlignerGapAffine aligner(4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh);
    
    // 执行 extension (local) 比对
    aligner.alignExtension(pattern, query);
    
    cout << "  PATTERN  " << pattern << "\n";
    cout << "  TEXT     " << query << "\n";
    cout << "  SCORE (RE)COMPUTED " << aligner.getAlignmentScore() << "\n";
    
    string cigar = aligner.getAlignment();
    cout << "  CIGAR    " << cigar << "\n";
    
    // 显示对齐结果（pretty print）
    pretty_print_alignment_display(pattern, query, cigar);
    
    // 直接基于 CIGAR 字符串处理，输出每个 X 左右 50 个字符内 M 的数量
    report_x_positions_window_simple(cigar, pattern_id, 50);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " input.fastq.gz" << endl;
        return 1;
    }
    
    const char* filename = argv[1];
    gzFile file = gzopen(filename, "r");
    if (!file) {
        cout << "Error opening file: " << filename << endl;
        return 1;
    }
    
    FastqRecord rec1, rec2;
    if (!read_record(file, rec1)) {
        cout << "No records found in file." << endl;
        gzclose(file);
        return 1;
    }
    
    // 主循环：每次用 rec1 与 rec2 配对
    while (read_record(file, rec2)) {
        string id1 = get_id(rec1.header);
        string id2 = get_id(rec2.header);
        
        if (valid_pair(id1, id2)) {
            cout << "Processing ZMW pair: " << id1 << " & " << id2 << "\n";
            // 对于合法对，fwd 记录直接使用，rev 记录取反向互补后再比对
            string rev_rc = reverse_complement(rec2.seq);
            align_wfa(rec1.seq, rev_rc, id1);
            if (!read_record(file, rec1)) break;
        } else {
            cout << "Skipping record with id: " << id1 
                 << " (pair invalid with id: " << id2 << ")\n";
            rec1 = rec2;
        }
    }
    
    gzclose(file);
    return 0;
}

