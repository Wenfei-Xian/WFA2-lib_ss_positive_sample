#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <zlib.h>
#include <cctype>
#include <vector>
#include <algorithm>
#include "bindings/cpp/WFAligner.hpp"

using namespace std;
using namespace wfa;

// ---------------------------
// FASTQ 记录结构体及辅助函数
// ---------------------------
struct FastqRecord {
    string header;
    string seq;
    string plus;
    string qual;
};

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

string get_id(const string &header) {
    return (header.empty() || header[0] != '@') ? header : header.substr(1);
}

bool valid_pair(const string &id1, const string &id2) {
    if (id1.size() < 3 || id2.size() < 3) return false;
    string base1 = id1.substr(0, id1.size()-3);
    string base2 = id2.substr(0, id2.size()-3);
    string suf1 = id1.substr(id1.size()-3);
    string suf2 = id2.substr(id2.size()-3);
    return (base1 == base2 && suf1 == "fwd" && suf2 == "rev");
}

// ---------------------------
// 序列及质量处理函数
// ---------------------------
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
        default:  return base;
    }
}

string reverse_complement(const string &seq) {
    string rc(seq.rbegin(), seq.rend());
    for (char &base : rc)
        base = complement(base);
    return rc;
}

string reverse_string(const string &s) {
    return string(s.rbegin(), s.rend());
}

// ---------------------------
// 对齐结果结构体及重构函数
// ---------------------------
struct AlignmentResult {
    string alignedPattern;       // 重构后的参考序列
    string alignedQuery;         // 重构后的查询序列
    string alignedPatternQual;   // 参考序列对应的质量字符串
    string alignedQueryQual;     // 查询序列对应的质量字符串
    string alignedOps;           // 每个对齐位置的操作（M, X, I, D）
};

// 根据原始序列、质量和 CIGAR 重构对齐结果（假设 CIGAR 为已展开形式，每个操作一个字符）
AlignmentResult reconstruct_alignment_with_quality(const string &pattern,
                                                     const string &query,
                                                     const string &patternQual,
                                                     const string &queryQual,
                                                     const string &cigar) {
    AlignmentResult res;
    size_t i = 0, j = 0;
    for (char op : cigar) {
        if (op == 'M' || op == 'X') {
            if (i < pattern.size() && j < query.size()) {
                res.alignedPattern.push_back(pattern[i]);
                res.alignedQuery.push_back(query[j]);
                res.alignedPatternQual.push_back(patternQual[i]);
                res.alignedQueryQual.push_back(queryQual[j]);
                res.alignedOps.push_back(op);
                i++; j++;
            }
        } else if (op == 'I') {
            if (j < query.size()) {
                res.alignedPattern.push_back('-');
                res.alignedQuery.push_back(query[j]);
                res.alignedPatternQual.push_back(' '); // 空白代表无质量
                res.alignedQueryQual.push_back(queryQual[j]);
                res.alignedOps.push_back(op);
                j++;
            }
        } else if (op == 'D') {
            if (i < pattern.size()) {
                res.alignedPattern.push_back(pattern[i]);
                res.alignedQuery.push_back('-');
                res.alignedPatternQual.push_back(patternQual[i]);
                res.alignedQueryQual.push_back(' ');
                res.alignedOps.push_back(op);
                i++;
            }
        }
    }
    return res;
}

// 辅助函数：将质量字符转换为数字（ASCII - 33）
int qual_to_num(char q) {
    return static_cast<int>(q) - 33;
}

// ---------------------------
// 对齐结果输出及错配质量检测
// ---------------------------
// 输出对齐结果，包括序列和对应的质量数值（数字间用空格分隔）
void print_alignment_with_quality(const AlignmentResult &ar) {
    cout << "  PATTERN  " << ar.alignedPattern << "\n";
    cout << "  TEXT     " << ar.alignedQuery << "\n";
    
    // 输出质量数值：对于参考序列
    cout << "  QUAL(P)  ";
    for (char q : ar.alignedPatternQual) {
        if(q != ' ')
            cout << qual_to_num(q) << " ";
        else
            cout << "  ";
    }
    cout << "\n";
    
    // 输出质量数值：对于查询序列
    cout << "  QUAL(Q)  ";
    for (char q : ar.alignedQueryQual) {
        if(q != ' ')
            cout << qual_to_num(q) << " ";
        else
            cout << "  ";
    }
    cout << "\n";
    
    cout << "  OPS      " << ar.alignedOps << "\n";
}

// 遍历对齐结果的操作串，针对每个 'X'（错配），
// 统计其在对齐结果中左侧最多50个字符和右侧最多50个字符内出现的 'M' 数量；
// 同时检查该位置在参考与查询的对齐质量值是否均等于93；
// 输出信息时使用 PATTERN 的对齐位置（1-based）。
void report_mismatch_with_quality(const AlignmentResult &ar, const string &pattern_id, int window = 50) {
    int n = ar.alignedOps.size();
    for (int pos = 0; pos < n; pos++) {
        if (ar.alignedOps[pos] == 'X') {
            // PATTERN index（1-based），忽略 gap（'-'）在参考中均已由重构时处理
            int pattern_index = pos + 1;
            
            // 统计左侧窗口中 'M' 数量（窗口范围 [max(0, pos-window), pos)）
            int leftM = 0;
            int left_start = max(0, pos - window);
            for (int k = left_start; k < pos; k++) {
                if (ar.alignedOps[k] == 'M')
                    leftM++;
            }
            // 统计右侧窗口中 'M' 数量（窗口范围 (pos, min(n, pos+window+1))）
            int rightM = 0;
            int right_end = min(n, pos + window + 1);
            for (int k = pos + 1; k < right_end; k++) {
                if (ar.alignedOps[k] == 'M')
                    rightM++;
            }
            
            // 获取该错配位置的质量值（参考与查询），转换为数字
            int qualP = (ar.alignedPatternQual[pos] != ' ' ? qual_to_num(ar.alignedPatternQual[pos]) : 0);
            int qualQ = (ar.alignedQueryQual[pos] != ' ' ? qual_to_num(ar.alignedQueryQual[pos]) : 0);
            
            cout << "Mismatch (X) at PATTERN index " << pattern_index
                 << ": left window M count = " << leftM
                 << ", right window M count = " << rightM
                 << ", qualities: (P) " << qualP << " (Q) " << qualQ << endl;
            
            // 如果左右窗口 M 数各等于 window 且两侧质量均为 93，则输出 #Good
            if (leftM == window && rightM == window && qualP == 93 && qualQ == 93) {
                cout << "#Good " << pattern_id << " " << pattern_index << endl;
            }
        }
    }
}

// ---------------------------
// 对齐并输出结果
// ---------------------------
// 使用 WFA2 C++ 绑定进行 local alignment (extension alignment)，
// 将 fwd 作为 reference，rev（经过反向互补及质量反转）作为 query；
// pattern_id 为 fwd 的 ID；同时传入 fwd 质量字符串和 rev 质量字符串（rev 质量须反转）
void align_wfa(const string &fwd, const string &rev, const string &fwdQual, const string &revQual, const string &pattern_id) {
    // 为 alignExtension 创建可修改的局部变量
    string pattern = fwd;  // reference
    string query = rev;    // query（已为反向互补）
    
    // 创建 Gap-Affine 模型比对对象：mismatch=4, gap-opening=6, gap-extension=6（与之前保持一致）
    WFAlignerGapAffine aligner(4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh);
    aligner.alignExtension(pattern, query);
    
    cout << "  PATTERN  " << pattern << "\n";
    cout << "  TEXT     " << query << "\n";
    cout << "  SCORE (RE)COMPUTED " << aligner.getAlignmentScore() << "\n";
    
    string cigar = aligner.getAlignment();
    cout << "  CIGAR    " << cigar << "\n";
    
    // 使用传入的原始质量字符串重构对齐结果
    AlignmentResult ar = reconstruct_alignment_with_quality(pattern, query, fwdQual, revQual, cigar);
    
    print_alignment_with_quality(ar);
    report_mismatch_with_quality(ar, pattern_id, 50);
}

// ---------------------------
// 主函数
// ---------------------------
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
            // 对于合法对，fwd 记录直接使用；对于 rev 记录，取反向互补及反转质量字符串
            string rev_rc = reverse_complement(rec2.seq);
            string revQual_rev = reverse_string(rec2.qual);
            align_wfa(rec1.seq, rev_rc, rec1.qual, revQual_rev, id1);
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

