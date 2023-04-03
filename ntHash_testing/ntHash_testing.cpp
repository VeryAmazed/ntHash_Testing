#include "nthash.hpp"
#include <fstream>

using namespace std;
#define pb push_back 
#define pii pair<int, int>
#define mp make_pair
#define f first
#define s second
int main(){
    string str1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    string str2 = "AGAGATTACGTCTGGTTGCAAGAGATCATAACAGGGGAAATTGATTGAAAATAAATATATCGCCAGCAGCACATGAACAAGTTTCGGAATGTGATCAATT"; // first 100 characters of some salmonella strand 
    string str3 = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTC"; // first 100 charcters of some covid-19 strand
    vector<string> vec;
    vec.pb(str1);
    vec.pb(str2);
    vec.pb(str3);

    fstream out;
    out.open("../output.txt", fstream::out);
    vector<vector<pii>> forwards(3, vector<pii>()); // vector to store all the k, pos combinations that cause the forward hashes to differ
    vector<vector<pii>> reverses(3, vector<pii>()); // vector to store all the k, pos combinations that cause the reverse hashes to differ
    uint64_t obj_fhash; // forward hash produced by the NtHash Object
    uint64_t obj_rhash; // reverse hash produced by the NtHash Object
    uint64_t ntf64_fhash; // forward hash produced by ntf64(seq, k)
    uint64_t ntr64_rhash; // reverse hash produced by ntr64(seq, k)

    for(unsigned k =1; k <= 64; k++){ // k-mer sizes that will be tested, 1 to 64
        for(size_t i =0; i < 10; i++){ // starting positions within the string to be tested
            for(int j =0; j < 3; j++){ // iterating across the 3 strings
                nthash::NtHash nt(vec[j], 1, k, i);
                nt.roll(); 
                obj_fhash = nt.get_forward_hash();
                obj_rhash = nt.get_reverse_hash();

                ntf64_fhash = nthash::ntf64((vec[j].c_str()) + i, k);
                ntr64_rhash = nthash::ntr64((vec[j].c_str()) + i, k);

                if(obj_fhash != ntf64_fhash){
                    forwards[j].pb(mp(k, i));
                    out << "Object/ntmc64 fhash != ntf64 fhash" << endl;
                    out << obj_fhash << endl;
                    out << ntf64_fhash << endl;
                    out << "k = " << k << " pos = " << i << " on string " << j+1 << endl;
                    out << endl;
                }

                if(obj_rhash != ntr64_rhash){
                    reverses[j].pb(mp(k, i));
                    out << "Object/ntmc64 rhash != ntr64 rhash" << endl;
                    out << obj_rhash << endl;
                    out << ntr64_rhash << endl;
                    out << "k = " << k << " pos = " << i << " on string " << j+1 << endl;
                    out << endl;
                }
            }
        }
    }
     out << endl;

    for(int i =0; i < 3; i++){
        out << "Forward Hash Differences String " << i+1 << endl;
        out << "Size: " << forwards[i].size() << endl;
        for(size_t j =0; j < forwards[i].size(); j++){
            out << "(" << forwards[i][j].f << ", " << forwards[i][j].s << ") ";
        }
        out << endl;
    }

    for(int i =0; i < 3; i++){
        out << "Reverse Hash Differences String " << i+1 << endl;
        out << "Size: " << reverses[i].size() << endl;
        for(size_t j =0; j < reverses[i].size(); j++){
            out << "(" << reverses[i][j].f << ", " << reverses[i][j].s << ") ";
        }
        out << endl;
    }
    out.close();
    //cout << "Hello, got here" << endl;
    return 0;
}