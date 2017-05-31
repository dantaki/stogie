#include <string>
#include <vector>
#include "api/BamReader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
using namespace std;
using namespace BamTools;
inline std::vector<std::string> split(std::string &s, char delim)
{
  std::vector<std::string> v;
  std::stringstream ss(s); std::string i;
  while(std::getline(ss, i, delim)) { v.push_back(i); }
  return v;
}
map <string, int> hashRef (const RefVector &ref)
{
        map<string, int> hash;
        for(uint16_t i=0; i!= ref.size(); ++i){ hash[ref[i].RefName]=i;}
        return hash;
} 
inline double  overlap(int32_t &s1, int32_t &e1, int32_t &s2, int32_t &e2)
{
        int s[2]={s1,s2}; int e[2] = {e1,e2};
        std::sort(s,s+2); std::sort(e,e+2);
        int ovr = e[0]-s[1]+1;
        float o[2] = { (float)ovr/(e2-s2+1), (float)ovr/(e1-s1+1)};
        std::sort(o,o+2);
        return o[0];
}
int main(int argc, char *argv[])
{
	string splash=	
				"\n8\"\"\"\"8    \"\"8\"\"    8\"\"\"88    8\"\"\"\"8    8     8\"\"\"\" \n"
				"8           8      8    8    8    \"    8     8     \n"
				"8eeeee      8e     8    8    8e        8e    8eeee \n"
				"    88      88     8    8    88  ee    88    88    \n"
				"e   88      88     8    8    88   8    88    88    \n"
				"8eee88      88     8eeee8    88eee8    88    88eee \n\n"	
				"stogie         Validate SV with CIGAR strings\n"
				"Version: 1.0	Author: Danny Antaki <dantaki@ucsd.edu>\n\n"
				"Usage: stogie -i <in.bam> -r <sv.bed> -x <FLOAT> -q <INT> -o <output.txt>\n\nOptions:\n"
				"    -i        Input: BAM filename\n"
				"    -r        Input: SV bed file\n"
				"    -x        Minimum reciprocal overlap [0.5]\n" 
				"    -q        Mapping quality threshold [10]\n"
				"    -o        Output: filename\n";
	string ifh; string bed; string ofh; int Q=10; double OVR=0.5;
	if ( (argc==1) ||
	     (argc==2 && string(argv[1]) == "-h") ||
	     (argc==2 && string(argv[1]) == "-help") ||
	     (argc==2 && string(argv[1]) == "--help"))
		{ cerr << splash << endl; return 1;}
	for(int i=1; i<argc; i++){
		if(string(argv[i]) == "-i" || string(argv[i]) == "--in" || string(argv[i]) == "-in"){ ifh = string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-l" || string(argv[i]) == "-r" || string(argv[i]) == "--l" || string(argv[i])=="--r" || string(argv[i])=="-bed" || string(argv[i])=="--bed"){ bed=string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-o" || string(argv[i]) == "--out" || string(argv[i]) == "-out"){ ofh = string(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-q" || string(argv[i]) == "--q"){ Q = atoi(argv[i+1]); i++; continue; }
		if(string(argv[i]) == "-x" || string(argv[i]) == "--x"){ OVR = ::atof(argv[i+1]); i++; continue; }
		cerr << "ERROR: Unknown option "<<string(argv[i])<< '\n' << splash <<endl;
		return 1;
	}
	if(ifh == "") { cerr << "ERROR: No BAM file given"<< endl; return 1; }
	if(bed == "") { cerr << "ERROR: No BED file given"<< endl; return 1; }
	if(ofh == ""){ ofh = ifh.substr(0,ifh.find_last_of('.'))+"_stogie.txt";  }		
	ifstream fin(bed.c_str());
	BamReader bam;
	if(!fin) { cerr << "ERROR: Cannot open " << ifh << endl; return 1; } 
	if (!bam.Open(ifh)){cerr << "ERROR: " << ifh << " could not be opened!"<< endl; return 1; }
	BamAlignment al;
	const RefVector refs = bam.GetReferenceData();
	map<string, int> chrom = hashRef(refs);
	bool chrFlag=0; if (!refs[0].RefName.compare(0, 3, "chr")){ chrFlag=1; } //chrFlag=true when reference is prefixed with "chr"  
	ofstream out(ofh.c_str());
        out << "CHR\tSTART\tEND\tLENGTH\tOVERLAP\tREADNAME\tSTRAND\tSV\tTYPE" << endl;
	string line; while(getline(fin,line,'\n'))
	{
		vector<string> pos = split(line,'\t');
		if(pos.size()<4){ cerr << "ERROR: Malformed BED record " << line << endl; } 
		int32_t start = atoi(pos[1].c_str()); int32_t end = atoi(pos[2].c_str()); 	
		string CHR=pos[0]; string TYPE=pos[3];
		if(TYPE.find("DEL")){ continue; }
		if(chrFlag==false && !CHR.compare(0,3,"chr")) { CHR.erase(0,3); }
		if(chrom.count(CHR)==0) { cerr << "ERROR: " << CHR << " not found in BAM file reference" << endl; return 1; }
		if(!bam.LocateIndex()) { cerr << "ERROR: Cannot find index for BAM" << endl; return 1; }
		if(!bam.SetRegion(chrom[CHR],start-5001,chrom[CHR],end+5001)) { cerr << "ERROR: Is your BAM file indexed?" << endl; return 1; } 
		while (bam.GetNextAlignmentCore(al))
		{
			if(al.MapQuality < Q || 
				al.IsDuplicate()==true || 
				al.IsFailedQC()==true || 
				al.IsMapped()==false )
				{ continue; } 	
			char ori='+'; if(al.IsReverseStrand()) { ori='-'; }
			uint16_t len=0; 
			vector<CigarOp> cigar = al.CigarData;
			al.BuildCharData();
			for(vector<CigarOp>::iterator it=cigar.begin(); it != cigar.end(); ++it){
                		if (it->Type == 'D' && it->Length > 19) {
                			int32_t s1 = al.Position+len; 
					int32_t e1 = al.Position+len+it->Length; 
					s1++;  
					double ovr = overlap(start,end,s1,e1);
					if (ovr >= OVR) { 
					out << CHR << '\t' << s1 << '\t' << e1 << '\t' << e1-s1+1 << '\t' << ovr << '\t' << al.Name << '\t'<<  ori << '\t' << CHR << ':' << start << '-' << end << '\t' << TYPE << endl; 
					}
				}
                		if(it->Type == 'D' || it->Type == 'M' || it->Type == '=' || it->Type == 'X') { len+=it->Length; }
			}	
		
		}
	}
	bam.Close(); fin.close(); out.close(); 
	return 0;
}
