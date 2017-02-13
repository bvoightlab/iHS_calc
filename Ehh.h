#include <vector>
#include <string.h>
#include <iostream>
#include <map>
#include <math.h>

using namespace std;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};


class Ehh {
  
  
 public:
  
  vector<float> ehh0;
  vector<float> ehh1;

  vector<int>  mhap0;
  vector<int>  mhap1;
   
  // statistic computed (interested)

  // uncorrected rho
  float rho0;
  float rho1;  
  
  // corrected rho
  float c_rho0;
  float c_rho1;

  // modified by BFV 07.20.12: retain exact EHH when dropping below THRESH
  float t_ehh0;
  float t_ehh1;

  // mod by BFV 07.20.12: also keep the index
  int t_index_ehh0;
  int t_index_ehh1;

  // longest common hapolotype length in 0 and 1
  float hap_len0;
  float hap_len1;
  
  // integral based on hap_len0 and hap_len1
  float int_ps0;
  float int_ps1;
  
  // the integral for 0 and 1
  float int_s0;
  float int_s1;
  
  // gap stats, snp density
  
  float max_ggap;       // maximum gap length in integrated region 
  
  int   total_marker;   // number of markers in integrated region 
  int   gap_count;      // number of gaps over threshold 
  float total_int_dist; // total integrated genetic distance

  
  Ehh(){
    thresh = 0.25;
    hapcount_thresh = 4;
  }
  
  void load_data(const vector<char *> &data, const vector<float> &map, const vector<float> &gap, vector<string> &warnings, const string &dir);

  void set_T_thresh(float value){
    thresh = value;
  }

  
 private:
  
  vector<char *> whole_data;  
  vector<float> gmap;
  vector<float>  gapv;  //correction factors for gaps in gmap

  vector<char *> data0;
  vector<char *> data1;
  
  float thresh;
  int   hapcount_thresh;


  void compute_ehh(vector<string> &warnings, const string &dir);
  float compute_ehh_score(int n0,int n1,int len, const vector<char *> &data,int &index);

  
};
