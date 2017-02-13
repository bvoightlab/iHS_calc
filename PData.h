#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <string.h>
#include "Ehh.h"

using namespace std;


class PData {
  
 public:

  PData(){
    rsize = 2500000;
    gap_psize = 100000; //modified BFV 10.28.08: changes from 20000 to 100000 for relaxed gap correction due to WGAS SNP density.
  }
  
  void load_info(char *infofile);
  int  load_data(char *datafile);  
  void compute_data(int index);
  
  int snp_num;


  vector<string> snp_list;

 private:
  
  vector<char *> data;


  vector<int>   phy_map;
  vector<float>  gen_map;
  vector<int>    sw_list;

  Ehh ehh;

  
  int rsize;
  int gap_psize;
  
  
};
