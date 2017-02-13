#include "PData.h"

int PData::load_data(char *datafile){
  
  
  FILE *fp = fopen(datafile,"r");
  if(fp==0)
    return -1;
    
  vector<char> line;
  
  while(1){
    char in = fgetc(fp); 
    if(isalnum(in)){
      line.push_back(in);
      continue;
    }
    if(in=='\n'||in=='\r'||in==EOF){   
      // convert line to single char array
      if(line.size()>0){
	char *pl = new char[line.size()+1];
	memset(pl,0,line.size()+1);
	for(int i=0;i<line.size();i++){
	  pl[i]=line[i];	  
	}
	data.push_back(pl);	
      }
      line.clear();
    }

    if(in==EOF)
      break;
  }
  
  if(data.size()<1)
    return -2;

  
  snp_num = strlen(data[0]);
  
  return 1;
  
}


void PData::load_info(char *infofile){
  
  FILE *fp = fopen(infofile,"r");
  if(fp==0)
    return;
  char line[128];


  while(fgets(line,128,fp)){
    
    char *str = strtok(line," ");
    string marker(str);
  
    snp_list.push_back(marker);
    
    str = strtok(NULL," ");
    phy_map.push_back(atoi(str));
    
    str = strtok(NULL," ");
    gen_map.push_back(atof(str));
    
    str = strtok(NULL," ");
    char c0[3];
    strcpy(c0,str);
    
    str = strtok(NULL," ");
    if(str[strlen(str)-1]=='\n')
      str[strlen(str)-1]=0;
    
    if(strcmp(c0,str)==0)
      sw_list.push_back(0);
    else if(strcmp(str,"NA")==0||strcmp(str,"?")==0)
      sw_list.push_back(-1);
    else
      sw_list.push_back(1);
    
  }
  fclose(fp);
      
}

void PData::compute_data(int index){

  //modified BFV 10.22.08: add a warning vector to capture potentials in calculations.
  vector<string> warning_vec;
  
  if(sw_list[index]==-1) {
    return;
  }

  int ipos = phy_map[index];
  int lpos = ipos-rsize;
  int rpos = ipos+rsize;
  
  int lbound = -1;
  int rbound = -1;
  
  for(int i=0;i<phy_map.size();i++){
    if(lbound==-1 && phy_map[i]>=lpos)
      lbound=i;
    if(rbound==-1&& phy_map[i]>rpos){
      rbound=i-1;
      break;
    }
  }
  
  if(rbound==-1)
    rbound = phy_map.size()-1;
  
  
  //modified BFV 10.29.08: removed to make sure all SNPs are printed [this was dropping the last one in the info file]
  //if(index-lbound+1<0 || rbound-index+1<2)
  //  return;
 
  //cout << "X"; //flag for edge issues.
  
  //modified BFV 10.22.08: warning flag to vector. -- delete
  //warningflag = "OK";
  //warning_vec.push_back(warningflag); //flag for edge issues.
  
  //get left data
  vector<char *> left_data;
  vector<float>  left_gmap;
  vector<float> lgap;
 
  for(int i=0;i<data.size();i++){
    int len = index-lbound+2;
    char *ldata = new char[len];
    memset(ldata,0,len);
    
    if(sw_list[index]==1){
      if(data[i][index]=='1') //modified BFV 10.25.08: "1" is derived allele, "0" is ancestral by default.
	ldata[0]='1';
      else
	ldata[0]='0';
    }else
      ldata[0]=data[i][index];
      
    
    int count=1;
    for(int j=index-1;j>=lbound;j--){
      
      // check for gap
      if(i==0){
	float gaps = phy_map[j+1]-phy_map[j];
	if( gaps >gap_psize){
	  
	  if(gaps > 3*gap_psize)
	    lgap.push_back(0);
	  else
	    lgap.push_back((float)gap_psize/gaps);
	}	
	else 
	  lgap.push_back(1);
      }
      
      ldata[count++]=data[i][j];

    }
    
    left_data.push_back(ldata);
  }
  
  for(int i=index;i>=lbound;i--)
    left_gmap.push_back(gen_map[i]);
   
  //computation is done here
  ehh.load_data(left_data,left_gmap,lgap,warning_vec, "l");

  
  // get results
  float l0 = ehh.rho0;
  float l1 = ehh.rho1;
  
  float cl0 = ehh.c_rho0;
  float cl1 = ehh.c_rho1;
  
  float hap_l0 = ehh.hap_len0;
  float hap_l1 = ehh.hap_len1;
  
  float int_ps0 = ehh.int_ps0;
  float int_ps1 = ehh.int_ps1;
  
  float int_s0 =  ehh.int_s0;
  float int_s1 =  ehh.int_s1;
  
  float max_ggap = ehh.max_ggap;
  
  int   gap_count = ehh.gap_count;
  int   total_marker = ehh.total_marker;
  float total_int_dist = ehh.total_int_dist;
  
  //mod by BFV 07.20.12: get more stuff
  float t_l_ehh0 = ehh.t_ehh0;
  float t_l_ehh1 = ehh.t_ehh1;

  int t_index_l_ehh0 = ehh.t_index_ehh0;
  int t_index_l_ehh1 = ehh.t_index_ehh1;

  //clean up left_data
  left_gmap.clear();
  lgap.clear();
  for(int i=0;i<left_data.size();i++)
    delete[](left_data[i]);
  left_data.clear();
  
  if(l0==-1|| l1==-1) {
    //cout << "ZL";    
    warning_vec.push_back("Size_EHH_l");
    //return;
  }  

  //get right data
  vector<char *> right_data;
  vector<float>  right_gmap;
  vector<float>    rgap;
  float freq = 0;
  for(int i=0;i<data.size();i++){
  
    char *rdata = new char[rbound-index+2];
    memset(rdata,0,rbound-index+2);
        
    if(sw_list[index]==1){
      if(data[i][index]=='1') //modified 10.25.08 BFV "1" is derived allele by default.
	rdata[0]='1';
      else
	rdata[0]='0'; 
    }else{
      rdata[0]=data[i][index];
    }
    
    if(rdata[0]=='1') //Calculate derived allele frequency
      freq++;
    
    int count=1;
    for(int j=index+1;j<=rbound;j++){
      // check for gap
      if(i==0){
	float gaps = phy_map[j]-phy_map[j-1];
	if(gaps>gap_psize){
	  if(gaps>3*gap_psize)
	    rgap.push_back(0);
	  else
	    rgap.push_back((float)gap_psize/gaps);
	}
	else
	  rgap.push_back(1);
      }
      rdata[count++]=data[i][j];
    }
    
        
    right_data.push_back(rdata);
  }
  
  for(int i=index;i<=rbound;i++)
    right_gmap.push_back(gen_map[i]);
     
  //cout << "Y2";

  //computation is done here
  ehh.load_data(right_data,right_gmap,rgap,warning_vec, "r");
  
  // get results
  float r0 = ehh.rho0;
  float r1 = ehh.rho1;
  
  float cr0 = ehh.c_rho0;
  float cr1 = ehh.c_rho1;
  
  hap_l0 += ehh.hap_len0;
  hap_l1 += ehh.hap_len1;
  
  int_ps0 += ehh.int_ps0;
  int_ps1 += ehh.int_ps1;
  
  int_s0 +=  ehh.int_s0;
  int_s1 +=  ehh.int_s1;
  
  //mod by BFV 07.20.12: get more stuff
  float t_r_ehh0 = ehh.t_ehh0;
  float t_r_ehh1 = ehh.t_ehh1;

  int t_index_r_ehh0 = ehh.t_index_ehh0;
  int t_index_r_ehh1 = ehh.t_index_ehh1;

  if(ehh.max_ggap > max_ggap)
    max_ggap = ehh.max_ggap;
    
  gap_count += ehh.gap_count;
  total_marker += ehh.total_marker;
  total_int_dist += ehh.total_int_dist;
  
  freq = freq/(data.size());
  
  //clean up right_data
  right_gmap.clear();
  for(int i=0;i<right_data.size();i++)
    delete[](right_data[i]);
  right_data.clear();
  
  if(r0==-1|| r1==-1) {
    //cout << "ZA";
    warning_vec.push_back("Size_EHH_r");
    //return;
  }

  //BFV 11.26.2014: obtain SNP ids where EHH drops below threshold
  string snpid_l_ehh0;
  string snpid_r_ehh0;
  string snpid_l_ehh1;
  string snpid_r_ehh1;

  int pos_l_ehh0;
  int pos_r_ehh0;
  int pos_l_ehh1;
  int pos_r_ehh1;

  //process 
  if (t_index_r_ehh0 == -1) { //edge case
    snpid_r_ehh0 = "NA";
    pos_r_ehh0 = -1;
  } else {
    snpid_r_ehh0 = snp_list[(index + t_index_r_ehh0)];
    pos_r_ehh0 = phy_map[(index + t_index_r_ehh0)];
  }

  if (t_index_l_ehh0 == -1) { //edge case
    snpid_l_ehh0 = "NA";
    pos_l_ehh0 = -1;
  } else {
    snpid_l_ehh0 = snp_list[(index - t_index_l_ehh0)];
    pos_l_ehh0 = phy_map[(index - t_index_l_ehh0)];
  }  

  if (t_index_r_ehh1 == -1) { //edge case
    snpid_r_ehh1 = "NA";
    pos_r_ehh1 = -1;
  } else {
    snpid_r_ehh1 = snp_list[(index + t_index_r_ehh1)];
    pos_r_ehh1 = phy_map[(index + t_index_r_ehh1)];
  }

  if (t_index_l_ehh1 == -1) { //edge case
    snpid_l_ehh1 = "NA";
    pos_l_ehh1 = -1;
  } else {
    snpid_l_ehh1 = snp_list[(index - t_index_l_ehh1)];
    pos_l_ehh1 = phy_map[(index - t_index_l_ehh1)];
  }

  //cout << t_index_l_ehh0 << " " << t_index_r_ehh0 << " " << t_index_l_ehh1 << " " << t_index_r_ehh1 << endl;
  //cout << snpid_l_ehh0 << " " << snpid_r_ehh0 << " " << snpid_l_ehh1 << " " << snpid_r_ehh1 << endl;
  //std::exit(0);

  //printf("%11s %10d %.3f   %7.3f    %7.2f %7.2f %7.2f %7.3f   %7.2f %7.2f %7.2f %7.3f   %5.2f %5.2f %d\n",snp_list[index].c_str(),phy_map[index],freq, log((cr0+cl0)/(cr1+cl1)),int_s0,int_s1,int_s0-int_s1,log(int_s0/int_s1),int_ps0,int_ps1,int_ps0-int_ps1,log(int_ps0/int_ps1),total_marker/total_int_dist,max_ggap,gap_count);  

  //modified BFV 10.22.08: reduced the net output of the analysis.
  float ehh_thresh = log((cr0+cl0)/(cr1+cl1)); //this is a statistic with EHH threshold fixed at (say) 0.25; rho is calculated at that value and log(EHH0_0.25/EHH1_0.25) is tabulated
  float ihs = log(int_s0/int_s1); //this appears to be iHS
  //float ihs_three = log(int_ps0/int_ps1); //this appears to be iHS on the most common haplotype

  if (warning_vec.size() > 0) {
    printf("%-11s %-10d %.4f  %8s  %7.2f %7.2f %4d  %7.4f %7.4f %7.4f %7.4f %10.2f %10.2f %10.2f %10.2f  %-11s %-11s %-11s %-11s %-10d %-10d %-10d %-10d  ",snp_list[index].c_str(),phy_map[index],freq,"NA",total_marker/total_int_dist,max_ggap,gap_count,t_l_ehh0,t_l_ehh1,t_r_ehh0,t_r_ehh1,l0,r0,l1,r1,snpid_l_ehh0.c_str(),snpid_r_ehh0.c_str(),snpid_l_ehh1.c_str(),snpid_r_ehh1.c_str(),pos_l_ehh0,pos_r_ehh0,pos_l_ehh1,pos_r_ehh1);
  } else {
    if (ihs >= 0) {
      printf("%-11s %-10d %.4f    % 4.3f  %7.2f %7.2f %4d  %7.4f %7.4f %7.4f %7.4f %10.2f %10.2f %10.2f %10.2f  %-11s %-11s %-11s %-11s %-10d %-10d %-10d %-10d  ",snp_list[index].c_str(),phy_map[index],freq,ihs,total_marker/total_int_dist,max_ggap,gap_count,t_l_ehh0,t_l_ehh1,t_r_ehh0,t_r_ehh1,l0,r0,l1,r1,snpid_l_ehh0.c_str(),snpid_r_ehh0.c_str(),snpid_l_ehh1.c_str(),snpid_r_ehh1.c_str(),pos_l_ehh0,pos_r_ehh0,pos_l_ehh1,pos_r_ehh1);
    } else {
      printf("%-11s %-10d %.4f    %4.3f  %7.2f %7.2f %4d  %7.4f %7.4f %7.4f %7.4f %10.2f %10.2f %10.2f %10.2f  %-11s %-11s %-11s %-11s %-10d %-10d %-10d %-10d  ",snp_list[index].c_str(),phy_map[index],freq,ihs,total_marker/total_int_dist,max_ggap,gap_count,t_l_ehh0,t_l_ehh1,t_r_ehh0,t_r_ehh1,l0,r0,l1,r1,snpid_l_ehh0.c_str(),snpid_r_ehh0.c_str(),snpid_l_ehh1.c_str(),snpid_r_ehh1.c_str(),pos_l_ehh0,pos_r_ehh0,pos_l_ehh1,pos_r_ehh1);
    }
  }

  if (warning_vec.size() == 0) {
    cout << "OK" << endl;
  } else {
    for (int k=0; k<warning_vec.size(); k++) {
      if (k != 0) {
	cout << ",";
      }
      cout << warning_vec[k];
    }
    cout << endl;
  }

}
  
int main(int argc, char **argv){
  
  PData pd;
  int start_loc;
  int end_loc;

  if (argc == 1) {
    cout << "usage: %>iHS_calc [compilecheck] myfile.info myfile.data [start_loc end_loc]" << endl;
    cout << "Version 1.4 [from BFV 11.26.14]" << endl;
    return(0);
  } else if (argc == 2) {
    cout << "OK!\n";
    return(0);
  } else if (argc == 3) {
    pd.load_info(argv[1]);
    pd.load_data(argv[2]);
    start_loc = 0;
    end_loc = pd.snp_num-1;
  } else if (argc == 5) {
    pd.load_info(argv[1]);
    pd.load_data(argv[2]);
    start_loc = atoi(argv[3]);
    end_loc = atoi(argv[4]);

    //perish if start loc is larger than size.
    if (start_loc > pd.snp_num) {
      cerr << "ERROR: Start location is larger than SNP size. Exiting." << endl;
      return(0);
    }

    //set size to last SNP if length is larger than end
    if (end_loc >= pd.snp_num) {
      cerr << "WARNING: End location is larger than SNP size. Setting to last SNP." << endl;      
      end_loc = pd.snp_num-1;
    }
  }
  
  //cout << start_loc << endl;
  //cout << end_loc << endl;
  //return(0);

  //modified BFV 10.22.08: added a header line to the output
  printf("%-11s %-10s %5s %9s  %7s %7s %4s  %7s %7s %7s %7s %10s %10s %10s %10s  %-11s %-11s %-11s %-11s %-10s %-10s %-10s %-10s  %8s\n","SNP","POS","FREQ_1","unstd_iHS","density","max_gap","ngap","L-EHH0","L-EHH1","R-EHH0","R-EHH1","L-EHH0_RHO","R-EHH0_RHO","L-EHH1_RHO","R-EHH1_RHO","L-EHH0_SNP","R-EHH0_SNP","L-EHH1_SNP","R-EHH1_SNP","L-EHH0_POS","R-EHH0_POS","L-EHH1_POS","R-EHH1_POS","Warnings");

  for(int i=start_loc; i<=end_loc; i++) 
    pd.compute_data(i);
  
  /*
  pd.load_info(argv[1]);
  pd.load_data(argv[2]);
  for(int i=0;i<pd.snp_num;i++){ 
    if(strcmp(pd.snp_list[i].c_str(),argv[3])==0){
      pd.compute_data(i);
      return 1;
    }
  }
  */

}
