#include "manage.hpp"
#include <iostream>
#include <fstream>
using namespace std;
extern vector<double> distance_chain_list;
int main(int argc, char *argv[]){
 Run(argv[1]);
 ofstream distance_chain_file;
 distance_chain_file.open("distance_chain.txt");
 for(int i=0;i<=distance_chain_list.size();i++)
 {distance_chain_file<<distance_chain_list[i]<<'\n';}
 distance_chain_file.close();
 return 0;
}