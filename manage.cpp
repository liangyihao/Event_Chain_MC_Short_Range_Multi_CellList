/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/
#include "manage.hpp"
#include "basic.hpp"
#include "Random_Number.hpp"
#include "Input_File_Parser.hpp"
#include "dcd_writer.hpp"
#include <iostream>
/*Used to manage systematic variables*/
double Lx=10,Ly=10,Lz=10;//0<=x<Lx...
vector<Instruction>Instruction_list;
int loop_times;


int Create_Type() {//Create Type, return it's id
    Bead_Type new_type;
    new_type.index=Types.size();
    new_type.Interactions_with_Beads.clear();
    new_type.Interactions_with_Types.clear();
    new_type.X.clear();
    Types.push_back(new_type);
    return new_type.index;
}

int Create_Bead(int type_id,double4 x){//Create Bead, return its secondary id. If failed, return -1.
    if(type_id>=Types.size()) {
        cout<<"Error Creating Bead: type-id incorrect"<<endl;
        return -1;
    }
    
    if((x.x>Lx)||(x.x<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.y>Ly)||(x.y<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    if((x.z>Lz)||(x.z<0)){
        cout<<"Error Creating Bead: position out of boundary"<<endl;
        return -1;
    }

    Types[type_id].X.push_back(x);
    vector<int3>temp;
    temp.clear();
    Types[type_id].Interactions_with_Beads.push_back(temp);
    return Types[type_id].X.size()-1;
}


void Run(char*InputFileName){
    rand_init(1);
    Input_File_Parser(InputFileName);
    Hard_Repulsion_Checker();

    Output_DCD_init(InputFileName);
    for(int l=0;l<loop_times;l++){
        for(int k=0;k<Instruction_list.size();k++){
            if(Instruction_list[k].Command==0){//Do ECMC
                Monte_Carlo(Instruction_list[k].Double_Para[0],Instruction_list[k].Int_Para[0]);
            }else if(Instruction_list[k].Command==1){//Do output
                if(l%Instruction_list[k].Int_Para[1]==0)cout<<l<<endl;
                if(l<Instruction_list[k].Int_Para[0])continue;
                if(l%Instruction_list[k].Int_Para[1]==0){
                    Output_DCD();
                    next_input_file_writer(InputFileName);
                }
            }
        }
    }
    Output_DCD_Close();
    next_input_file_writer(InputFileName);
}
