/*
Author: Yihao Liang
liangyihaosjtu@gmail.com
This code is for Event Chain Monte Carlo for pairwise interacting many body system
*/

#ifndef BASIC_ECMC
#define BASIC_ECMC

#include <vector>
#include "public.hpp"
#include "Short_Range_Interaction_Server.hpp"

using namespace std;

extern vector<Bead_Type> Types;
extern int2 Active_Bead;
extern vector<double(*)(double4,double4,int,double*,double)> Event_Time_Generator_List_For_Bonds;
extern vector<Parameter_List> Param_Lists_For_Bonds;
extern vector<Short_Range_Interaction_Between_Types*>Short_Range_Interaction_Between_Types_List;
void Get_Event(double&time, int2&id_next_active_bead, int axis);
void Monte_Carlo(double Stop_Clock,int axis);
#endif
