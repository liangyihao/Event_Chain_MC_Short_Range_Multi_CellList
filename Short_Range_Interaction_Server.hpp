//deal with short-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#ifndef SHORT_RANGE_INTERACTION_SERVER
#define SHORT_RANGE_INTERACTION_SERVER

#include <vector>
#include <iostream>
#include <cmath>
#include "public.hpp"
#include "CellList.hpp"
#include "Short_Range_potentials_Basic.hpp"

using namespace std;
extern vector<Parameter_List> Parameter_List_For_Short_Range_Interaction;
class Short_Range_Interaction_Between_Types{
	private:
	double(*Event_Time_Generator)(double4,double4,int,double*,double);
	double*Params;//Params[0]:Lx, Params[1]:Ly, Params[2]:Lz
	double SHORT_INTERACTION_RANGE;
	int type_id_1,type_id_2;
	vector<double4>*X_type_1;
	vector<double4>*X_type_2;

	/*For each type, create Cell List seperately.
	Sometimes, if the number of beads in such type is small, 
	users can choose turn-off the corresponding cell list.
	*/
	CellList* Cell_List_Pointer_1;
	CellList* Cell_List_Pointer_2;
	bool Using_CellList_1,Using_CellList_2;
	public:
	Short_Range_Interaction_Between_Types(int type_id_1, int type_id_2, vector<double4>*X_type_1, vector<double4>*X_type_2, double(*Event_Time_Generator)(double4,double4,int,double*,double), double*Params, double SHORT_INTERACTION_RANGE,bool Using_CellList_1,bool Using_CellList_2);
	void Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void Update(int2 ids, double4 New_X);
	
	double (*get_Gen())(double4,double4,int,double*,double){
		return Event_Time_Generator; 
	}
	double get_cutoff_r(){
		return SHORT_INTERACTION_RANGE;
	}
	int2 Interacting_Types(){
		int2 temp;
		temp.x=type_id_1;
		temp.y=type_id_2;
		return temp;
	}
	void Reconstruction(int axis,double New_Length){
		if(abs(axis)==1)Params[0]=New_Length;
		if(abs(axis)==2)Params[1]=New_Length;
		if(abs(axis)==3)Params[2]=New_Length;
		if(Using_CellList_1)Cell_List_Pointer_1->Reconstruct(axis,New_Length);
		if(Using_CellList_2)Cell_List_Pointer_2->Reconstruct(axis,New_Length);
	}
	int2 get_types(){
		int2 temp;
		temp.x=type_id_1;
		temp.y=type_id_2;
		return temp;
	}
	void Get_Left_Wrap_Event_Hard_Sphere(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void check_overlap_for_Hard_Core(double Lx,double Ly,double Lz){
		if(Event_Time_Generator!=Event_Time_Hard_Sphere)return;
		vector<int>bead_ids;
		for(int bead_id1=0;bead_id1<X_type_1->size();bead_id1++){
			if(Using_CellList_2){
				Cell_List_Pointer_2->particles_in_range((*X_type_1)[bead_id1],SHORT_INTERACTION_RANGE,bead_ids);
			}else{
				double dx,dy,dz,dr2;
				for(int g=0;g<X_type_2->size();g++){
					dx=(*X_type_1)[bead_id1].x-(*X_type_2)[g].x;
					dy=(*X_type_1)[bead_id1].y-(*X_type_2)[g].y;
					dz=(*X_type_1)[bead_id1].z-(*X_type_2)[g].z;
					if(dx<-Lx/2)dx+=Lx;if(dx>Lx/2)dx-=Lx;
					if(dy<-Ly/2)dy+=Ly;if(dy>Ly/2)dy-=Ly;
					if(dz<-Lz/2)dz+=Lz;if(dz>Lz/2)dz-=Lz;
					dr2=dx*dx+dy*dy+dz*dz;
					if(dr2<(1.0-100*EPSILON)*SHORT_INTERACTION_RANGE*SHORT_INTERACTION_RANGE)bead_ids.push_back(g);
				}
			}
			if(type_id_1==type_id_2){//remove itself from bead_ids list
				for(int l=0;l<bead_ids.size();l++)if(bead_ids[l]==bead_id1){
					bead_ids[l]=bead_ids[bead_ids.size()-1];
					bead_ids.pop_back();
				}
			}
			if(bead_ids.size()>0){
				for(int l=0;l<bead_ids.size();l++){
					cout<<"Error, Hard Sphere Overlap, they are"<<endl;
					cout<<type_id_1<<' '<<bead_id1<<" Position:("<<(*X_type_1)[bead_id1].x<<','<<(*X_type_1)[bead_id1].y<<','<<(*X_type_1)[bead_id1].z<<')'<<endl;
					cout<<type_id_2<<' '<<bead_ids[l]<< " Position:("<<(*X_type_2)[bead_ids[l]].x<<','<<(*X_type_2)[bead_ids[l]].y<<','<<(*X_type_2)[bead_ids[l]].z<<')'<<endl;
				}
				exit(0);
			}
		}
	}
};
#endif
