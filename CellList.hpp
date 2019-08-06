//deal with short-range interactions
//range:[0,Lx],[0,Ly],[0,Lz]
#ifndef CELL_LIST
#define CELL_LIST

#include <vector>
#include <iostream>
#include "public.hpp"
using namespace std;

typedef struct CellStruct{
	int type_id;//type id of particles in this cell
	vector<int>particle_list;//secondary id of particles in this cell
	CellStruct(int type_id){
		this->type_id=type_id;
		particle_list.clear();
	}

	void Insert(int bead_id) {
		particle_list.push_back(bead_id);
		return;
}
	void Delete(int bead_id) {
		for(int k=0;k<particle_list.size();k++)
			if(particle_list[k]==bead_id) {
				particle_list[k]=particle_list[particle_list.size()-1];
				particle_list.pop_back();
				return;
			}
		cout<<"Deletion failed. The particle ("<<type_id<<","<<bead_id<<") not found."<<endl;
		return;
}

}Cell;

class CellList{
private:
//One cell list only store 1 type of bead and serve only one short range interaction instance
	int type_id;
	double (*Gen)(double4,double4,int,double*,double);
	double *Params;

	vector<double4>*X_pointer;
	double Lx,Ly,Lz;
	double dL;//size of each cell
	int NC_x,NC_y,NC_z;//Number of cells per dimension
	vector<int3>InWhichCell;//InWhichCell[bead_id] stores the in which cell the bead (type_id,bead_is) is
	vector<vector<vector<Cell> > >Cells;
	/*
	CellList[i][j][k] access the cell in 
	i*Lx<= x < (i+1)*Lx
	j*Ly<= y < (j+1)*Ly
	k*Lz<= z < (k+1)*Lz
	*/
	void Event_with_Cell(TwoBody_Event&Event, int3 Cell_IJK, int2 Active_Bead, double4 X_Active_Bead, int axis);
	int3 In_Which_Cell(int2 ids,double4 X_ids) {
		if(ids.x==type_id){
			return InWhichCell[ids.y];
		}else{
			int3 temp;
			temp.x=X_ids.x/dL;
			temp.y=X_ids.y/dL;
			temp.z=X_ids.z/dL;
			
			if(temp.x==NC_x)temp.x=NC_x-1;
			if(temp.y==NC_y)temp.y=NC_y-1;
			if(temp.z==NC_z)temp.z=NC_z-1;

			return temp;
		}
	}
	void Event_with_Cell_Hard_Wrap(TwoBody_Event&Event, int3 Cell_IJK, int2 Active_Bead, double4 X_Active_Bead, int axis);
public:
	CellList(double Lx, double Ly, double Lz, double dL, int type_id, vector<double4>*X_pointer, double (*Gen)(double4,double4,int,double*,double),double *Params);
	void Get_Event(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void Update(int2 ids, double4 New_X);
	void print();
	void Reconstruct(int axis, double New_Length);
	void Get_Left_Wrap_Event_Hard_Sphere(TwoBody_Event&Event, int2 Active_Bead, double4 X_Active_Bead, int axis);
	void particles_in_range(double4 X, double range, vector<int>&bead_ids){
			int3 IWC;
			IWC.x=X.x/dL;
			IWC.y=X.y/dL;
			IWC.z=X.z/dL;
			
			if(IWC.x==NC_x)IWC.x=NC_x-1;
			if(IWC.y==NC_y)IWC.y=NC_y-1;
			if(IWC.z==NC_z)IWC.z=NC_z-1;

			int dI;
			dI=range/dL+1;
			bead_ids.clear();
			double dx,dy,dz,dr2;
			for(int i=IWC.x-dI;i<IWC.x+dI;i++)
				for(int j=IWC.y-dI;j<IWC.y+dI;j++)
					for(int k=IWC.z-dI;k<IWC.z+dI;k++){
						int I,J,K;
						I=i;J=j;K=k;
						while(I<0)I+=NC_x;while(I>=NC_x)I-=NC_x;
						while(J<0)J+=NC_y;while(J>=NC_y)J-=NC_y;
						while(K<0)K+=NC_z;while(K>=NC_z)K-=NC_z;
						for(int f=0;f<Cells[I][J][K].particle_list.size();f++){
							int g;
							g=Cells[I][J][K].particle_list[f];
							dx=X.x-(*X_pointer)[g].x;
							dy=X.y-(*X_pointer)[g].y;
							dz=X.z-(*X_pointer)[g].z;
							if(dx<-Lx/2)dx+=Lx;if(dx>Lx/2)dx-=Lx;
							if(dy<-Ly/2)dy+=Ly;if(dy>Ly/2)dy-=Ly;
							if(dz<-Lz/2)dz+=Lz;if(dz>Lz/2)dz-=Lz;

							dr2=dx*dx+dy*dy+dz*dz;
							if(dr2<(1.0-100*EPSILON)*range*range)bead_ids.push_back(g);
						}
					}
	}
};
#endif
