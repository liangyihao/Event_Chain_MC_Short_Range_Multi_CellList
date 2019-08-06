#include "manage.hpp"
#include "basic.hpp"
#include "Random_Number.hpp"
#include "Input_File_Parser.hpp"
#include "dcd_writer.hpp"
#include "Short_Range_potentials_Basic.hpp"
#include <iostream>
#include <cmath>

void Monte_Carlo_No_Wrap(int2 First_Active_Bead,int axis) {
    double Clock=0;
    double time,exe_time;
    double temp;
    Active_Bead=First_Active_Bead;

    if((abs(axis)>3)||(axis==0)){cout<<"Error, axis id should be -3,-2,-1,1,2,3"<<endl;return;}
    int2 id_next_active_bead;
    bool go_ahead=true;
    while(go_ahead) {
        Get_Event(time,id_next_active_bead,axis);

        if(abs(axis)==1)exe_time=min(Lx,time);//axis==1: on +x direction, axis==-1: on -x direction
        if(abs(axis)==2)exe_time=min(Ly,time);//axis==2: on +y direction, axis==-2: on -y direction
        if(abs(axis)==3)exe_time=min(Lz,time);//axis==3: on +z direction, axis==-3: on -z direction
        
        double4 New_X;
        New_X=Types[Active_Bead.x].X[Active_Bead.y];
        if(abs(axis)==1){
            New_X.x+=exe_time*(axis/abs(axis));
            if(New_X.x>Lx){New_X.x=Lx;go_ahead=false;}
            if(New_X.x<0 ){New_X.x= 0;go_ahead=false;}
        }
        if(abs(axis)==2){
            New_X.y+=exe_time*(axis/abs(axis));
            if(New_X.y>Ly){New_X.y=Ly;go_ahead=false;}
            if(New_X.y<0 ){New_X.y= 0;go_ahead=false;}
        }
        if(abs(axis)==3){
            New_X.z+=exe_time*(axis/abs(axis));
            if(New_X.z>Lz){New_X.z=Lz;go_ahead=false;}
            if(New_X.z<0 ){New_X.z= 0;go_ahead=false;}
        }
        for(int k=0;k<Types[Active_Bead.x].Interactions_with_Types.size();k++){
            Short_Range_Interaction_Between_Types_List[Types[Active_Bead.x].Interactions_with_Types[k]]->Update(Active_Bead,New_X);
        }
        Types[Active_Bead.x].X[Active_Bead.y]=New_X;
        Clock+=exe_time;
        Active_Bead=id_next_active_bead;
    }
}

bool Compress_Sys(int axis, double target_length, double init_active_layer_thickness){//if length successfully reach target length, return true;else return false.
    cout<<"Old size:"<<Lx<<','<<Ly<<','<<Lz<<endl;
    for(int l=0;l<Short_Range_Interaction_Between_Types_List.size();l++)//For debug
        Short_Range_Interaction_Between_Types_List[l]->check_overlap_for_Hard_Core();


    axis=abs(axis);
    double L_axis;
    if(axis==1)L_axis=Lx;
    if(axis==2)L_axis=Ly;
    if(axis==3)L_axis=Lz;
    if(L_axis-target_length<EPSILON*target_length)return true;
    cout<<"Try compressing system in ";
    if(axis==1)cout<<"Lx ";if(axis==2)cout<<"Ly ";if(axis==3)cout<<"Lz ";
    cout<<"Direction"<<endl;
    int2 First_Active_Bead;
    for(int type_id=0;type_id<Types.size();type_id++)
        for(int bead_id=0;bead_id<Types[type_id].X.size();bead_id++){
            double h;
            if(axis==1)h=Lx-Types[type_id].X[bead_id].x;
            if(axis==2)h=Ly-Types[type_id].X[bead_id].y;
            if(axis==3)h=Lz-Types[type_id].X[bead_id].z;

            if(h>init_active_layer_thickness)continue;
            First_Active_Bead.x=type_id;
            First_Active_Bead.y=bead_id;
            Monte_Carlo_No_Wrap(First_Active_Bead,-axis);
        }
    /////////////////////////Now determine the new L-axis
    double min_collision=max(Lx,max(Ly,Lz));
    TwoBody_Event Event;
    int2 AB,min_collision_A,min_collision_B;
    double4 X_AB,X_TG;
    for(int f=0;f<Short_Range_Interaction_Between_Types_List.size();f++){
        if(Short_Range_Interaction_Between_Types_List[f]->get_Gen()!=Event_Time_Hard_Sphere)continue;
        int2 types;
        int type_id1,type_id2;
        types=Short_Range_Interaction_Between_Types_List[f]->get_types();
        type_id1=types.x;
        type_id2=types.y;

        for(int bead_id=0;bead_id<Types[type_id1].X.size();bead_id++){
            AB.x=type_id1;
            AB.y=bead_id;
            X_AB=Types[type_id1].X[bead_id];
            Event.Event_Time=2*L_axis;
            Event.Target_Bead=AB;
            Short_Range_Interaction_Between_Types_List[f]->Get_Left_Wrap_Event_Hard_Sphere(Event,AB,X_AB,axis);
            X_TG=Types[Event.Target_Bead.x].X[Event.Target_Bead.y];
            double dx_axis;
            if(axis==1)dx_axis=X_TG.x-X_AB.x;
            if(axis==2)dx_axis=X_TG.y-X_AB.y;
            if(axis==3)dx_axis=X_TG.z-X_AB.z;
            if(dx_axis>0)continue;
            if(Event.Event_Time<min_collision){
                min_collision=Event.Event_Time;
                min_collision_A=AB;
                min_collision_B=Event.Target_Bead;
            }
        }
        
        for(int bead_id=0;bead_id<Types[type_id2].X.size();bead_id++){
            AB.x=type_id2;
            AB.y=bead_id;
            X_AB=Types[type_id2].X[bead_id];
            Event.Event_Time=2*L_axis;
            Event.Target_Bead=AB;
            Short_Range_Interaction_Between_Types_List[f]->Get_Left_Wrap_Event_Hard_Sphere(Event,AB,X_AB,axis);
            X_TG=Types[Event.Target_Bead.x].X[Event.Target_Bead.y];
            double dx_axis;
            if(axis==1)dx_axis=X_TG.x-X_AB.x;
            if(axis==2)dx_axis=X_TG.y-X_AB.y;
            if(axis==3)dx_axis=X_TG.z-X_AB.z;
            if(dx_axis>0)continue;
            if(Event.Event_Time<min_collision){
                min_collision=Event.Event_Time;
                min_collision_A=AB;
                min_collision_B=Event.Target_Bead;
            }
        }
    }
    cout<<"Min_Collision: "<<min_collision<<" between beads:("<<min_collision_A.x<<','<<min_collision_A.y<<") and ("<<min_collision_B.x<<","<<min_collision_B.y<<")"<<endl;
    ////////////////////////Do compress//////////////////////////
    double New_Length;
    New_Length=L_axis-min_collision;
    if(New_Length<target_length)New_Length=target_length;
    for(int l=0;l<Short_Range_Interaction_Between_Types_List.size();l++)
        Short_Range_Interaction_Between_Types_List[l]->Reconstruction(axis,New_Length);
    
    for(int l=0;l<Param_Lists_For_Bonds.size();l++){
        if(abs(axis)==1)Param_Lists_For_Bonds[l].data[0]=New_Length;
        if(abs(axis)==2)Param_Lists_For_Bonds[l].data[1]=New_Length;
        if(abs(axis)==3)Param_Lists_For_Bonds[l].data[2]=New_Length;
    }
    if(axis==1)Lx=New_Length;
    if(axis==2)Ly=New_Length;
    if(axis==3)Lz=New_Length;

    cout<<"New size:"<<New_Length<<endl;
    for(int l=0;l<Short_Range_Interaction_Between_Types_List.size();l++)//For debug
        Short_Range_Interaction_Between_Types_List[l]->check_overlap_for_Hard_Core();

    if(abs(New_Length-target_length)<EPSILON*target_length){
        cout<<"Target length reached"<<endl;
        cout<<"Current Length:"<<New_Length<<endl;
        return true;
    }
    cout<<"Try finished, current length "<<New_Length<<endl;
    return false;
}