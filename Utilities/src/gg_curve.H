/****************************************************************/
/*    GRIDGEN  -->> .rea                                        */
/*    Author: Igor Pivkin, Brown University, RI                 */
/*    email: piv@cfm.brown.edu                                  */
/****************************************************************/

#ifndef _gg2rea_H_
#define _gg2rea_H_


#define EPSILON 0.0000001
#define MAXBC 5
#define MAXBCLINE 50
#define NONE -99999

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

using namespace std;

class Mesh;

typedef set<int> Int_Set;
typedef Int_Set::iterator Int_Set_Iter;

class Vertex{
private:
public:
  double x,y,z;

  double nx, ny, nz; /*normal vector contributed by all the faces sharing this vertex */
  double nf;  /* number of faces sharing this vertex or weigth */

  Int_Set elmts_vertisin;//elements vertex is in
  Vertex(){};
  Vertex(const double x_, const double y_, const double z_){
    x = x_; y = y_; z = z_;
    nx = 0.0; ny = 0.0; nz = 0.0;
    nf = 0.0;
  };
  ~Vertex(){
    //need to fill in destructor
  };
  inline void Add_elmts_vertisin(const int eid_){
    elmts_vertisin.insert(eid_);
  };
  inline int Get_number_of_elmts_vertisin(){
    return elmts_vertisin.size();
  };
  inline void Clear_elmts_vertisin(){
    elmts_vertisin.clear();
  };
};//class Vertex

enum {
  HEX = 0,
  TET = 1,
  PRI = 2,
  PYR = 3
};


class Element{
protected:
public:
  int type;
  int myid;
  int vert_num;
  int* vert_ids;
  int face_num;
  Mesh* mesh; //pointer to the mesh
  int* conn_element;
  int* conn_face;
  Element(){ mesh = NULL; vert_ids = NULL; conn_element = NULL; conn_face = NULL;};
  ~Element(){};//fill in destructor
  void Set_myid(const int myid_);
  int Get_myid(){ return myid; };
  virtual void Orient(){
    fprintf(stderr,"virtual Element::orientate() is being called\n");
  };
  void Get_mesh(Mesh* mesh_);
  int Get_number_of_vertices();
  int Get_vert_id(const int vid_);
  int Get_number_of_faces();
  virtual int Find_match_face(const int v0,const int v1,const int v2, const int v3){
    fprintf(stderr,"virtual Element::Match_face(...) is being called\n");
    return NONE;
  };
  virtual void Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3){
    fprintf(stderr,"virtual Element::Get_face_ids(...) is being called\n");
  };
  void Set_vert_id(const int id_, const int vid_);
};//class Element

class Tet: public Element{
private:
public:
  Tet() : Element(){
    type = TET;
    vert_num = 4;
    vert_ids = new int[vert_num];
    face_num = 4;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
  };
  Tet(const int v0, const int v1, const int v2, const int v3) : Element() {
    type = TET;
    vert_num = 4;
    vert_ids = new int[vert_num];
    face_num = 4;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
    Set_vert_ids(v0,v1,v2,v3);
  };
  ~Tet(){};//fill in destructor
  void Orient();
  void Set_vert_ids(const int v0, const int v1, const int v2, const int v3);
  void Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3);
  int Find_match_face(const int v0,const int v1,const int v2, const int v3);
};//class Tet

class Hex: public Element{
private:
public:
  Hex() : Element(){
    type = HEX;
    vert_num = 8;
    vert_ids = new int[vert_num];
    face_num = 6;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
  };
  Hex(const int v0, const int v1, const int v2, const int v3,
      const int v4, const int v5, const int v6, const int v7) : Element(){
    type = HEX;
    vert_num = 8;
    vert_ids = new int[vert_num];
    face_num = 6;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
    Set_vert_ids(v0,v1,v2,v3,v4,v5,v6,v7);
  };
  ~Hex(){};//fill in destructor
  void Orient();
  void Set_vert_ids(const int v0, const int v1, const int v2, const int v3,
		    const int v4, const int v5, const int v6, const int v7);
  void Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3);
  int Find_match_face(const int v0,const int v1,const int v2, const int v3);
};//class Hex

class Pri: public Element{
private:
public:
  Pri() : Element(){
    type = PRI;
    vert_num = 6;
    vert_ids = new int[vert_num];
    face_num = 5;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
  };
  Pri(const int v0, const int v1, const int v2, const int v3,
      const int v4, const int v5) : Element(){
    type = PRI;
    vert_num = 6;
    vert_ids = new int[vert_num];
    face_num = 5;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
    Set_vert_ids(v0,v3,v5,v1,v2,v4);
  };
  ~Pri(){};//fill in destructor
  void Orient();
  void Set_vert_ids(const int v0, const int v1, const int v2, const int v3, 
		    const int v4, const int v5);
  void Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3);
  int Find_match_face(const int v0,const int v1,const int v2, const int v3);
};//class Pri

class Pyr: public Element{
private:
public:
  Pyr() : Element(){
    type = PYR;
    vert_num = 5;
    vert_ids = new int[vert_num];
    face_num = 5;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
  };
  Pyr(const int v0, const int v1, const int v2, const int v3,
      const int v4) : Element(){
    type = PYR;
    vert_num = 5;
    vert_ids = new int[vert_num];
    face_num = 5;
    conn_element = new int[face_num];
    conn_face = new int[face_num];
    Set_vert_ids(v0,v1,v2,v3,v4);
  };
  ~Pyr(){};//fill in destructor
  void Orient();
  void Set_vert_ids(const int v0, const int v1, const int v2, const int v3, 
		    const int v4);
  void Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3);
  int Find_match_face(const int v0,const int v1,const int v2, const int v3);
};//class Pyr

class BFace{
private:
public:
  int vert_num;
  int* vert_ids;
  BFace(const int vert_num_, const int v0, const int v1, const int v2, const int v3){
    if((vert_num_ > 4)||(vert_num_ < 3)){
      fprintf(stderr,"Invalid # of vertices in BFace!\n");
      exit(-1);
    }
    else{
      vert_num = vert_num_;
      vert_ids = new int[vert_num];
      vert_ids[0] = v0;
      vert_ids[1] = v1;
      vert_ids[2] = v2;
      if(vert_num == 4) vert_ids[3] = v3;
    };
  };
  ~BFace(){ 
    vert_num=0;
    delete [] vert_ids; 
  };
  void Get_face_ids(int &v0, int &v1, int &v2, int &v3){
    v0 = vert_ids[0];
    v1 = vert_ids[1];
    v2 = vert_ids[2];
    if(vert_num == 4){
      v3 = vert_ids[3];
    }
    else{
      v3 = NONE;
    }
  }
};

typedef vector<BFace*> BFace_List;

class Boundary{
private:
  char type;
  char *name;
  int lines_num;
  int curved;
  char curve_type;
  char** lines;
  char* inline_params;
public:
  BFace_List bface_list;
  Boundary(){
    int j;
    name = (char*) calloc(30,sizeof(char));
    lines = (char**) calloc(MAXBC,sizeof(char*));
    for(j=0;j<MAXBC;++j) lines[j] = (char*) calloc(MAXBCLINE,sizeof(char));
    curved = -1;
  };
  ~Boundary(){};//fill in destructor
  inline void Set_name(const char* name_){ name = strdup(name_); };
  inline char* Get_name(){ return name; };
  inline void Set_curved(){ curved = 1; };
  inline int Get_curved(){ return curved; };
  inline void Set_curve_type(const char curve_type_){ curve_type = curve_type_; };
  inline char Get_curve_type(){ return curve_type; };
  inline void Set_type(const char type_){ type = type_; };
  inline char Get_type(){ return type; };
  inline void Set_number_of_lines(const int lines_num_){ 
    if(lines_num_ <= MAXBC){
      lines_num = lines_num_; 
    }
    else{
      fprintf(stderr,"\nWARNING: The number of BC lines is greater than MAXBC !!\n");
      fprintf(stderr,"WARNING: Only %ld lines will be used for %s\n",MAXBC,strdup(Get_name()));
      lines_num = MAXBC;
    };
    if(lines_num<0) lines_num = 0;
  };
  inline int Get_number_of_lines(){ return lines_num; };
  inline void Set_lines(int linenum_, const char* line_){ lines[linenum_] = strdup(line_); };
  inline char* Get_lines(const int linenum_){ return lines[linenum_]; };
  inline void Set_inline_params(const char* line_){ inline_params = strdup(line_); };
  inline char* Get_inline_params(){ return inline_params; };
  inline void Add_BFace(const int vert_num_, const int v0, const int v1, const int v2, const int v3){
    BFace* bface = new BFace(vert_num_,v0,v1,v2,v3);
    bface_list.push_back(bface);
  };
  inline int Get_number_of_bfaces(){ return bface_list.size(); };
};


typedef vector<Vertex*> Vertex_List;
typedef vector<Element*> Element_List;
typedef vector<Boundary*> Boundary_List;


class Periodic_Info{
private:
public:
  double xpmin, xpmax, ypmin, ypmax, zpmin, zpmax;
  Periodic_Info(){};
  ~Periodic_Info(){};
};



class Mesh{
private:
  Vertex_List vertex_list;
  Element_List element_list;
  Boundary_List boundary_list;
  char title[BUFSIZ];//project title 
  int tet_num, hex_num, pri_num, pyr_num;
public:
  int **tetvnum, **hexvnum, **privnum, **pyrvnum;
  Periodic_Info* periodic_info;
  Mesh(){
    int i;
    tet_num = hex_num = pri_num = pyr_num = 0;
    periodic_info = new Periodic_Info();
    
    tetvnum = new int*[4]; for(i=0;i<4;++i) tetvnum[i] = new int[4];
    tetvnum[0][0]=0; tetvnum[0][1]=1; tetvnum[0][2]=2; tetvnum[0][3]=-1;
    tetvnum[1][0]=0; tetvnum[1][1]=1; tetvnum[1][2]=3; tetvnum[0][3]=-1;
    tetvnum[2][0]=1; tetvnum[2][1]=2; tetvnum[2][2]=3; tetvnum[0][3]=-1;
    tetvnum[3][0]=0; tetvnum[3][1]=2; tetvnum[3][2]=3; tetvnum[0][3]=-1;
    
    privnum = new int*[5]; for(i=0;i<5;++i) privnum[i] = new int[4];
    privnum[0][0]=0; privnum[0][1]=1; privnum[0][2]=2; privnum[0][3]=3;
    privnum[1][0]=0; privnum[1][1]=1; privnum[1][2]=4; privnum[1][3]=-1;
    privnum[2][0]=1; privnum[2][1]=2; privnum[2][2]=5; privnum[2][3]=4;
    privnum[3][0]=3; privnum[3][1]=2; privnum[3][2]=5; privnum[3][3]=-1;
    privnum[4][0]=0; privnum[4][1]=3; privnum[4][2]=5; privnum[4][3]=4;
    
    pyrvnum = new int*[5]; for(i=0;i<5;++i) pyrvnum[i] = new int[4];
    pyrvnum[0][0]=0; pyrvnum[0][1]=1; pyrvnum[0][2]=2; pyrvnum[0][3]=3;
    pyrvnum[1][0]=0; pyrvnum[1][1]=1; pyrvnum[1][2]=4; pyrvnum[1][3]=-1;
    pyrvnum[2][0]=1; pyrvnum[2][1]=2; pyrvnum[2][2]=4; pyrvnum[2][3]=-1;
    pyrvnum[3][0]=3; pyrvnum[3][1]=2; pyrvnum[3][2]=4; pyrvnum[3][3]=-1;
    pyrvnum[4][0]=0; pyrvnum[4][1]=3; pyrvnum[4][2]=4; pyrvnum[4][3]=-1;
    
    hexvnum = new int*[6]; for(i=0;i<6;++i) hexvnum[i] = new int[4];
    hexvnum[0][0]=0; hexvnum[0][1]=1; hexvnum[0][2]=2; hexvnum[0][3]=3;
    hexvnum[1][0]=0; hexvnum[1][1]=1; hexvnum[1][2]=5; hexvnum[1][3]=4;
    hexvnum[2][0]=1; hexvnum[2][1]=2; hexvnum[2][2]=6; hexvnum[2][3]=5;
    hexvnum[3][0]=3; hexvnum[3][1]=2; hexvnum[3][2]=6; hexvnum[3][3]=7;
    hexvnum[4][0]=0; hexvnum[4][1]=3; hexvnum[4][2]=7; hexvnum[4][3]=4;
    hexvnum[5][0]=4; hexvnum[5][1]=5; hexvnum[5][2]=6; hexvnum[5][3]=7;
  };
  ~Mesh(){};//fill in destructor
  void Read_input_files();
  void Orient_elements();
  void Compute_connectivity(){};
  void Dump_rea_files(){};
  
  inline void Add_vertex(const double x_, const double y_, const double z_);
  
  int Get_number_of_vertices();
  
  void Add_tet(const int v0, const int v1, const int v2, const int v3){
    Tet* tet = new Tet(v0,v1,v2,v3);
    tet->Get_mesh(this);
    element_list.push_back(tet);
    ++tet_num;
  };
  
  void Add_hex(const int v0, const int v1, const int v2, const int v3,
		      const int v4, const int v5, const int v6, const int v7){
    Hex* hex = new Hex(v0,v1,v2,v3,v4,v5,v6,v7);
    hex->Get_mesh(this);
    element_list.push_back(hex);
    ++hex_num;
  };
  
  void Add_pri(const int v0, const int v1, const int v2, const int v3,
		      const int v4, const int v5){
    Pri* pri = new Pri(v0,v1,v2,v3,v4,v5);
    pri->Get_mesh(this);
    element_list.push_back(pri);
    ++pri_num;
  };

  void Add_pyr(const int v0, const int v1, const int v2, const int v3,
		      const int v4){
    Pyr* pyr = new Pyr(v0,v1,v2,v3,v4);
    pyr->Get_mesh(this);
    element_list.push_back(pyr);
    ++pyr_num;
  };

  inline int Get_number_of_elements(){ return element_list.size(); };
  inline int Get_number_of_tets(){ return tet_num; };
  inline int Get_number_of_hexs(){ return hex_num; };
  inline int Get_number_of_pris(){ return pri_num; };
  inline int Get_number_of_pyrs(){ return pyr_num; };
  inline int Get_number_of_bcs(){ return boundary_list.size(); };
  
  inline void Add_boundary(const char* name_){
    Boundary* bc = new Boundary();
    bc->Set_name(name_);
    boundary_list.push_back(bc);
  };
  
  void Find_match_element_face(int &eid_, int &fid_, const int v0, const int v1, const int v2, const int v3);
  int Find_match_vertex(const double x, const double y, const double z);
  double Determinant(int v0, int v1, int v2, int v3);
  void Connect();
  void Set_Boundary_Conditions();
  void Dump_rea();
  void Dump_sum_rea();
  void Dump_prod_rea();
  void Dump_plt();
  void Dump_walls();
  int curvedfaces_check();
  void curvedfaces_fix();
};



static char *findSection (char *name, char *buf, FILE *fp){
  char *p;
  while (p = fgets (buf, BUFSIZ, fp))
    if (strstr (p, name))
      break;
  return p;
};



#endif

