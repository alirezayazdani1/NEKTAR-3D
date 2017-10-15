/****************************************************************/
/*    GRIDGEN  -->> NEKTAR3D                                    */
/*    Author: Igor Pivkin, Brown University, RI                 */
/*    updated by Leopold Grinberg, Brown University, RI         */
/*    email: piv@@cfm.brown.edu                                  */
/*    08/04/02                                                  */
/*    g++ -O3 gg2rea3d.C -o gg2rea3d                            */
/****************************************************************/


#include "gg_curve.H"
#define V15

void Tet::Orient(){
  int j,v0,v1,v2,v3;
  v1 = max(vert_ids[0], vert_ids[1]);
  v2 = max(vert_ids[2], vert_ids[3]);
  v0 = max(v1,v2);
  
  v1 = -1000;
  // find next max vert id
  for(j=0;j<4;++j)
    if(vert_ids[j] < v0 && vert_ids[j] > v1)
      v1 = vert_ids[j];
  
  v2 = -1000;
  for(j=0;j<4;++j)
    if(vert_ids[j] < v1 && vert_ids[j] > v2)
      v2 = vert_ids[j];
  
  v3 = -1000;
  for(j=0;j<4;++j)
    if(vert_ids[j] < v2 && vert_ids[j] > v3)
      v3 = vert_ids[j];
  
  if(v0 == -1000) fprintf(stderr, " v0 f.d\n");
  if(v1 == -1000) fprintf(stderr, " v1 f.d\n");
  if(v2 == -1000) fprintf(stderr, " v2 f.d\n");
  if(v3 == -1000){
    fprintf(stderr, " v3 f.d: %d %d %d %d \n",v0,v1,v2,v3);
    fprintf(stderr, " tet_ids: %d %d %d %d \n",
	    vert_ids[0], vert_ids[1],vert_ids[2], vert_ids[3]);
  };
  
  if(mesh->Determinant(v3,v2,v1,v0) < 0){
    vert_ids[0] = v2; vert_ids[1] = v3; vert_ids[2] = v1; vert_ids[3] = v0;
  }
  else{
    vert_ids[0] = v3; vert_ids[1] = v2; vert_ids[2] = v1; vert_ids[3] = v0;
  };
};

void Hex::Orient(){
  int v0,v1,v2,v3,v4,v5,v6,v7;
  v0 = vert_ids[0]; v1 = vert_ids[1]; v2 = vert_ids[2]; v3 = vert_ids[3];
  v4 = vert_ids[4]; v5 = vert_ids[5]; v6 = vert_ids[6]; v7 = vert_ids[7];
  if(mesh->Determinant(v0,v1,v3,v4) < 0.0){
    vert_ids[0] = v4; vert_ids[1] = v5; vert_ids[2] = v6; vert_ids[3] = v7; 
    vert_ids[4] = v0; vert_ids[5] = v1; vert_ids[6] = v2; vert_ids[7] = v3; 
  }
};

void Pyr::Orient(){
  int v0,v1,v2,v3,v4;
  v0 = vert_ids[4]; v1 = vert_ids[3]; v2 = vert_ids[0]; v3 = vert_ids[1];
  v4 = vert_ids[2];
  if(mesh->Determinant(v3,v2,v1,v0) > 0.0){
    vert_ids[0] = v3; vert_ids[1] = v2; vert_ids[2] = v1; vert_ids[3] = v4; 
    vert_ids[4] = v0;
  }
  else{
    vert_ids[0] = v1; vert_ids[1] = v2; vert_ids[2] = v3; vert_ids[3] = v4; 
    vert_ids[4] = v0;
  };
};

void Pri::Orient(){
  int v0,v1,v2,v3,v4,v5;
  v0 = max(vert_ids[0], vert_ids[1]); v1 = max(v0, vert_ids[2]);
  v2 = max(vert_ids[3], vert_ids[4]); v3 = max(v2, vert_ids[5]);
  if(v1>v3){
    v0 = vert_ids[0]; v1 = vert_ids[1]; v2 = vert_ids[2];
    vert_ids[1]=vert_ids[3]; vert_ids[0]=vert_ids[4]; vert_ids[2]=vert_ids[5];
    vert_ids[4] = v0; vert_ids[3] = v1; vert_ids[5] = v2;
  }
  v0 = vert_ids[0]; v1 = vert_ids[1]; v2 = vert_ids[2]; 
  v3 = vert_ids[3]; v4 = vert_ids[4]; v5 = vert_ids[5];
  if(v3 > v4 && v3 > v5){
    vert_ids[0] = v2; vert_ids[1] = v1; vert_ids[2] = v4; 
    vert_ids[3] = v5; vert_ids[4] = v0; vert_ids[5] = v3;
  }
  else if(v4 > v3 && v4 > v5){
    vert_ids[0] = v0; vert_ids[1] = v2; vert_ids[2] = v5; 
    vert_ids[3] = v3; vert_ids[4] = v1; vert_ids[5] = v4;
  }
  else if(v5 > v3 && v5 > v4){
    vert_ids[0] = v1; vert_ids[1] = v0; vert_ids[2] = v3; 
    vert_ids[3] = v4; vert_ids[4] = v2; vert_ids[5] = v5;
  }
};


void Element::Get_mesh(Mesh* mesh_){ mesh = mesh_;};

int Element::Get_number_of_vertices(){ return vert_num; };

int Element::Get_vert_id(const int vid_){
  if((vid_<0)||(vid_>vert_num)){ 
    //fprintf(stderr,"Invalid vertex id in Get_vert_id!");
    return NONE;
    //exit(-1);
  }
  return vert_ids[vid_];
};
int Element::Get_number_of_faces(){ return face_num; };

void Element::Set_myid(const int myid_){ myid = myid_; };

void Element::Set_vert_id(const int id_, const int vid_){
  if((id_>vert_num)||(id_<0)){
    fprintf(stderr,"Invalid vertex id in Set_vert_id!");
  }
  else{
    vert_ids[id_]=vid_;
  }
};

void Tet::Set_vert_ids(const int v0, const int v1, const int v2, const int v3){
  vert_ids[0]=v0; vert_ids[1]=v1; vert_ids[2]=v2; vert_ids[3]=v3;
};

void Pyr::Set_vert_ids(const int v0, const int v1, const int v2, const int v3, const int v4){
  vert_ids[0]=v0; vert_ids[1]=v1; vert_ids[2]=v2; vert_ids[3]=v3; 
  vert_ids[4]=v4;
};  

void Pri::Set_vert_ids(const int v0, const int v1, const int v2, const int v3, const int v4, const int v5){
  vert_ids[0]=v0; vert_ids[1]=v1; vert_ids[2]=v2; vert_ids[3]=v3;
  vert_ids[4]=v4; vert_ids[5]=v5;
};

void Hex::Set_vert_ids(const int v0, const int v1, const int v2, const int v3, const int v4, const int v5, const int v6, const int v7){
  vert_ids[0]=v0; vert_ids[1]=v1; vert_ids[2]=v2; vert_ids[3]=v3;
  vert_ids[4]=v4; vert_ids[5]=v5; vert_ids[6]=v6; vert_ids[7]=v7;
};

void Hex::Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3){
  if(fid_>=face_num){
    fprintf(stderr,"Invalid fid in Hex::Get_face_ids!!\n");
  }
  else{
    v0 = Get_vert_id(mesh->hexvnum[fid_][0]); 
    v1 = Get_vert_id(mesh->hexvnum[fid_][1]);
    v2 = Get_vert_id(mesh->hexvnum[fid_][2]); 
    v3 = Get_vert_id(mesh->hexvnum[fid_][3]);
  }
};

void Tet::Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3){
  if(fid_>=face_num){
    fprintf(stderr,"Invalid fid in Hex::Get_face_ids!!\n");
  }
  else{
    v0 = Get_vert_id(mesh->tetvnum[fid_][0]); 
    v1 = Get_vert_id(mesh->tetvnum[fid_][1]);
    v2 = Get_vert_id(mesh->tetvnum[fid_][2]); 
    v3 = NONE;
  }
};

void Pri::Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3){
  if(fid_>=face_num){
    fprintf(stderr,"Invalid fid in Pri::Get_face_ids!!\n");
  }
  else{
    v0 = Get_vert_id(mesh->privnum[fid_][0]); 
    v1 = Get_vert_id(mesh->privnum[fid_][1]);
    v2 = Get_vert_id(mesh->privnum[fid_][2]); 
    if(mesh->privnum[fid_][3] < 0) {
      v3 = NONE;
    }
    else{
      v3 = Get_vert_id(mesh->privnum[fid_][3]);
    }
  }
};

void Pyr::Get_face_ids(const int fid_, int &v0, int &v1, int &v2, int &v3){
  if(fid_>=face_num){
    fprintf(stderr,"Invalid fid in Pyr::Get_face_ids!!\n");
  }
  else{
    v0 = Get_vert_id(mesh->pyrvnum[fid_][0]); 
    v1 = Get_vert_id(mesh->pyrvnum[fid_][1]);
    v2 = Get_vert_id(mesh->pyrvnum[fid_][2]); 
    if(mesh->pyrvnum[fid_][3] < 0) {
      v3=NONE;
    }
    else{
      v3 = Get_vert_id(mesh->pyrvnum[fid_][3]);
    }
  }
};

int Tet::Find_match_face(const int v0,const int v1,const int v2, const int v3){
  //return the face id which match the vertices
  if(v3 >= 0){
    return NONE;
  }
  else{
    int fid,vid,cnt;
    for(fid=0;fid<face_num;++fid){
      cnt=0;
      for(vid=0;vid<3;++vid){
	if(v0 == Get_vert_id(mesh->tetvnum[fid][vid])) ++cnt;
	if(v1 == Get_vert_id(mesh->tetvnum[fid][vid])) ++cnt;
	if(v2 == Get_vert_id(mesh->tetvnum[fid][vid])) ++cnt;
      }
      if(cnt == 3) return fid;
    }
    return NONE;
  }
};

int Hex::Find_match_face(const int v0,const int v1,const int v2, const int v3){
  //return the face id which match the vertices
  if(v3 < 0){
    return NONE;
  }
  else{
    int fid,vid,cnt;
    for(fid=0;fid<face_num;++fid){
      cnt=0;
      for(vid=0;vid<4;++vid){
	if(v0 == Get_vert_id(mesh->hexvnum[fid][vid])) ++cnt;
	if(v1 == Get_vert_id(mesh->hexvnum[fid][vid])) ++cnt;
	if(v2 == Get_vert_id(mesh->hexvnum[fid][vid])) ++cnt;
	if(v3 == Get_vert_id(mesh->hexvnum[fid][vid])) ++cnt;
      }
      if(cnt == 4) return fid;
    }
    return NONE;
  }
};

int Pri::Find_match_face(const int v0,const int v1,const int v2, const int v3){
  //return the face id which match the vertices
  int fid,vid,cnt;
  
/*  fprintf(stderr,"input vertices: %ld %ld %ld %ld\n",v0,v1,v2,v3);
  fprintf(stderr,"my faces:\n");*/
  for(fid=0;fid<face_num;++fid){
    for(vid=0;vid<4;++vid){
      fprintf(stderr,"%ld  ",Get_vert_id(mesh->privnum[fid][vid]));
    }
    fprintf(stderr,"\n");
  }
  
  if(v3 < 0){
    for(fid=1;fid<face_num;fid=fid+2){
      cnt=0;
      for(vid=0;vid<3;++vid){
	if(v0 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
	if(v1 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
	if(v2 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
      }
      if(cnt == 3) return fid;
    }
    return NONE;
  }
  else{
    for(fid=0;fid<face_num;fid=fid+2){
      cnt=0;
      for(vid=0;vid<4;++vid){
	if(v0 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
	if(v1 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
	if(v2 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
	if(v3 == Get_vert_id(mesh->privnum[fid][vid])) ++cnt;
      }
      if(cnt == 4) return fid;
    }
    return NONE;
  }
};

int Pyr::Find_match_face(const int v0,const int v1,const int v2, const int v3){
  //return the face id which match the vertices
  int fid,vid,cnt;
  if(v3 < 0){
    for(fid=1;fid<face_num;++fid){
      cnt=0;
      for(vid=0;vid<3;++vid){
	if(v0 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
	if(v1 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
	if(v2 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
      }
      if(cnt == 3) return fid;
    }
    return NONE;
  }
  else{
    fid=0;
    cnt=0;
    for(vid=0;vid<4;++vid){
      if(v0 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
      if(v1 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
      if(v2 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
      if(v3 == Get_vert_id(mesh->pyrvnum[fid][vid])) ++cnt;
    }
    if(cnt == 4) return fid;
    return NONE;
  }
};



double Mesh::Determinant(int v0, int v1, int v2, int v3){
  //returns the determinant for 4 given vertices
  double ax = vertex_list[v1]->x - vertex_list[v0]->x;
  double ay = vertex_list[v1]->y - vertex_list[v0]->y;
  double az = vertex_list[v1]->z - vertex_list[v0]->z;
  
  double bx = vertex_list[v2]->x - vertex_list[v0]->x;
  double by = vertex_list[v2]->y - vertex_list[v0]->y;
  double bz = vertex_list[v2]->z - vertex_list[v0]->z;
  
  double cx = vertex_list[v3]->x - vertex_list[v0]->x;
  double cy = vertex_list[v3]->y - vertex_list[v0]->y;
  double cz = vertex_list[v3]->z - vertex_list[v0]->z;    
  
  double dx = ay*bz - az*by;
  double dy = az*bx - ax*bz;
  double dz = ax*by - ay*bx;
  
  return cx*dx + cy*dy + cz*dz;
};

void Mesh::Orient_elements(){
  for(int i=0;i<Get_number_of_elements();++i)
    element_list[i]->Orient();
};

void Mesh::Add_vertex(const double x_, const double y_, const double z_){
    Vertex* v = new Vertex(x_,y_,z_);
    vertex_list.push_back(v);
  };




int Mesh::Get_number_of_vertices(){ return vertex_list.size(); };

void Mesh::Read_input_files(){
  char buf[BUFSIZ];
  char bufa[BUFSIZ];
  char string[100];

  FILE *griin;
  FILE *infoin;
  double x_,y_,z_;
  int v0,v1,v2,v3,v4,v5,v6,v7;
  int i,j;
  int Npt,Ntmp1,Ntmp2,Nbc,inl;

  double s1[3], s2[3];
  double s1cs2[3], slength, weigth, mulfactor ;
  double xscale = 1.0;
  double yscale = 1.0;
  double zscale = 1.0;

  char   btype;
  int    curved;
  int    bcalcnormal = 0;
  char   *inlineparam;

  fprintf(stderr, "++++++++++++++++++++\n"); 
  fprintf(stderr,"GRIDGEN 3D -> NEKTAR\n");
  fprintf(stderr, "Enter project name: "); 
  fgets(title, BUFSIZ, stdin);
  sprintf(bufa, "%s.inp", strtok(title, "\n"));
  if((griin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    exit(-1);
  }
  else{
    fprintf(stderr, "Reading %s file\n",bufa);
  }
  
  fprintf(stderr, "Enter Scale in X,Y, and Z directions (Default 1. 1. 1.): ");
  fgets(buf, BUFSIZ, stdin);
  if(sscanf(buf, "%lf %lf %lf", &xscale, &yscale, &zscale)!=3)
  {
	fprintf(stderr, "Scales in all direction will be 1.0\n");
  }
  else
  {
	fprintf(stderr, "Scales will be %lf %lf %lf\n", xscale, yscale, zscale);

  }


  rewind(griin);
  findSection ("Nodes",buf,griin);
#ifndef V15
  fgets(buf,  BUFSIZ, griin);
  sscanf(buf, "%d", &Npt); 
#else
  sscanf(buf, "%s %d", string, &Npt);  
#endif


  // read nodes
  fprintf(stderr, "Reading node coords ---> ");
  
  for(i=0;i<Npt;++i){
    fscanf(griin, "%lf%lf%lf", &x_,&y_,&z_);
    x_ = xscale*x_;
    y_ = yscale*y_;
    z_ = zscale*z_;
    Add_vertex(x_,y_,z_);
  }

  fprintf(stderr,"Number of nodes: %ld\n",Get_number_of_vertices());

  fprintf(stderr, "Reading elements ---> ");
  rewind(griin);
  findSection ("Elements",buf,griin);
  while(!strstr(buf,"Variables")){
    fgets(buf, BUFSIZ, griin);
    if(!strstr(buf,"Variables")){
      sscanf(buf, "%d%d%d%d%d%d%d%d%d%d",&Ntmp1,&Ntmp2,&v0,&v1,&v2,&v3,&v4,&v5,&v6,&v7);
      switch(Ntmp1){
      case 1:
	Add_tet(v0-1,v1-1,v2-1,v3-1);
	break;
      case 2:
	Add_hex(v0-1,v1-1,v3-1,v2-1,v4-1,v5-1,v7-1,v6-1);
	break;
      case 3:
	Add_pri(v0-1,v1-1,v2-1,v3-1,v4-1,v5-1);
	break;
      case 4:
	Add_pyr(v0-1,v1-1,v2-1,v3-1,v4-1);
	break;
      default:
	fprintf(stderr, "Unknown type of element!!\n");
      }
    }
  }
  for(i=0;i<Get_number_of_elements();++i) 
    element_list[i]->Set_myid(i);

  fprintf(stderr, "Number of elements: %d\n",Get_number_of_elements());
  
  //Boundary Conditions
  fprintf(stderr, "Reading boundary conditions\n");
  rewind(griin);
  findSection ("Boundary Table",buf,griin);

#ifndef V15
  fgets(buf,  BUFSIZ, griin);

  sscanf(buf, "%d", &Nbc);
  for(i=0;i<Nbc;++i){
    fgets(buf,BUFSIZ, griin);
    Add_boundary(strdup(buf));
  }
#else
  sscanf(buf, "%s %s %d", string, string, &Nbc);
  for(i=0;i<Nbc;++i){
    fgets(buf,BUFSIZ, griin);
    Add_boundary(strdup(buf));
  }
#endif



  //Create .info file if not exist
  int new_info = 0;
  sprintf(bufa, "%s.info", strtok(title, "\n"));
  if((infoin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    fprintf(stderr, "File %s will be generated automatically now.\n",bufa);
    new_info=1; 
  }
  else{
    fprintf(stderr, "Do you want to use current version of %s file (y or n)?: ",bufa);
    fgets(buf,BUFSIZ, stdin);
    if(buf[0]=='n') new_info=1;
  }
  if(new_info == 1){//generating new .info file
    infoin = fopen(bufa, "w");
    fprintf(infoin, "***BOUNDARY CONDITIONS***\n");
    for(i=0;i<Get_number_of_bcs();++i){
      fprintf(infoin, "+++ do NOT modify the NEXT line +++\n");
      fprintf(infoin,strdup(boundary_list[i]->Get_name()));
      fprintf(infoin, "0 CURVED\n");
      fprintf(infoin, "a CURVE_TYPE\n");
      fprintf(infoin, "W TYPE\n");
      fprintf(infoin, "0 LINES\n");
      fprintf(infoin, "\n");
    }

    fprintf(infoin, "\n+++ do NOT modify the NEXT line +++\n");
    fprintf(infoin, "*** CURVED SIDE DATA ***\n");
    fprintf(infoin, "0 Number of curve types\n");
    fprintf(infoin, "Sphere\n");
    fprintf(infoin, "0.0 0.0 0.0 1.0 A\n");
    

    fprintf(infoin, "\n+++ do NOT modify the NEXT line +++\n");
    fprintf(infoin, "*** PERIODIC DATA ***\n");
    fprintf(infoin, "0.0 XPERMIN\n");
    fprintf(infoin, "0.0 XPERMAX\n");
    fprintf(infoin, "0.0 YPERMIN\n");
    fprintf(infoin, "0.0 YPERMAX\n");
    fprintf(infoin, "0.0 ZPERMIN\n");
    fprintf(infoin, "0.0 ZPERMAX\n");

    fprintf(infoin, "\n+++ do NOT modify the NEXT line +++\n");
    fprintf(infoin, "*** HEAD ***\n");
    fprintf(infoin, "22 LINES\n");
    fprintf(infoin, "****** PARAMETERS *****\n");
    fprintf(infoin, "GRIDGEN 3D -> NEKTAR\n");
    fprintf(infoin, "3 DIMENSIONAL RUN\n");
    fprintf(infoin, "15 PARAMETERS FOLLOW\n");
    fprintf(infoin, "0.01     KINVIS\n"); 
    fprintf(infoin, "4        MODES\n"); 
    fprintf(infoin, "0.001    DT\n"); 
    fprintf(infoin, "10000    NSTEPS\n");
    fprintf(infoin, "1        EQTYPE\n");
    fprintf(infoin, "2        INTYPE\n");
    fprintf(infoin, "1000     IOSTEP\n");
    fprintf(infoin, "10       HISSTEP\n");
    fprintf(infoin, "3        PRECON\n");
    fprintf(infoin, "1        NINLETS\n");
    fprintf(infoin, "1        NOUTLETS\n");
    fprintf(infoin, "1.0      DPSCAL\n");
    fprintf(infoin, "1.0      DVSCAL\n");
    fprintf(infoin, "4        NPODORDER\n");
    fprintf(infoin, "2        NPPODORDER\n");
    fprintf(infoin, "0  Lines of passive scalar data follows\n");
    fprintf(infoin, "0  LOGICAL SWITCHES FOLLOW\n");
    fprintf(infoin, "Dummy line from old nekton file\n"); 

    fprintf(infoin, "\n+++ do NOT modify the NEXT line +++\n");
    fprintf(infoin, "*** TAIL ***\n");
    fprintf(infoin, "16 LINES\n");

    fprintf(infoin, "***** NO THERMAL BOUNDARY CONDITIONS *****\n");
    fprintf(infoin, "4 INITIAL CONDITIONS *****\n");
    fprintf(infoin, "Given\n");
    fprintf(infoin, "u = 0.0\n");
    fprintf(infoin, "v = 0.0\n");
    fprintf(infoin, "w = 0.0\n");
    fprintf(infoin, "***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q\n");
    fprintf(infoin, "0 Lines of Drive force data follow\n");
    fprintf(infoin, "***** Variable Property Data ***** Overrrides Parameter data.\n");
    fprintf(infoin, "1 Lines follow.\n");
    fprintf(infoin, "0 PACKETS OF DATA FOLLOW\n");
    fprintf(infoin, "***** HISTORY AND INTEGRAL DATA *****\n");
    fprintf(infoin, "0 POINTS.  Hcode, I,J,H,IEL\n");
    fprintf(infoin, "UVWP  H  1 0 0 10\n");
    fprintf(infoin, "***** OUTPUT FIELD SPECIFICATION *****\n");
    fprintf(infoin, "0 SPECIFICATIONS FOLLOW\n");
    fclose(infoin); 
    
    fprintf(stderr, "Please modify this file and press ENTER  ");
    fgets(buf,BUFSIZ, stdin);
  }
  
  //Read Boundary Conditions
  sprintf(bufa, "%s.info", strtok(title, "\n"));
  if((infoin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    exit(-1);
  }
  else{
    fprintf(stderr, "Reading %s file\n",bufa);
  }
  
  for (i=0;i<Get_number_of_bcs();++i){ 
    rewind(infoin);
    findSection (strdup(boundary_list[i]->Get_name()),buf,infoin);
    if (!strstr(buf,strdup(boundary_list[i]->Get_name())))
      {
	fprintf(stderr,"Reading .info file: do not see any boundary conditions for %s", strdup(boundary_list[i]->Get_name()) );
	exit(-1);
      }
    fgets(buf,BUFSIZ, infoin);
    sscanf(buf, "%d", &inl);
    if(inl==1) boundary_list[i]->Set_curved();
    fgets(buf,BUFSIZ, infoin);
    boundary_list[i]->Set_curve_type(buf[0]);
    fgets(buf,BUFSIZ, infoin);
    boundary_list[i]->Set_type(buf[0]);
    fgets(buf,  BUFSIZ, infoin);
    sscanf(buf, "%d", &inl);

    if (inl>=0){
      boundary_list[i]->Set_number_of_lines(inl);
      for (j=0;j<boundary_list[i]->Get_number_of_lines();++j){
        fgets(buf,BUFSIZ, infoin);
       boundary_list[i]->Set_lines(j,buf);
     }
    
      //Read INLINE PARAMETERS if they are there
      fgets(buf,BUFSIZ, infoin);
      if(strstr(buf,"INLINE")){
        fgets(buf,BUFSIZ, infoin);
        boundary_list[i]->Set_inline_params(buf);
      }
      else{
        sprintf(buf,"0. 0. 0.\n");
        boundary_list[i]->Set_inline_params(buf);
      }
    }
    else{//modified by Leopold Grinberg;
         //in order to set outflow boundary cond for a press. 
      boundary_list[i]->Set_number_of_lines(0);
      fgets(buf,BUFSIZ, infoin);
      boundary_list[i]->Set_inline_params(buf);
    }
    //End of reading INLINE PARAMETERS
  }
  
  fprintf(stderr,"Reading periodic data\n");
  rewind(infoin);
  findSection ("*** PERIODIC DATA ***",buf,infoin);
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->xpmin));
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->xpmax));
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->ypmin));
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->ypmax));
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->zpmin));
  fgets(buf,  BUFSIZ, infoin);
  sscanf(buf, "%lf", &(periodic_info->zpmax));

  fclose(infoin);

  fprintf(stderr,"Reading boundary faces\n");
  int Nbf,BC_type,vert_num_;
  rewind(griin);
  findSection ("Boundary Faces",buf,griin);
#ifndef V15
  fgets(buf,  BUFSIZ, griin);
  sscanf(buf, "%d", &Nbf);
#else
  sscanf(buf, "%s %s %d", string, string, &Nbf);
#endif


  double mid_edge_0_x,mid_edge_0_y,mid_edge_0_z,
         mid_edge_1_x,mid_edge_1_y,mid_edge_1_z,
         mid_edge_2_x,mid_edge_2_y,mid_edge_2_z;
  double norm_v0,norm_v1,norm_v2;

  for(i=0;i<Nbf;++i){
    fgets(buf, BUFSIZ, griin);
    sscanf(buf, "%d%d%d%d%d%d",&BC_type,&vert_num_,&v0,&v1,&v2,&v3);

//    calculate normal  here using the vertex

     btype =  boundary_list[BC_type-1]->Get_type();
     curved = boundary_list[BC_type-1]->Get_curved();

     bcalcnormal = 0;

     if ((curved) && (btype =='W'))
     {
        inlineparam = boundary_list[BC_type-1]->Get_inline_params();
        if(strstr(inlineparam,"SI"))
        {        
        bcalcnormal = 1;
        }

     }

     if ((vert_num_ ==3) && (bcalcnormal))
     {
         s1[0] = vertex_list[v2-1]->x - vertex_list[v0-1]->x;
         s1[1] = vertex_list[v2-1]->y - vertex_list[v0-1]->y;
         s1[2] = vertex_list[v2-1]->z - vertex_list[v0-1]->z;
         
         s2[0] = vertex_list[v1-1]->x - vertex_list[v0-1]->x;
         s2[1] = vertex_list[v1-1]->y - vertex_list[v0-1]->y;
         s2[2] = vertex_list[v1-1]->z - vertex_list[v0-1]->z;

	 s1cs2[0]  = s1[1]*s2[2]-s1[2]*s2[1];
	 s1cs2[1]  = s1[2]*s2[0]-s1[0]*s2[2];
	 s1cs2[2]  = s1[0]*s2[1]-s1[1]*s2[0];
 
         slength = sqrt(s1cs2[0]*s1cs2[0]+ s1cs2[1]*s1cs2[1]+ s1cs2[2]*s1cs2[2]);
         weigth = 2.0/slength;
         mulfactor = weigth/slength;

#if 1
        /* mulfactor is a normalization */
        /* define mulfactor as 1/distance to the center of the opposite edge */
         mid_edge_0_x = 0.5*(vertex_list[v1-1]->x+vertex_list[v2-1]->x);
         mid_edge_0_y = 0.5*(vertex_list[v1-1]->y+vertex_list[v2-1]->y);
         mid_edge_0_z = 0.5*(vertex_list[v1-1]->z+vertex_list[v2-1]->z);

         mid_edge_1_x = 0.5*(vertex_list[v0-1]->x+vertex_list[v2-1]->x);
         mid_edge_1_y = 0.5*(vertex_list[v0-1]->y+vertex_list[v2-1]->y);
         mid_edge_1_z = 0.5*(vertex_list[v0-1]->z+vertex_list[v2-1]->z);

         mid_edge_2_x = 0.5*(vertex_list[v1-1]->x+vertex_list[v0-1]->x);
         mid_edge_2_y = 0.5*(vertex_list[v1-1]->y+vertex_list[v0-1]->y);
         mid_edge_2_z = 0.5*(vertex_list[v1-1]->z+vertex_list[v0-1]->z);

         norm_v0 = 1.0/sqrt((mid_edge_1_x-vertex_list[v0-1]->x)*(mid_edge_1_x-vertex_list[v0-1]->x)+
                            (mid_edge_1_y-vertex_list[v0-1]->y)*(mid_edge_1_y-vertex_list[v0-1]->y)+
                            (mid_edge_1_z-vertex_list[v0-1]->z)*(mid_edge_1_z-vertex_list[v0-1]->z));


         norm_v1 = 1.0/sqrt((mid_edge_2_x-vertex_list[v1-1]->x)*(mid_edge_2_x-vertex_list[v1-1]->x)+
                            (mid_edge_2_y-vertex_list[v1-1]->y)*(mid_edge_2_y-vertex_list[v1-1]->y)+
                            (mid_edge_2_z-vertex_list[v1-1]->z)*(mid_edge_2_z-vertex_list[v1-1]->z));

         norm_v2 = 1.0/sqrt((mid_edge_0_x-vertex_list[v2-1]->x)*(mid_edge_0_x-vertex_list[v2-1]->x)+
                            (mid_edge_0_y-vertex_list[v2-1]->y)*(mid_edge_0_y-vertex_list[v2-1]->y)+
                            (mid_edge_0_z-vertex_list[v2-1]->z)*(mid_edge_0_z-vertex_list[v2-1]->z));


        vertex_list[v0-1]->nx  += (s1cs2[0]*norm_v0);
        vertex_list[v0-1]->ny  += (s1cs2[1]*norm_v0);
        vertex_list[v0-1]->nz  += (s1cs2[2]*norm_v0);

        vertex_list[v1-1]->nx  += (s1cs2[0]*norm_v1);
        vertex_list[v1-1]->ny  += (s1cs2[1]*norm_v1);
        vertex_list[v1-1]->nz  += (s1cs2[2]*norm_v1);

        vertex_list[v2-1]->nx  += (s1cs2[0]*norm_v2);
        vertex_list[v2-1]->ny  += (s1cs2[1]*norm_v2);
        vertex_list[v2-1]->nz  += (s1cs2[2]*norm_v2);

#else



        vertex_list[v0-1]->nx  += (s1cs2[0]*mulfactor);
        vertex_list[v0-1]->ny  += (s1cs2[1]*mulfactor);
        vertex_list[v0-1]->nz  += (s1cs2[2]*mulfactor);

        vertex_list[v1-1]->nx  += (s1cs2[0]*mulfactor);
        vertex_list[v1-1]->ny  += (s1cs2[1]*mulfactor);
        vertex_list[v1-1]->nz  += (s1cs2[2]*mulfactor);

        vertex_list[v2-1]->nx  += (s1cs2[0]*mulfactor);
        vertex_list[v2-1]->ny  += (s1cs2[1]*mulfactor);
        vertex_list[v2-1]->nz  += (s1cs2[2]*mulfactor);

#endif

        vertex_list[v0-1]->nf  += 1;
        vertex_list[v1-1]->nf  += 1;
        vertex_list[v2-1]->nf  += 1;
     }

    boundary_list[BC_type-1]->Add_BFace(vert_num_,v0-1,v1-1,v2-1,v3-1);

  }

  int VN = Get_number_of_vertices();

  for(int i=0; i<VN; i++)
  {
    if(vertex_list[i]->nf >0.0)
    {

        slength = sqrt( vertex_list[i]->nx*vertex_list[i]->nx  +  vertex_list[i]->ny*vertex_list[i]->ny + vertex_list[i]->nz*vertex_list[i]->nz);

        vertex_list[i]->nx = (vertex_list[i]->nx)/slength;
        vertex_list[i]->ny = (vertex_list[i]->ny)/slength;
        vertex_list[i]->nz = (vertex_list[i]->nz)/slength;
         
    }
  } 


  printf(" done with boundary list \n ");

};

void Mesh::Find_match_element_face(int &eid_, int &fid_, const int v0, const int v1, const int v2, const int v3){
    //search for an element and face with given vertices find an
    //element(need elmts_vertisin)
  int eid,cnt;
  int e=NONE;
  Int_Set_Iter iter1;

  for(iter1=vertex_list[v0]->elmts_vertisin.begin();iter1!=vertex_list[v0]->elmts_vertisin.end();++iter1){
    eid=*iter1;
    cnt=1;
    
    if(eid != eid_){
      if(find(vertex_list[v1]->elmts_vertisin.begin(),vertex_list[v1]->elmts_vertisin.end(),eid) != vertex_list[v1]->elmts_vertisin.end()) ++cnt;
      if(find(vertex_list[v2]->elmts_vertisin.begin(),vertex_list[v2]->elmts_vertisin.end(),eid) != vertex_list[v2]->elmts_vertisin.end()) ++cnt;
      if(v3 >= 0)
	if(find(vertex_list[v3]->elmts_vertisin.begin(),vertex_list[v3]->elmts_vertisin.end(),eid) != vertex_list[v3]->elmts_vertisin.end()) ++cnt;

      
      if(v3 >= 0){
	if(cnt == 4) e=eid;
      }
      else{
	if(cnt == 3) e=eid;
      }
    }
  }
  
  eid_=e;//the element # is found
  
  if(eid_ == NONE){
    fid_ = NONE;
    return;
  }
  else{//call the Element:Match_face to get the face id
    fid_=element_list[eid_]->Find_match_face(v0,v1,v2,v3);
    if(fid_ == NONE){
      fprintf(stderr,"Smth is wrong in Mesh::Find_match_element_face(...)\n");
      exit(-1);
    };
  };
};

int Mesh::Find_match_vertex(const double x_, const double y_, const double z_){
  int vid;
  for(vid=0;vid<Get_number_of_vertices();++vid){
    if (fabs(x_ - vertex_list[vid]->x) < EPSILON)
      if (fabs(y_ - vertex_list[vid]->y) < EPSILON)
	if (fabs(z_ - vertex_list[vid]->z) < EPSILON)
	  return vid;
  }
  return NONE;
};

int Mesh::curvedfaces_check(){
  int flag=1;//1 is OK, i.e. no 2 curved face
  int vid,eid,j,bid,bfid,v0,v1,v2,v3,eid_,fid_;
  int* num_curved_sides = new int[Get_number_of_elements()];

  for(vid=0;vid<Get_number_of_vertices();++vid)
    vertex_list[vid]->Clear_elmts_vertisin();
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(j=0;j<element_list[eid]->Get_number_of_vertices();++j){
      vid = element_list[eid]->Get_vert_id(j);
      vertex_list[vid]->Add_elmts_vertisin(eid);
    }
  }
  
  for(eid=0;eid<Get_number_of_elements();++eid) num_curved_sides[eid]=0;
  
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1){
      for(bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid_=NONE; fid_=NONE;
	Find_match_element_face(eid_,fid_,v0,v1,v2,v3);
	if(eid_ != NONE) ++num_curved_sides[eid_];
             if(num_curved_sides[eid_]>1)
            {
              fprintf(stderr," V0 = %d  V1 = %d  V2 = %d  V3 = %d , Bfaces = %d\n", v0, v1, v2, v3, boundary_list[bid]->Get_number_of_bfaces());
              fprintf(stderr,"eid_ = %d, fid_ = %d\n", eid_, fid_);

            }
        

      }
    }
  }
  
  for(eid=0;eid<Get_number_of_elements();++eid) 
    if(num_curved_sides[eid] > 1){
      flag = 0;
      fprintf(stderr,"WARNING: %ld curved faces in element %ld\n",num_curved_sides[eid],eid+1);
    }
  free(num_curved_sides);
  for(vid=0;vid<Get_number_of_vertices();++vid)
    vertex_list[vid]->Clear_elmts_vertisin();
  return flag;
};

void Mesh::curvedfaces_fix(){
  //the purpose of this function is to get rid of tets with 2 curved faces
  //nektar can handle only one curved face per element
  //Done for tets only!!
  //Idea: split every element in edge shell into 2
  int vid,eid,fid,i,bid,bfid,v0,v1,v2,v3,e0,e1,flag,v_new,eid_;
  double x_new,y_new,z_new;
  int orig_number_of_elements=Get_number_of_elements();
  int* num_curved_sides = new int[Get_number_of_elements()];
  int** faceverts = (int**) calloc(orig_number_of_elements,sizeof(int*));
  for(i=0;i<orig_number_of_elements;++i)
    faceverts[i] = new int[6];
  Int_Set shell;
  Int_Set_Iter shell_iter = shell.begin();

  for(vid=0;vid<Get_number_of_vertices();++vid)
    vertex_list[vid]->Clear_elmts_vertisin();
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(i=0;i<element_list[eid]->Get_number_of_vertices();++i){
      vid = element_list[eid]->Get_vert_id(i);
      vertex_list[vid]->Add_elmts_vertisin(eid);
    }
  }
  
  for(eid=0;eid<Get_number_of_elements();++eid) num_curved_sides[eid]=0;
  //get number of curved faces for each element
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1){
      for(bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	if(eid != NONE){
	  if(num_curved_sides[eid] < 2){
	    faceverts[eid][3*num_curved_sides[eid]+0]=v0;
	    faceverts[eid][3*num_curved_sides[eid]+1]=v1;
	    faceverts[eid][3*num_curved_sides[eid]+2]=v2;
	  } 
	  ++num_curved_sides[eid];
	}
      }
    }
  }
  //find an edge to split for each element
  for(eid=0;eid<orig_number_of_elements;++eid){
    if(num_curved_sides[eid] == 2){
      shell.clear();
      for(i=0;i<3;++i)
	if((faceverts[eid][i] != faceverts[eid][3])
	   &&(faceverts[eid][i] != faceverts[eid][4])
	   &&(faceverts[eid][i] != faceverts[eid][5]))
	  e0=faceverts[eid][i];
      for(i=3;i<6;++i)
	if((faceverts[eid][i] != faceverts[eid][0])
	   &&(faceverts[eid][i] != faceverts[eid][1])
	   &&(faceverts[eid][i] != faceverts[eid][2]))
	  e1=faceverts[eid][i];
      //e0 and e1 are the ends of the edge
      //build a shell around the (e0,e1) edge
      set_intersection(vertex_list[e0]->elmts_vertisin.begin(),
		       vertex_list[e0]->elmts_vertisin.end(),
		       vertex_list[e1]->elmts_vertisin.begin(),
		       vertex_list[e1]->elmts_vertisin.end(),
		       inserter(shell,shell_iter));
      /*
      fprintf(stderr,"eid:: %ld  ==> ",eid);
      fprintf(stderr,"e0 = %ld, e1 = %ld  --> ",e0,e1);
      for(shell_iter=shell.begin();shell_iter!=shell.end();++shell_iter)
      fprintf(stderr," %ld ",*shell_iter);
      fprintf(stderr,"\n");
      */
      //check that all elems in the shell are tetrahedras
      flag=1;
      for(shell_iter=shell.begin();shell_iter!=shell.end();++shell_iter)
	if((element_list[*shell_iter]->type != TET)
	   ||(num_curved_sides[*shell_iter]==-1))
	  flag=0;
      
      if(flag != 0){
	//get coords of new vertex in the middle of (e0,e1) edge
	x_new = 0.5*(vertex_list[e0]->x+vertex_list[e1]->x);
	y_new = 0.5*(vertex_list[e0]->y+vertex_list[e1]->y);
	z_new = 0.5*(vertex_list[e0]->z+vertex_list[e1]->z);
	Add_vertex(x_new,y_new,z_new);
	v_new=Get_number_of_vertices()-1;//new vertex is the last one
	for(shell_iter=shell.begin();shell_iter!=shell.end();++shell_iter){
	  //split elems in the shell
	  eid_=*shell_iter;
	  v0=element_list[eid_]->Get_vert_id(0);
	  v1=element_list[eid_]->Get_vert_id(1);
	  v2=element_list[eid_]->Get_vert_id(2);
	  v3=element_list[eid_]->Get_vert_id(3);
	  for(i=0;i<4;++i)
	    if(element_list[eid_]->Get_vert_id(i) == e0)
	      element_list[eid_]->Set_vert_id(i,v_new);
	  num_curved_sides[eid_]=-1;//don't want to split an element twice
	  Add_tet(v0,v1,v2,v3);//add new element
	  eid_=Get_number_of_elements()-1;//new element is the last one
	  for(i=0;i<4;++i)
	    if(element_list[eid_]->Get_vert_id(i) == e1)
	      element_list[eid_]->Set_vert_id(i,v_new);
	}
      }
      else{
	//fprintf(stderr,"There are NON-tetrahedras in the shell\n");
	//may want to split only 1 element(if this elem is tetrahedra)??
	//can be useful for hybrid meshes
      }
    }
  }

  shell.clear();
  free(num_curved_sides);
  for(i=0;i<orig_number_of_elements;++i) 
    free(faceverts[i]);
  free(faceverts);
  for(vid=0;vid<Get_number_of_vertices();++vid)
    vertex_list[vid]->Clear_elmts_vertisin();
};


void Mesh::Connect(){
  int eid,fid,j,vid,v0,v1,v2,v3,eid_,fid_;
  fprintf(stderr,"Doing conn. job ...");

  for(vid=0;vid<Get_number_of_vertices();++vid)
    vertex_list[vid]->Clear_elmts_vertisin();
    
  for(eid=0;eid<Get_number_of_elements();++eid){
    //if(!(eid%100)) fprintf(stderr, ".");
    for(j=0;j<element_list[eid]->Get_number_of_vertices();++j){
      vid = element_list[eid]->Get_vert_id(j);
      vertex_list[vid]->Add_elmts_vertisin(eid);
    }
  }
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    //if(!(eid%100)) fprintf(stderr, ".");
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      element_list[eid]->Get_face_ids(fid,v0,v1,v2,v3);
      eid_=eid; fid_=fid;
      Find_match_element_face(eid_,fid_,v0,v1,v2,v3);
      element_list[eid]->conn_element[fid] = eid_;
      element_list[eid]->conn_face[fid] = fid_;
    }
  }
  fprintf(stderr, "\n");
};

void Mesh::Set_Boundary_Conditions(){
  int bid,bfid,v0,v1,v2,v3,eid_,fid_,eid,fid;
  double dx,dy,dz;
  char btype;
  fprintf(stderr,"Setting boundary conditions ... \n");
  for(bid=0;bid<Get_number_of_bcs();++bid){
    for(bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
      boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
      eid_=NONE; fid_=NONE;
      Find_match_element_face(eid_,fid_,v0,v1,v2,v3);
      btype=boundary_list[bid]->Get_type();
      if((btype != 'X') && (btype != 'x') &&
	 (btype != 'Y') && (btype != 'y') &&
	 (btype != 'Z') && (btype != 'z')) {
	element_list[eid_]->conn_element[fid_] = -(bid+1);//use this for BC
	element_list[eid_]->conn_face[fid_] = -(bid+1);//use this for BC
      }
      else{//we have periodic bc
	dx = dy = dz = 0.0;
	switch(btype){
	case 'x': case 'X':
	  dx = periodic_info->xpmax - periodic_info->xpmin;
	  if(fabs(periodic_info->xpmax - vertex_list[v0]->x) < EPSILON) dx = -dx;
	  break;
	case 'y': case 'Y':
	  dy = periodic_info->ypmax - periodic_info->ypmin;
	  if(fabs(periodic_info->ypmax - vertex_list[v0]->y) < EPSILON) dy = -dy;
	  break;
	case 'z': case 'Z':
	  dz = periodic_info->zpmax - periodic_info->zpmin;
	  if(fabs(periodic_info->zpmax - vertex_list[v0]->z) < EPSILON) dz = -dz;
	  break;
	}
	v0 = Find_match_vertex(vertex_list[v0]->x+dx,vertex_list[v0]->y+dy,vertex_list[v0]->z+dz);
	v1 = Find_match_vertex(vertex_list[v1]->x+dx,vertex_list[v1]->y+dy,vertex_list[v1]->z+dz);
	v2 = Find_match_vertex(vertex_list[v2]->x+dx,vertex_list[v2]->y+dy,vertex_list[v2]->z+dz);
	if(v3>=0) 
	  v3 = Find_match_vertex(vertex_list[v3]->x+dx,vertex_list[v3]->y+dy,vertex_list[v3]->z+dz);
	if(v0<0 || v1<0 || v2<0){
	  fprintf(stderr,"\n\nERROR in periodic info: Can't find periodic vertex!\n");
	  fprintf(stderr,"Make sure that your mesh is periodic, check periodic info.\n");
	  exit(-1);
	}
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	element_list[eid_]->conn_element[fid_] = eid;
	element_list[eid_]->conn_face[fid_] = fid;
	if((eid<0)||(fid<0)){
	  fprintf(stderr,"\n\nERROR in periodic info: Can't find periodic element!\n");
	  fprintf(stderr,"Make sure that your mesh is periodic, check periodic info.\n");
	  exit(-1);
	}
      }
    }
  }
};

void Mesh::Dump_sum_rea(){
  int eid,vid,fid,bid;
  FILE *reaout;
  char bufa[BUFSIZ];
  
  reaout = fopen("sum.rea", "w");
  
  //fprintf(stderr, "Writing sum.rea file\n");
  
  fprintf(reaout, "****** PARAMETERS *****\n");
  fprintf(reaout, " GRIDGEN 3D -> NEKTAR (nekscal2d -i -u \"x+y+z\" -n4 sum.rea)\n");
  fprintf(reaout, " 3 DIMENSIONAL RUN\n");
  fprintf(reaout, " 5 PARAMETERS FOLLOW\n");
  fprintf(reaout, " 0.001  DT\n"); 
  fprintf(reaout, " 0000. NSTEPS\n");
  fprintf(reaout, " 2. INTYPE\n");
  fprintf(reaout, " 1. EQTYPE\n");
  fprintf(reaout, " 1. LAMBDA\n");
  fprintf(reaout, " 0  Lines of passive scalar data follows\n");
  fprintf(reaout, " 0  LOGICAL SWITCHES FOLLOW\n");
  fprintf(reaout, "Dummy line from old nekton file\n"); 
  fprintf(reaout, "**MESH DATA** x,y,z, values of vertices 1,2,3,4.\n");  
  fprintf(reaout, "%d 	3	1   NEL NDIM NLEVEL\n", Get_number_of_elements());
  
  //fprintf(stderr, "Writing coords\n");
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    switch(element_list[eid]->type){
    case TET: fprintf(reaout, "Element %d Tet\n", eid+1); break;
    case HEX: fprintf(reaout, "Element %d Hex\n", eid+1); break;
    case PRI: fprintf(reaout, "Element %d Prism\n", eid+1); break;
    case PYR: fprintf(reaout, "Element %d Pyr\n", eid+1); break;
    default:  fprintf(stderr,"Unknown type of element!!"); break;
    }
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->x);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->y);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->z);
    fprintf(reaout, "\n");	      
  }

  fprintf(reaout, "***** CURVED SIDE DATA ***** \n");
  FILE *infoin;
  int Nunk,i;
  sprintf(bufa, "%s.info", strtok(title, "\n"));
  if((infoin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    exit(-1);
  }  
  rewind(infoin);
  findSection ("*** CURVED SIDE DATA ***",bufa,infoin);
  fgets(bufa,BUFSIZ, infoin);
  sscanf(bufa, "%d", &Nunk);
  fprintf(reaout, "%d Number of curve types\n", Nunk);
  for(i=0;i<Nunk;++i){
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa)); 
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa));
  }
  fclose(infoin);

  int numcur=0;
  int v0,v1,v2,v3;
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      numcur = numcur + boundary_list[bid]->Get_number_of_bfaces();
  }
  fprintf(reaout, "%d Curved sides follow\n", numcur);
#if 1
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      for(int bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	if(eid==NONE || fid==NONE) 
	  fprintf(stderr,"smth is wrong in Mesh::Dump_walls()!!\n");
	fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
      }
  }
#else
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid < 0){
	bid=-bid-1;//get real boundary condition id
	if(boundary_list[bid]->Get_curved() == 1)
	  fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
      }
    }
  }
#endif  

  fprintf(reaout, "***** BOUNDARY CONDITIONS ***** \n");
  fprintf(reaout, "***** FLUID BOUNDARY CONDITIONS ***** \n");
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid >= 0){
	fprintf(reaout, "E %d %d %d %d\n", eid+1,fid+1,bid+1,element_list[eid]->conn_face[fid]+1);
      }
      else{//we have boundary conditions here
	fprintf(reaout, "v %d %d 0. 0. 0.\n", eid+1, fid+1);
	fprintf(reaout, "u = x+y+z\n");
      }
    }
  }
 
  fprintf(reaout, "***** NO THERMAL BOUNDARY CONDITIONS *****\n");
  fprintf(reaout, "0        INITIAL CONDITIONS *****\n");
  fprintf(reaout, "***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q\n");
  fprintf(reaout, "1                 Lines of Drive force data follow\n");
  fprintf(reaout, "u = -LAMBDA*(x+y+z)\n");
  fprintf(reaout, "***** Variable Property Data ***** Overrrides Parameter data.\n");
  fprintf(reaout, " 1 Lines follow.\n");
  fprintf(reaout, " 0 PACKETS OF DATA FOLLOW\n");
  fprintf(reaout, "***** HISTORY AND INTEGRAL DATA *****\n");
  fprintf(reaout, " 0   POINTS.  Hcode, I,J,H,IEL\n");
  fprintf(reaout, " ***** OUTPUT FIELD SPECIFICATION *****\n");
  fprintf(reaout, "  0 SPECIFICATIONS FOLLOW\n");
  
  fclose(reaout);
};

void Mesh::Dump_prod_rea(){
  int eid,vid,fid,bid;
  FILE *reaout;
  char bufa[BUFSIZ];
  
  reaout = fopen("prod.rea", "w");
  
  //fprintf(stderr, "Writing prod.rea file\n");
  
  fprintf(reaout, "****** PARAMETERS *****\n");
  fprintf(reaout, " GRIDGEN 3D -> NEKTAR (nekscal2d -i -u \"x*y*z\" -n4 prod.rea)\n");
  fprintf(reaout, " 3 DIMENSIONAL RUN\n");
  fprintf(reaout, " 5 PARAMETERS FOLLOW\n");
  fprintf(reaout, " 0.001  DT\n"); 
  fprintf(reaout, " 0000. NSTEPS\n");
  fprintf(reaout, " 2. INTYPE\n");
  fprintf(reaout, " 1. EQTYPE\n");
  fprintf(reaout, " 1. LAMBDA\n");
  fprintf(reaout, " 0  Lines of passive scalar data follows\n");
  fprintf(reaout, " 0  LOGICAL SWITCHES FOLLOW\n");
  fprintf(reaout, "Dummy line from old nekton file\n"); 
  fprintf(reaout, "**MESH DATA** x,y,z, values of vertices 1,2,3,4.\n");  
  fprintf(reaout, "%d 	3	1   NEL NDIM NLEVEL\n", Get_number_of_elements());
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    switch(element_list[eid]->type){
    case TET: fprintf(reaout, "Element %d Tet\n", eid+1); break;
    case HEX: fprintf(reaout, "Element %d Hex\n", eid+1); break;
    case PRI: fprintf(reaout, "Element %d Prism\n", eid+1); break;
    case PYR: fprintf(reaout, "Element %d Pyr\n", eid+1); break;
    default:  fprintf(stderr,"Unknown type of element!!"); exit(-1); break;
    }
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->x);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->y);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",vertex_list[element_list[eid]->Get_vert_id(vid)]->z);
    fprintf(reaout, "\n");	      
  }

  fprintf(reaout, "***** CURVED SIDE DATA ***** \n");
  FILE *infoin;
  int Nunk,i;
  sprintf(bufa, "%s.info", strtok(title, "\n"));
  if((infoin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    exit(-1);
  }  
  rewind(infoin);
  findSection ("*** CURVED SIDE DATA ***",bufa,infoin);
  fgets(bufa,BUFSIZ, infoin);
  sscanf(bufa, "%d", &Nunk);
  fprintf(reaout, "%d Number of curve types\n", Nunk);
  for(i=0;i<Nunk;++i){
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa)); 
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa));
  }
  fclose(infoin);

  int numcur=0;
  int v0,v1,v2,v3;
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      numcur = numcur + boundary_list[bid]->Get_number_of_bfaces();
  }
  fprintf(reaout, "%d Curved sides follow\n", numcur);
#if 1
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      for(int bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	if(eid==NONE || fid==NONE) 
	  fprintf(stderr,"smth is wrong in Mesh::Dump_walls()!!\n");
	fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
      }
  }
#else
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid < 0){
	bid=-bid-1;//get real boundary condition id
	if(boundary_list[bid]->Get_curved() == 1)
	  fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
      }
    }
  }
#endif  
  
  fprintf(reaout, "***** BOUNDARY CONDITIONS ***** \n");
  fprintf(reaout, "***** FLUID BOUNDARY CONDITIONS ***** \n");
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid >= 0){
	fprintf(reaout, "E %d %d %d %d\n", eid+1,fid+1,bid+1,element_list[eid]->conn_face[fid]+1);
      }
      else{//we have boundary conditions here
	fprintf(reaout, "v %d %d 0. 0. 0.\n", eid+1, fid+1);
	fprintf(reaout, "u = x*y*z\n");
      }
    }
  }
 
  fprintf(reaout, "***** NO THERMAL BOUNDARY CONDITIONS *****\n");
  fprintf(reaout, "0        INITIAL CONDITIONS *****\n");
  fprintf(reaout, "***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q\n");
  fprintf(reaout, "1                 Lines of Drive force data follow\n");
  fprintf(reaout, "u = -LAMBDA*(x*y*z)\n");
  fprintf(reaout, "***** Variable Property Data ***** Overrrides Parameter data.\n");
  fprintf(reaout, " 1 Lines follow.\n");
  fprintf(reaout, " 0 PACKETS OF DATA FOLLOW\n");
  fprintf(reaout, "***** HISTORY AND INTEGRAL DATA *****\n");
  fprintf(reaout, " 0   POINTS.  Hcode, I,J,H,IEL\n");
  fprintf(reaout, " ***** OUTPUT FIELD SPECIFICATION *****\n");
  fprintf(reaout, "  0 SPECIFICATIONS FOLLOW\n");
  
  fclose(reaout);
};

void Mesh::Dump_rea(){
  int eid,vid,fid,bid,bfid,i;
  FILE *reaout;
  char bufa[BUFSIZ];
  FILE *infoin;
  int Nunk;
  char btype;
  char *inlineparam;
  
  sprintf(bufa, "%s.rea", strtok(title, "\n"));
  reaout = fopen(bufa, "w");
  
  fprintf(stderr, "Writing .rea files\n");
  
  sprintf(bufa, "%s.info", strtok(title, "\n"));
  if((infoin = fopen(bufa, "r")) == (FILE*) NULL) {
    fprintf(stderr, "ATTENTION: unable to open the input file -- %s.\n",bufa);
    exit(-1);
  } 
  rewind(infoin);
  findSection ("*** HEAD ***",bufa,infoin);
  if(!strstr(bufa,"*** HEAD ***")){
    fprintf(reaout, "****** PARAMETERS *****\n");
    fprintf(reaout, "GRIDGEN 3D -> NEKTAR\n");
    fprintf(reaout, "3 DIMENSIONAL RUN\n");
    fprintf(reaout, "6 PARAMETERS FOLLOW\n");
    fprintf(reaout, "0.01     KINVIS\n"); 
    fprintf(reaout, "4        MODES\n"); 
    fprintf(reaout, "0.001    DT\n"); 
    fprintf(reaout, "10000    NSTEPS\n");
    fprintf(reaout, "1        EQTYPE\n");
    fprintf(reaout, "2        INTYPE\n");
    fprintf(reaout, "0  Lines of passive scalar data follows\n");
    fprintf(reaout, "0  LOGICAL SWITCHES FOLLOW\n");
    fprintf(reaout, "Dummy line from old nekton file\n"); 
  }
  else{
    fgets(bufa,BUFSIZ, infoin);
    sscanf(bufa, "%d", &Nunk);
    for(i=0;i<Nunk;++i){
      fgets(bufa,BUFSIZ, infoin);
      fprintf(reaout,strdup(bufa)); 
    }
  }

  fprintf(reaout, "**MESH DATA** x,y,z, values of vertices 1,2,3,4.\n");
  fprintf(reaout, "%d 	3	1   NEL NDIM NLEVEL\n", Get_number_of_elements());
  fprintf(stderr, "Writing node coords\n");

  double scal=1.0;
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    switch(element_list[eid]->type){
    case TET: fprintf(reaout, "Element %d Tet\n", eid+1); break;
    case HEX: fprintf(reaout, "Element %d Hex\n", eid+1); break;
    case PRI: fprintf(reaout, "Element %d Prism\n", eid+1); break;
    case PYR: fprintf(reaout, "Element %d Pyr\n", eid+1); break;
    default:  fprintf(stderr,"Unknown type of element!!"); break;
    }
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",scal*vertex_list[element_list[eid]->Get_vert_id(vid)]->x);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",scal*vertex_list[element_list[eid]->Get_vert_id(vid)]->y);
    fprintf(reaout, "\n");
    for(vid=0;vid<element_list[eid]->Get_number_of_vertices();++vid)
      fprintf(reaout, "%lf ",scal*vertex_list[element_list[eid]->Get_vert_id(vid)]->z);
    fprintf(reaout, "\n");	      
  }

  fprintf(stderr, "Writing curved side data\n");
  fprintf(reaout, "***** CURVED SIDE DATA ***** \n");

  rewind(infoin);
  findSection ("*** CURVED SIDE DATA ***",bufa,infoin);
  fgets(bufa,BUFSIZ, infoin);
  sscanf(bufa, "%d", &Nunk);
  fprintf(reaout, "%d Number of curve types\n", Nunk);
  for(i=0;i<Nunk;++i){
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa)); 
    fgets(bufa,BUFSIZ, infoin);
    fprintf(reaout,strdup(bufa));
  }

  int numcur=0;
  int v0,v1,v2,v3;
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      numcur = numcur + boundary_list[bid]->Get_number_of_bfaces();
  }
  fprintf(reaout, "%d Curved sides follow\n", numcur);
#if 0//ordered by element number
  for(bid=0;bid<Get_number_of_bcs();++bid){
    if(boundary_list[bid]->Get_curved() == 1)
      for(bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	if(eid==NONE || fid==NONE) 
	  fprintf(stderr,"smth is wrong in curved info!!\n");
	fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
      }
  }
#else//ordered by BCs number
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid < 0){
	bid=-bid-1;//get real boundary condition id
	if(boundary_list[bid]->Get_curved() == 1)
        {
	  fprintf(reaout,"%ld %ld %c\n",fid+1,eid+1,boundary_list[bid]->Get_curve_type());
        

         btype  =  boundary_list[bid]->Get_type();
         if( (btype =='W') )
         {
            inlineparam = boundary_list[bid]->Get_inline_params();
            if(strstr(inlineparam,"SI"))
            {
               element_list[eid]->Get_face_ids(fid,v0,v1,v2,v3);
               fprintf(reaout, "%lf %lf %lf\n", vertex_list[v0]-> nx, vertex_list[v1]->nx , vertex_list[v2]->nx);
               fprintf(reaout, "%lf %lf %lf\n", vertex_list[v0]-> ny, vertex_list[v1]->ny , vertex_list[v2]->ny);
               fprintf(reaout, "%lf %lf %lf\n", vertex_list[v0]-> nz, vertex_list[v1]->nz , vertex_list[v2]->nz);

            }
         }
       }
      }
    }
  }
#endif  

  fprintf(stderr, "Writing boundary conditions\n");
  fprintf(reaout, "***** BOUNDARY CONDITIONS ***** \n");
  fprintf(reaout, "***** FLUID BOUNDARY CONDITIONS ***** \n");
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(fid=0;fid<element_list[eid]->Get_number_of_faces();++fid){
      bid=element_list[eid]->conn_element[fid];
      if(bid >= 0){
	fprintf(reaout, "E %d %d %d %d\n", eid+1,fid+1,bid+1,element_list[eid]->conn_face[fid]+1);
      }
      else{//we have -(bid+1) boundary conditions here
	bid=-bid-1;//get real boundary condition id
	fprintf(reaout, "%c %d %d %s",boundary_list[bid]->Get_type(),eid+1,fid+1,boundary_list[bid]->Get_inline_params());
	for(i=0;i<boundary_list[bid]->Get_number_of_lines();++i){
	  fprintf(reaout,strdup(boundary_list[bid]->Get_lines(i)));
	}
      }
    }
  }
 
  rewind(infoin);
  findSection ("*** TAIL ***",bufa,infoin);
  if(!strstr(bufa,"*** TAIL ***")){
    fprintf(reaout, "***** NO THERMAL BOUNDARY CONDITIONS *****\n");
    fprintf(reaout, "4        INITIAL CONDITIONS *****\n");
    fprintf(reaout, "Given\n");
    fprintf(reaout, "u = 1.0\n");
    fprintf(reaout, "v = 0.0\n");
    fprintf(reaout, "w = 0.0\n");
    fprintf(reaout, "***** DRIVE FORCE DATA ***** PRESSURE GRAD, FLOW, Q\n");
    fprintf(reaout, "0                 Lines of Drive force data follow\n");
    fprintf(reaout, "***** Variable Property Data ***** Overrrides Parameter data.\n");
    fprintf(reaout, "1 Lines follow.\n");
    fprintf(reaout, "0 PACKETS OF DATA FOLLOW\n");
    fprintf(reaout, "***** HISTORY AND INTEGRAL DATA *****\n");
    fprintf(reaout, "1   POINTS.  Hcode, I,J,H,IEL\n");
    fprintf(reaout, "UVWP  H  1 0 0 10\n");
    fprintf(reaout, "***** OUTPUT FIELD SPECIFICATION *****\n");
    fprintf(reaout, "0 SPECIFICATIONS FOLLOW\n");
  }
  else{
    fgets(bufa,BUFSIZ, infoin);
    sscanf(bufa, "%d", &Nunk);
    for(i=0;i<Nunk;++i){
      fgets(bufa,BUFSIZ, infoin);
      fprintf(reaout,strdup(bufa)); 
    }
  }

  fclose(infoin);
  fclose(reaout);
};

void Mesh::Dump_plt(){
  FILE *tecout;
  char bufa[BUFSIZ];
  int eid,k,i,pr;
  int vid;
  int tetvid[] = {0,3,1,3,2,3,2,3};
  int hexvid[] = {0,1,2,3,4,5,6,7};
  int privid[] = {0,1,2,3,4,4,5,5};
  int pyrvid[] = {0,4,1,4,2,4,3,4};//need to test pyrs
  
  sprintf(bufa, "%s.plt", strtok(title, "\n"));
  tecout = fopen(bufa, "w");
  
  fprintf(tecout, "TITLE = \"  \"\n");
  fprintf(tecout, "VARIABLES = \"X\", \"Y\", \"Z\", \"FUNCTION\"\n");
  fprintf(tecout, "ZONE  N=%d, E=%d, F=FEPOINT, ET=BRICK \n",
	  Get_number_of_elements()*8,Get_number_of_elements());
  
  for(eid=0;eid<Get_number_of_elements();++eid){
    for(i=0;i<8;++i){
      switch(element_list[eid]->type){
      case TET: vid=tetvid[i]; break;
      case HEX: vid=hexvid[i]; break;
      case PRI: vid=privid[i]; break;
      case PYR: vid=pyrvid[i]; break;
      default:  fprintf(stderr,"Unknown type of element!!"); break;
      }
      fprintf(tecout, "%lf,%lf,%lf,%d\n", 
	      vertex_list[element_list[eid]->Get_vert_id(vid)]->x,
	      vertex_list[element_list[eid]->Get_vert_id(vid)]->y,
	      vertex_list[element_list[eid]->Get_vert_id(vid)]->z,
	      0);
    }
  }
  for (k=0;k<Get_number_of_elements();++k){
    pr=8*k;
    for(i=1;i<=8;++i) fprintf(tecout, "%d, ",pr+i);
    fprintf(tecout, "\n ");
  }
  fclose(tecout);
}

void Mesh::Dump_walls(){
  FILE *wallsout;
  char bufa[BUFSIZ];
  int bid,bfid,v0,v1,v2,v3,eid,fid,numwalls;
  char btype;
  
  sprintf(bufa, "%s.walls", strtok(title, "\n"));
  wallsout = fopen(bufa, "w");
  
  numwalls=0;
  fprintf(wallsout,"Body\n");
  for(bid=0;bid<Get_number_of_bcs();++bid){
    btype=boundary_list[bid]->Get_type();
    if(btype == 'W')
      numwalls = numwalls + boundary_list[bid]->Get_number_of_bfaces();
  }
  fprintf(wallsout,"%ld\n",numwalls);
  
  for(bid=0;bid<Get_number_of_bcs();++bid){
    btype=boundary_list[bid]->Get_type();
    if(btype == 'W')
      for(bfid=0;bfid<boundary_list[bid]->Get_number_of_bfaces();++bfid){
	boundary_list[bid]->bface_list[bfid]->Get_face_ids(v0,v1,v2,v3);
	eid=NONE; fid=NONE;
	Find_match_element_face(eid,fid,v0,v1,v2,v3);
	if(eid==NONE || fid==NONE) 
	  fprintf(stderr,"smth is wrong in Mesh::Dump_walls()!!\n");
	fprintf(wallsout,"%ld  %ld\n",eid+1,fid+1);
      }
  }
  fclose(wallsout);
}

main (int argc, char *argv[]){
  Mesh* mesh = new Mesh();
  mesh->Read_input_files();
#if 1
  char  buf[BUFSIZ];
  if(mesh->curvedfaces_check() == 0){
    fprintf(stderr,"\n\nNektar allows only one curved face per element.\n");
    fprintf(stderr,"Your mesh has elements with more than 2 curved faces! \n");
    fprintf(stderr,"Do you want to split these elements (y or n)?: ");
    fgets(buf,BUFSIZ, stdin);
    if(buf[0]!='n'){
      fprintf(stderr,"Splitting elements ...\n");
      
      int fix_cnt=0;
      while((mesh->curvedfaces_check() != 1)&&(fix_cnt != 6)){
	mesh->curvedfaces_fix();
	++fix_cnt;
      }
      
      if(mesh->curvedfaces_check() == 1){
	fprintf(stderr,"Fixed!\n\n\n");
      }
      else{
	fprintf(stderr,"There are some elements with 2 curved faces left!\n");
	fprintf(stderr,"This algorithm works for tetrahedras only.\n\n\n");
      }
    }
  }
#endif
  mesh->Orient_elements();
  mesh->Connect();
  mesh->Dump_sum_rea();
  mesh->Dump_prod_rea();
  mesh->Set_Boundary_Conditions();
  mesh->Dump_rea();
  mesh->Dump_walls();
  mesh->Dump_plt();
}
