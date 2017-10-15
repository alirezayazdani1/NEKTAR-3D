/*---------------------------------------------------------------------------*
 *                        RCS Information                                    *
 *                                                                           *
 * $Source: 
 * $Revision:
 * $Date: 
 * $Author:
 * $State: 
 *---------------------------------------------------------------------------*/

#include "nektar.h"
#include "nekstruct.h"
#include "pbc_1d.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <veclib.h>


#ifdef PBC_1D

/* external function */
double BoundaryArea(Bndry *Bc, char *label);
double BoundaryAreaInlet(Bndry *Bc, char *label);
void M2Ptransform(int Nmodes, double om, my_dcmplx *MODES, int N, double Time_period, double *f);
void M2Ptransform(int Nmodes, double om, my_dcmplx* MODES , int N, double Time_period, int index_start, int index_end, double* f);
void FilterDFT(int N, int Nmodes, double  *f, double *Acos_sin, int FLAG_averag);
void Convolve(int N, double* f1, double* f2, double* ans, int time_index);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int it);
void ParallelConvolve(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);
void ParallelConvolve_test(int N, double* f1, double* f2, double* ans, int index_start, int index_end, int it);

#define N_MAX_OUTLET 26 

int PBC1D::Create(Bndry *Ubc){

  int i,j;
  Bndry *Bc;
  int LABEL_FLAGS[N_MAX_OUTLET];
  memset(&LABEL_FLAGS[0],'\0',N_MAX_OUTLET*sizeof(int));

  standard_labels[0] = "bcA(x,y,z)\n";
  standard_labels[1] = "bcB(x,y,z)\n";
  standard_labels[2] = "bcC(x,y,z)\n";
  standard_labels[3] = "bcD(x,y,z)\n";
  standard_labels[4] = "bcE(x,y,z)\n";
  standard_labels[5] = "bcF(x,y,z)\n";
  standard_labels[6] = "bcG(x,y,z)\n";
  standard_labels[7] = "bcH(x,y,z)\n";
  standard_labels[8] = "bcI(x,y,z)\n";
  standard_labels[9] = "bcJ(x,y,z)\n";
  standard_labels[10] = "bcK(x,y,z)\n";
  standard_labels[11] = "bcL(x,y,z)\n";
  standard_labels[12] = "bcM(x,y,z)\n";
  standard_labels[13] = "bcN(x,y,z)\n";
  standard_labels[14] = "bcO(x,y,z)\n";
  standard_labels[15] = "bcP(x,y,z)\n";
  standard_labels[16] = "bcQ(x,y,z)\n";
  standard_labels[17] = "bcR(x,y,z)\n";
  standard_labels[18] = "bcS(x,y,z)\n";
  standard_labels[19] = "bcT(x,y,z)\n";
  standard_labels[20] = "bcU(x,y,z)\n";
  standard_labels[21] = "bcV(x,y,z)\n";
  standard_labels[22] = "bcW(x,y,z)\n";
  standard_labels[23] = "bcX(x,y,z)\n";
  standard_labels[24] = "bcY(x,y,z)\n";
  standard_labels[25] = "bcZ(x,y,z)\n";

  Nout = iparam("NOUTLETS");
  Ninl = iparam("NINLETS");
  

  /* STEP 1 check for consistency */

 /* create vector of integers "LABEL_FLAGS" 
    initialize with zero 
    loop over Bondaries look for ALL labels
    if label of outlet was found set corresponding flag to 1;
    sum results from all pertitions if parallel
    if Nout != sum (LABEL_FLAGS) -> setup is incosistent
       error message and exit.
    else
      deallocate  "LABEL_FLAGS" and continue setup
   */
  
   for (i = 0; i < N_MAX_OUTLET; ++i){

     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          LABEL_FLAGS[i] = 1;
          break;
       }
     }
   }

   DO_PARALLEL{
      int *iwork;
      iwork = ivector(0,N_MAX_OUTLET-1);
      memset(iwork,'\0',N_MAX_OUTLET*sizeof(int));
      gisum (LABEL_FLAGS, N_MAX_OUTLET, iwork);
      free(iwork);
   }

   j = 0;
   for (i = 0; i < N_MAX_OUTLET; ++i){
     if (LABEL_FLAGS[i] != 0)
       j++;
   }
   if (j != Nout){
     fprintf(stderr, "Error in setup outlet BC. NOUTLETS in REA file = %d actual outlet count %d \n ",Nout,j);
     exit(-1);
   }

   for (i = 0; i < Nout; ++i){
     if (LABEL_FLAGS[i] == 0){
        fprintf(stderr, "Error in setup outlet BC. NOUTLETS in REA file = %d, outlet No. %d has no label \n",Nout,i);
        exit(-1);
     }
   }

   if (Nout == 0) 
    return 0;


   /* STEP 2 define parameters and allocate memory*/
  double dt = dparam("DELT");
  double t  = dparam("t");

  Nimpedance_modes  =  iparam("NIMPEDMODES");  
  Tperiod           =  dparam("DTIMEPERIOD");
 
  if (Tperiod == 0)
    Tperiod = dt*10.0;

  Nsteps_per_cycle = int ((Tperiod + dt*0.1)/dt)+1;

  j = numnodes();
  if ((Nsteps_per_cycle/j) > 100 && j > 1)
     parallel_convolution = 1;
  else
     parallel_convolution = 0;

  parallel_convolution = 0;

  /* for parallel convolution split "impedance" in Nproc. parts */

   if (parallel_convolution){

     Imp_length_local = (int) (Nsteps_per_cycle/j);
     i = Nsteps_per_cycle - Imp_length_local*j;
     if (mynode() < i)
       Imp_length_local++;

    /* to partition the global vector "impedance" find 
        first index corresponding to local processor */
  
     int *itemp;
     itemp = new int[j*2];
     memset(itemp,'\0',j*2*sizeof(int));
     DO_PARALLEL{
       itemp[mynode()] = Imp_length_local; 
       gisum(itemp,j,&itemp[j]);
     }

     Imp_index_start = 0;
     for (i = 1; i <= mynode(); i++) 
       Imp_index_start += itemp[i-1];

     free(itemp);
  }
  t = t-((int)((t+dt*0.1)/Tperiod))*Tperiod;
     time_step_in_cycle = int ((t+0.1*dt)/dt);

  FR_history = dmatrix(0,Nout-1,0,Nsteps_per_cycle-1);
  if (!parallel_convolution)
    impedance  = dmatrix(0,Nout-1,0,Nsteps_per_cycle-1);
  else
    impedance  = dmatrix(0,Nout-1,0,Imp_length_local);
  
  A_cos_sin  = dmatrix(0,Nout-1,0,Nimpedance_modes*2);
  impedance_modes = new my_dcmplx*[Nout];
  for (i = 0; i < Nout; ++i)
   impedance_modes[i] = new my_dcmplx[Nimpedance_modes*2+1];

  Pressure   = dvector(0,Nout+Ninl);
  Area       = dvector(0,Nout+Ninl);
  memset(Pressure,'\0',(Nout+Ninl)*sizeof(double));
  memset(Area,    '\0',(Nout+Ninl)*sizeof(double));

  for (i = 0; i < Nout; ++i){
    memset(FR_history[i],'\0',Nsteps_per_cycle*sizeof(double));
    if (!parallel_convolution)
      memset(impedance[i], '\0',Nsteps_per_cycle*sizeof(double));
    else
      memset(impedance[i], '\0',(Imp_length_local+1)*sizeof(double));
    memset(A_cos_sin[i], '\0',(Nimpedance_modes*2+1)*sizeof(double));
  }


  /* STEP 3 */
  /*  compute impedance */

  double rad,r_temp; 
  
  for (i = 0; i < Nout; ++i){
     rad =  BoundaryArea(Ubc, standard_labels[i]);

     DO_PARALLEL{
       gdsum(&rad,1,&r_temp);
     }
     Area[i] = rad;
     ROOTONLY
      printf("Area[outlet No. %d]  = %f \n",i,Area[i]);

     /*  radius must be passed in dimensional units! [cm] !!!! */ 
     rad = sqrt(rad/M_PI)*0.1;
     {
       SMALL_VESSEL root_vessel(rad);
       /* set omega  */
       ROOTONLY fprintf(stdout,"Ws = %2.16f \n",root_vessel.Ws);    
       int index = Nimpedance_modes+1;
       for (j = 1; j <= Nimpedance_modes; j++){
           impedance_modes[i][index] = root_vessel.GetZ0(j);
           impedance_modes[i][Nimpedance_modes-j].my_dcmplx_conj(impedance_modes[0][index]);
           ROOTONLY
             fprintf(stdout,"rank %d: impedance_modes[out = %d][mode = %d] = (%e,%e) \n",mynode(),
                             i,index,impedance_modes[i][index].real,impedance_modes[i][index].imag);
           index++;
       }
       impedance_modes[i][Nimpedance_modes] = root_vessel.GetZ0(0);
       ROOTONLY
         fprintf(stdout,"rank %d:impedance_modes[out = %d][mode = %d] = (%e,%e) \n",mynode(),
                               i,Nimpedance_modes,impedance_modes[i][Nimpedance_modes].real,
                                              impedance_modes[i][Nimpedance_modes].imag);

       ROOTONLY
         printf("SMALL_VESSEL - done, r = %f \n", rad );

      if (!parallel_convolution)
         M2Ptransform(Nimpedance_modes,omega_small_atree,impedance_modes[i],
                      Nsteps_per_cycle,2.0*M_PI/omega_small_atree,impedance[i]);
      else
         M2Ptransform(Nimpedance_modes,omega_small_atree,impedance_modes[i],
                      Nsteps_per_cycle,2.0*M_PI/omega_small_atree,
                      Imp_index_start,(Imp_index_start+Imp_length_local),impedance[i]);
    }
  }

  if (mynode() == (numnodes()-1))
    Imp_length_local--;

  /* compute area at inlet */
  for (i = 0; i < Ninl; ++i){
    rad =  BoundaryAreaInlet(Ubc, standard_labels[i]);
    r_temp = 0.0;
    DO_PARALLEL
     gdsum(&rad,1,&r_temp);
    Area[Nout+i] = rad;
    ROOTONLY
      printf("Area - inlet[%d] = %f \n",i,Area[Nout+i]);
  }

  /* compute number of faces with type='O' for each outlet in each partition */
  Nfaces_per_outlet = new int[Nout];
  Nfaces_per_inlet  = new int[Ninl];
  int face_index;

  for (i = 0; i < Nout; ++i){
     Nfaces_per_outlet[i] = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 )
          Nfaces_per_outlet[i]++;
     }
   } 
  

  for (i = 0; i < Ninl; ++i){
     Nfaces_per_inlet[i] = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'V' || Bc->type == 'v') && strcmp(Bc->blabel,standard_labels[i]) == 0 )
          Nfaces_per_inlet[i]++;
     }
   }

 
   /*  for each outlet get store ID's of faces with type == 'O' */
   ID_faces_per_outlet = new int*[Nout];
   for (i = 0; i < Nout; ++i)
     ID_faces_per_outlet[i] = new int[Nfaces_per_outlet[i]]; 

   for (i = 0; i < Nout; ++i){
     j = 0;
     face_index = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'O' || Bc->type == 'o') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          ID_faces_per_outlet[i][face_index] = j;
          face_index++;
       }
       j++;
     }
   }

   ID_faces_per_inlet = new int*[Ninl];
   for (i = 0; i < Ninl; ++i)
     ID_faces_per_inlet[i] = new int[Nfaces_per_inlet[i]];

   for (i = 0; i < Ninl; ++i){
     j = 0;
     face_index = 0;
     for(Bc=Ubc;Bc;Bc=Bc->next){
       if ( (Bc->type == 'V' || Bc->type == 'v') && strcmp(Bc->blabel,standard_labels[i]) == 0 ){
          ID_faces_per_inlet[i][face_index] = j;
          face_index++;
       }
       j++;
     }
   }

   Nnodes_inlet  = new int[Ninl];
   Nnodes_outlet = new int[Nout];

#ifdef PARALLEL
   /* create communicator for outlets */
//   create_comm_BC(Nout, Nfaces_per_outlet);
   create_comm_BC_inlet_outlet(Ninl, Nfaces_per_inlet, Nnodes_inlet, Nout, Nfaces_per_outlet, Nnodes_outlet);
#endif

  return 0;
}

void PBC1D::SetGeofac(Bndry *Ubc, Bndry *Vbc, Bndry *Wbc){
  
   int index;
   Bndry *BcU, *BcV, *BcW; 
  
   for (BcU=Ubc,BcV=Vbc,BcW=Wbc;BcU;BcU = BcU->next,BcV = BcV->next,BcW = BcW->next){
     if ( BcU->type == 'o' || BcU->type == 'O' || BcU->type == 'v' || BcU->type == 'V' ){
               BcU->elmt->Surface_geofac(BcU);
               BcV->elmt->Surface_geofac(BcV);
               BcW->elmt->Surface_geofac(BcW);
     }      
   }
} 
  

void PBC1D::SetRC(char *name){

    /* set-up RC boundary condition */
    R1 = new double[Nout*2];
    C1 = R1+Nout;
    memset(R1,'\0',Nout*2*sizeof(double));
    flowrate_RCR_old = new double[Nout];
    memset(flowrate_RCR_old,'\0',Nout*sizeof(double));
    /* check if RCfile exist if yes read values of R1 and C1,
       otherwise compute R1 and C1, then create the RCfile */

    int i;
    FILE *pRC_File;
    char fname_RC[BUFSIZ];
    sprintf (fname_RC, "%s.RC", name );
    pRC_File = fopen(fname_RC,"r");
    if (pRC_File==NULL){
      ROOTONLY
        pRC_File = fopen(fname_RC,"w");

      for (i = 0; i < Nout; ++i){
        R1[i] = impedance_modes[i][0].real;// /sqrt(Area[i]/M_PI);
        R1[i] *= (0.1*0.1)/(ni_small_atree*density_small_atree); //scaling  L^2/mu
        C1[i]  = 0.2/R1[i];
        ROOTONLY
          fprintf(pRC_File,"%2.16f  %2.16f \n",R1[i],C1[i]);
      }
      ROOTONLY
         fclose(pRC_File);
    }
    else{
      for (i = 0; i < Nout; ++i)
        fscanf(pRC_File," %lf %lf ",&R1[i],&C1[i]);
      fclose(pRC_File);
    }
    /* print summary */

    ROOTONLY{
      for (i = 0; i < Nout; ++i)
        fprintf(stdout,"PBC1D::setRC -- outlet %d: R1 = %f  C1 = %f \n",i,R1[i],C1[i]);
    }
}

void PBC1D::ResetRC(char *name){
    int i;
    FILE *pRC_File;
    char fname_RC[BUFSIZ];
    sprintf (fname_RC, "%s.RC", name );
    pRC_File = fopen(fname_RC,"r");
    if (pRC_File==NULL){
      ROOTONLY
         fprintf(stdout,"PBC1D::resetRC -- can not open RCfile \n");
    }
    else{
      for (i = 0; i < Nout; ++i)
        fscanf(pRC_File," %lf %lf ",&R1[i],&C1[i]);
      fclose(pRC_File);
      ROOTONLY{
        for (i = 0; i < Nout; ++i)
          fprintf(stdout,"PBC1D::resetRC -- outlet = %d R1 = %f  C1 = %f \n",i,R1[i],C1[i]);
      }
    }
}

void PBC1D::ReadFlowrateHistory(char *name){
 /* if history of flow rate exists - get it */
  FILE *pFile;
  char fname[BUFSIZ];

  /* open name.imp file for reading
     if file is not found return,
     flow rate history remains zero
     if file is present analyse it. */
  sprintf (fname, "%s.imp", name);
  pFile = fopen(fname,"r");
  if (pFile==NULL){
    ROOTONLY
      fprintf(stdout,
      "PBC1D::readFlowrateHistoryNew -- imp file was not found, setting flow rate history to 0.0 \n");
    return;
  }        
    
  /* analyse the data in name.imp file */
  int i,j;
  char buf[BUFSIZ];
  int Nsteps_modes; /* number of steps if space=='P' , number of modes if space=='M'*/ 
  int Noutlets_local;
  char space,Convolve_Fix;
  
  /* 1. check in if the flow rate history is saved in modal or physical space */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%c",&space);

  /* 2. read the number of outlets */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%d",&Noutlets_local);

  /* 3. read the number of steps if flow rate is given in physical space 
  or number of modes if it is given in modal space */
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%d",&Nsteps_modes);

  /* 4. read instruction on how to proceed, 
        Convolve_Fix =='F' compute pressure from given flow rate history
        and computed Impedance and fix it
        Convolve_Fix =='C' means keep updating flow rate every time step 
        and use convolution of flow rate and impedance to compute Pressure.*/
  fgets (buf, BUFSIZ, pFile);
  sscanf(buf,"%c",&Convolve_Fix);

 
  
  /* if space = modal we have enough information to proceed */
  
  if (Noutlets_local != Nout){ // inconsistent file
     ROOTONLY
       fprintf(stdout,
       "PBC1D::readFlowrateHistoryNew -- Number of outlets (%d) does not match number of outlet in *imp file (%d) \n",
       Nout,Noutlets_local);
     fclose(pFile);  
     return;  
  }
  
  char *p;
  rewind(pFile);	
  while (p = fgets (buf, BUFSIZ, pFile)){
    if (strstr (p, "FlowRateData")){
        break;
    }
  }
  if (space=='M'){	
    for (i = 0; i < Nout; i++){
    	fscanf(pFile,"%lf ",&A_cos_sin[i][0]);
    	for (j = 1; j <= Nsteps_modes; j++)
          fscanf(pFile,"%lf %lf",&A_cos_sin[i][2*j-1],&A_cos_sin[i][2*j]);
    }	
  }
  else{
  
     if (Nsteps_modes == Nsteps_per_cycle){
       for (i = 0; i < Nout; i++){
       	 for (j = 0; j < Nsteps_per_cycle; j++)	
    	    fscanf(pFile,"%lf", &FR_history[i][j]);
       }
     } 
     else{ // need to interpolate
       double dt = dparam("DT");  
       double dt_local = Tperiod/(Nsteps_modes-1.0);
       double tj, inv_dt_local_p2,coef_a,coef_b,coef_c;
       int J0,J1,J2;
       inv_dt_local_p2 = 1.0/(dt_local*dt_local);
       double *FR_local;
       FR_local = new double[Nsteps_modes];
       
       for (i = 0; i < Nout; i++){
       	
       	 for (j = 0; j < Nsteps_modes; j++)	
    	    fscanf(pFile,"%lf", &FR_local[j]);
       
       	 for (j = 0; j < Nsteps_per_cycle; j++){	
            tj = j*dt;
#if defined (__blrts__)            
            J1 = (int) round(tj/dt_local);
#else
            J1 = (int) (tj/dt_local);
#endif 
            if ( (J1 > 0) && J1 < (Nsteps_modes-1) ){
               J0 = J1-1;
               J2 = J1+1;
            }
            else{
                if (J1 <= 0){
                  J0 = 0;
                  J1 = 1;
                  J2 = 2;   
                }
                else{ 
                  J2 = Nsteps_modes-1;
                  J1 = J2-1;
                  J0 = J2-2;
                }
            }
            /*  compute interpolation coefficients use second order Lagrange */
            coef_a = (tj-dt_local*J1)*(tj-dt_local*J2)*(0.5*inv_dt_local_p2);
            coef_b = (tj-dt_local*J0)*(tj-dt_local*J2)*(-inv_dt_local_p2);
            coef_c = (tj-dt_local*J1)*(tj-dt_local*J0)*(0.5*inv_dt_local_p2);  
            FR_history[i][j] = FR_local[J0]*coef_a+FR_local[J1]*coef_b+FR_local[J2]*coef_c;
         } // end of for (j = 0; ...
       }    // end of for (i = 0; ...
       delete[] FR_local;
     }
   }

   rewind(pFile);
   j = 0;
   while (p = fgets (buf, BUFSIZ, pFile)){
     if (strstr (p, "FlowRateModal")){
       fscanf(pFile,"%d", &j); // total number of fourier modes per outlet        
       for (i = 0; i < Nout; i++){
         for (j = 0; j < (Nimpedance_modes*2+1); j++)
           fscanf(pFile,"%lf", &A_cos_sin[i][j]);
       }
       j = 1; //flag telling that data in modal space is also provided
       break;
     }
   }

   fclose(pFile);

     /* if data in modal space is provided - skip the next section */
    double *FR_history_temp;
    if (j == 1) 
      goto skip_Ftransform;
    
    /* if time > Tcycle do fourier transform of flow-rate history, obtain fourier coef.  */
    /* if Flow-rate is given in Physical space and  time > Tcycle 
        obtain fourier coefficient */
    /* if Flow-rate is given in modal space transform it into physical space.  */
    FR_history_temp = dvector(0,Nsteps_per_cycle-2);    

    if ((dparam("t") > Tperiod) && (space=='P')){
      for (i = 0; i < Nout; i++){
        memcpy(FR_history_temp,FR_history[i],(Nsteps_per_cycle-1)*sizeof(double));
        FilterDFT(Nsteps_per_cycle-1, Nimpedance_modes, FR_history_temp,A_cos_sin[i],0);
      }
    }
    free(FR_history_temp);

    skip_Ftransform:

    if ( space=='M'){
      double dt = dparam("DT");
      double arg;
      int k;
      
      for (i = 0; i < Nout; i++){
        for (j = 0; j < Nsteps_per_cycle; j++){
          arg = omega_small_atree*j*dt;
          FR_history[i][j] = A_cos_sin[i][0];
          for (k=1; k <= Nsteps_modes; k++)
          FR_history[i][j] += A_cos_sin[i][k*2-1]*cos(arg*k) + A_cos_sin[i][k*2]*sin(arg*k);
        }
      }
    }
}

void PBC1D::SaveFlowrateHistory(char *name){

  int i,j;
  FILE *pFile;
  char fname[BUFSIZ];

  sprintf (fname, "%s.imp", name);
  pFile = fopen(fname,"w");

  fprintf(pFile," P /* P - physical, M - modal space */ \n");
  fprintf(pFile,"%d /* number of outlets */ \n",Nout);
  fprintf(pFile,"%d /* Number of steps   */  \n",Nsteps_per_cycle);  
  fprintf(pFile,"C  /* instruction: C - P(t) = conv. (F,Z) , F - fix P(t) */ \n");
  fprintf(pFile,"FlowRateData \n");
  for (i = 0; i < Nout; i++){
    for (j = 0; j < Nsteps_per_cycle; j++)
      fprintf(pFile," %.10f \n",FR_history[i][j]);
  }
  fprintf(pFile,"FlowRateModal \n");
  fprintf(pFile,"%d  /* number of fourier modes */ \n",Nimpedance_modes*2+1);
  for (i = 0; i < Nout; i++){
    for (j = 0; j < (Nimpedance_modes*2+1); j++)
      fprintf(pFile," %.10f \n",A_cos_sin[i][j]);
  }
  fclose(pFile);
}


double PBC1D::GetPval(Bndry *Pbc){

   int i;
   for (i = 0; i < Nout; ++i){
     if (strcmp(Pbc->blabel,standard_labels[i]) == 0)
        return Pressure[i];
   }  
   /* if no match for standard_labels was found return zero */
   return 0.0;
}


void PBC1D::UpdateTimestepCycle(){
   double t = dparam("t");
   double dt = dparam("DELT");
   /* parameter "t" = time at the end of the cuurent time step,
      that is t^(n+1), thus we need to substruct "dt" from "t" to get 
      the time_step_in_cycle correctly */

   t -= dt; 

   t = t-((int)((t+dt*0.1)/Tperiod))*Tperiod;
     time_step_in_cycle = int ((t+0.1*dt)/dt);
}

void PBC1D::UpdateFRHistory(double *flowrate){

  /* INPUT: values of flowd-rate at outlets at current time 
     function updates arrays FR_history where flow-rate history 
     over a cycle is stored 
     if (t == Tperiod-dt) filter FR_history using Fourier transform 
     cut all frecuencies higher then  Nimpedance_modes 
     store values of Fourier coefficients.
   
     starting from the end of second cycle average values of Fourier
     coefficients with those from previous cycle
  
*/

  int i,FLAG_filter;
  for (i = 0; i < Nout; ++i)
    FR_history[i][time_step_in_cycle] = flowrate[i];


   /* predict flow-rate for the rest of the first time period */
  if   (dparam("t") < (Tperiod-dparam("DELT"))){
    double factor;
    int j;

    for (i = 0; i < Nout; ++i){
      factor = FR_history[i][time_step_in_cycle];
      for (j = time_step_in_cycle; j < Nsteps_per_cycle; j++)
        FR_history[i][j] = factor*cos(M_PI*0.5*(j-time_step_in_cycle)/(Nsteps_per_cycle-time_step_in_cycle));
    }
  }

  if (dparam("t") < (1.5*Tperiod))
     FLAG_filter = 0;
  else
     FLAG_filter = 1;
 
  if  (time_step_in_cycle == (Nsteps_per_cycle-2) ) {
    fprintf(stdout,"PBC1D::update_FR_history, FLAG_filter = %d time_step_in_cycle = %d Nsteps_per_cycle = %d \n",FLAG_filter,time_step_in_cycle,Nsteps_per_cycle);
    for (i = 0; i < Nout; ++i)
      FilterDFT(Nsteps_per_cycle-1, Nimpedance_modes, FR_history[i],A_cos_sin[i],FLAG_filter);
  }
}

void PBC1D::ComputePressureSteady(double *flowrate){

  double Po=0.0,alpha;
  double inv_dt = 1.0/dparam("DELT");
  int i;
  double DPSCAL = dparam("DPSCAL");

  for (i = 0; i < Nout; ++i){
     alpha = R1[i]*C1[i]*inv_dt;
     Pressure[i] = DPSCAL*(R1[i]*flowrate[i]+Po+Pressure[i]*alpha)/(1.0+alpha);
     flowrate_RCR_old[i] = flowrate[i];
  }
}

void PBC1D::ComputePressureRCNonsteady(double *flowrate){

  double Po=0.0,alpha;
  double inv_dt = 1.0/dparam("DELT");
  int i;
  double DPSCAL = dparam("DPSCAL");
  static int INIT_FLAG=0;
  static double **Asin_coef, **Acos_coef;
  static double *Ao;
  int Number_rc_modes = 25;
  if (INIT_FLAG == 0){

     Ao = dvector(0,Nout-1);
     Asin_coef = dmatrix(0,Nout-1,0,Number_rc_modes-1);
     Acos_coef = dmatrix(0,Nout-1,0,Number_rc_modes-1);
     memset(Ao,'\0',Nout*sizeof(double));

     for (i = 0; i < Nout; ++i){
       memset(Asin_coef[i],'\0',Number_rc_modes*sizeof(double));
       memset(Acos_coef[i],'\0',Number_rc_modes*sizeof(double));
     }

/*   ICA peak follows ECA peak */
/*
     Ao[1]           = 1.531376625418332;
     Asin_coef[1][0] = -0.1989430208905198;
     Asin_coef[1][1] = -0.0742076020935493;
     Asin_coef[1][2] = -0.1319826486529148;
     Asin_coef[1][3] = 0.0419093159711077;
     Asin_coef[1][4] = 0.3877176977820486;
     Asin_coef[1][5] = 0.0232159698060311;
     Asin_coef[1][6] = -0.0262147642389971;
     Asin_coef[1][7] = 0.2045538297607134;
     Asin_coef[1][8] = -0.0075767679009015;
     Asin_coef[1][9] = -0.0093384094198911;
     Asin_coef[1][10] = 0.0673670044939792;
     Asin_coef[1][11] = -0.0081829522598334;
     Asin_coef[1][12] = 0.0613477553174892;
     Asin_coef[1][13] = -0.0303619733810808;
     Asin_coef[1][14] = 0.0142090613526252;
     Asin_coef[1][15] = 0.0425480821357585;
     Asin_coef[1][16] = 0.0135161983852775;
     Asin_coef[1][17] = -0.0028540616978910;
     Asin_coef[1][18] = -0.0063017352548369;
     Asin_coef[1][19] = 0.0305215998651123;
     Asin_coef[1][20] = 0.0088068128602008;
     Asin_coef[1][21] = -0.0101996601574633;
     Asin_coef[1][22] = -0.0125411052135368;
     Asin_coef[1][23] = -0.0033854748144212;
     Asin_coef[1][24] = 0.0061052164312761;
     Asin_coef[1][25] = -0.0082433018261473;
     Asin_coef[1][26] = -0.0125181216549941;
     Asin_coef[1][27] = -0.0023898113172775;

     Acos_coef[1][0] = -0.0602519488910583;
     Acos_coef[1][1] = -0.3278505713662204;
     Acos_coef[1][2] = 0.5841706200860648;
     Acos_coef[1][3] = 0.1104895873535616;
     Acos_coef[1][4] = 0.0805509153630266;
     Acos_coef[1][5] = -0.0092731071733014;
     Acos_coef[1][6] = -0.0363197215390623;
     Acos_coef[1][7] = 0.0027833401869274;
     Acos_coef[1][8] = -0.0338576111883325;
     Acos_coef[1][9] = -0.0379486162982755;
     Acos_coef[1][10] = -0.0046348340050727;
     Acos_coef[1][11] = 0.0153143714282117;
     Acos_coef[1][12] = -0.0651553715602252;
     Acos_coef[1][13] = -0.0172588509471056;
     Acos_coef[1][14] = 0.0368305682492296;
     Acos_coef[1][15] = -0.0395710730063579;
     Acos_coef[1][16] = -0.0175947796234016;
     Acos_coef[1][17] = -0.0298682561384160;
     Acos_coef[1][18] = -0.0021733679416264;
     Acos_coef[1][19] = -0.0116731854913752;
     Acos_coef[1][20] = -0.0345565452677156;
     Acos_coef[1][21] = -0.0185414481637263;
     Acos_coef[1][22] = -0.0062416534810577;
     Acos_coef[1][23] = 0.0026100811234523;
     Acos_coef[1][24] = -0.0063449196505491;
     Acos_coef[1][25] = -0.0138819538597907;
     Acos_coef[1][26] = 0.0030197962723845;
     Acos_coef[1][27] = 0.0063134900380622;
*/

/*   ICA peak preceeds ECA peak */

     Ao[1] = 1.491129133777125;
     Asin_coef[1][0] = -0.2307293645462276;
     Asin_coef[1][1] = 0.0673947478743045;
     Asin_coef[1][2] = 0.1234666688951011;
     Asin_coef[1][3] = 0.1207507745138338;
     Asin_coef[1][4] = 0.1766623174274239;
     Asin_coef[1][5] = -0.2649975418627119;
     Asin_coef[1][6] = -0.1840905783438354;
     Asin_coef[1][7] = 0.0900406071895171;
     Asin_coef[1][8] = -0.1442740737984572;
     Asin_coef[1][9] = -0.0718543502039253;
     Asin_coef[1][10] = 0.0328518197155938;
     Asin_coef[1][11] = -0.1141498607343912;
     Asin_coef[1][12] = -0.0574531034700849;
     Asin_coef[1][13] = -0.0918472012510603;
     Asin_coef[1][14] = -0.0219902119451843;
     Asin_coef[1][15] = 0.0112862075060890;
     Asin_coef[1][16] = 0.0324544483650855;
     Asin_coef[1][17] = 0.0594224545040913;
     Asin_coef[1][18] = 0.0330739610841887;
     Asin_coef[1][19] = 0.0244314293513142;
     Asin_coef[1][20] = -0.0081478802975255;
     Asin_coef[1][21] = -0.0101407819828567;
     Asin_coef[1][22] = -0.0087021462221979;
     Asin_coef[1][23] = -0.0099920432306531;
     Asin_coef[1][24] = -0.0015217006827334;

     Acos_coef[1][0] = -0.0337974595816419;
     Acos_coef[1][1] = -0.2160202643735548;
     Acos_coef[1][2] = 0.4670312618426233;
     Acos_coef[1][3] = -0.2925355336526577;
     Acos_coef[1][4] = -0.2510428244255957;
     Acos_coef[1][5] = -0.0697772519851484;
     Acos_coef[1][6] = -0.0050110370638509;
     Acos_coef[1][7] = -0.0061134590871459;
     Acos_coef[1][8] = -0.0007178761923128;
     Acos_coef[1][9] = 0.0232338720924626;
     Acos_coef[1][10] = -0.0086846036211113;
     Acos_coef[1][11] = 0.0013297313216929;
     Acos_coef[1][12] = -0.0077327261415910;
     Acos_coef[1][13] = 0.0733270436054532;
     Acos_coef[1][14] = 0.1093802405979988;
     Acos_coef[1][15] = 0.0452343906395702;
     Acos_coef[1][16] = 0.0794167402990195;
     Acos_coef[1][17] = 0.0150675507207043;
     Acos_coef[1][18] = -0.0142026525177695;
     Acos_coef[1][19] = -0.0260577734710864;
     Acos_coef[1][20] = -0.0224402024063286;
     Acos_coef[1][21] = 0.0005784871650041;
     Acos_coef[1][22] = 0.0000962867558026;
     Acos_coef[1][23] = 0.0065276078059094;
     Acos_coef[1][24] = 0.0055496383708904;
  
     INIT_FLAG = 1;
  }

     /* fix R1 compute R2,  R2/R1 ratio is given in fourier space, then R2 = R1 * ratio */
  for (i = 1; i < Nout; ++i)
    R1[i] = R1[0]*Ao[i];

  double t = dparam("t"),sin_alpha,cos_alpha;
  double omega_t = 2.0*3.14159265358979*t/Tperiod;
  int j;
  
  for (i = 0; i < Number_rc_modes; ++i){
    alpha = omega_t*(i+1);
    sin_alpha = sin(alpha);
    cos_alpha = cos(alpha);
    for (j = 1; j < Nout; ++j)
      R1[j] += R1[0]*(Asin_coef[j][i]*sin_alpha+Acos_coef[j][i]*cos_alpha);
  }

  for (i = 0; i < Nout; ++i){
     alpha = R1[i]*C1[i]*inv_dt;
     Pressure[i] = DPSCAL*(R1[i]*flowrate[i]+Po+Pressure[i]*alpha)/(1.0+alpha);
     flowrate_RCR_old[i] = flowrate[i];
  }
}

void PBC1D::ComputePressure(){

/* compute pressure using convolution 
   presure is computed in [cm,sec,gr], so scale it to nondim units */



   int i,time_step_in_cycle_advanced;
   double DPSCAL = dparam("DPSCAL");

   if (time_step_in_cycle == (Nsteps_per_cycle-2)) 
       time_step_in_cycle_advanced = 0;
   else
       time_step_in_cycle_advanced = time_step_in_cycle+1;


   for (i = 0; i < Nout; ++i){
      if (!parallel_convolution)
         ParallelConvolve(Nsteps_per_cycle-1, FR_history[i], impedance[i], &Pressure[i], time_step_in_cycle_advanced);
      else{
       ParallelConvolve(Nsteps_per_cycle-1, FR_history[i], impedance[i], &Pressure[i], Imp_index_start, (Imp_index_start+Imp_length_local), time_step_in_cycle_advanced);

       DO_PARALLEL{
         double temp = 0.0;
         gdsum(&Pressure[i],1,&temp);
       }
     }
     Pressure[i] *= DPSCAL*0.1/(Tperiod*density_small_atree); 
  }

}

#endif
