/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#include <stdio.h>
#define MY_PREFIX fix_shake_cuda
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"
#include "fix_shake_cuda_cu.h"
#include "cuda_pair_virial_kernel_nc.cu"

#define _shake_atom           MY_AP(shake_atom)
#define _shake_type           MY_AP(shake_type)
#define _shake_flag           MY_AP(shake_flag)
#define _xshake               MY_AP(xshake)
#define _dtfsq                MY_AP(dtfsq)
#define _bond_distance        MY_AP(bond_distance)
#define _angle_distance       MY_AP(angle_distance)
#define _max_iter			  MY_AP(max_iter)
#define _tolerance			  MY_AP(tolerance)
__device__ __constant__ int* _shake_atom;
__device__ __constant__ int* _shake_type;
__device__ __constant__ int* _shake_flag;
__device__ __constant__ X_FLOAT3* _xshake;
__device__ __constant__ F_FLOAT _dtfsq;
__device__ __constant__ X_FLOAT* _bond_distance;
__device__ __constant__ X_FLOAT* _angle_distance;
__device__ __constant__ int _max_iter;
__device__ __constant__ X_FLOAT _tolerance;

#include "fix_shake_cuda_kernel.cu"

void Cuda_FixShakeCuda_UpdateNmax(cuda_shared_data* sdata)
{
		cudaMemcpyToSymbol(MY_CONST(x)       , & sdata->atom.x    .dev_data, sizeof(X_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(v)       , & sdata->atom.v    .dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(f)       , & sdata->atom.f    .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(mask)    , & sdata->atom.mask .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(tag)     , & sdata->atom.tag  .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal        , sizeof(int)      );
		cudaMemcpyToSymbol(MY_CONST(nmax)    , & sdata->atom.nmax          , sizeof(int)      );
		cudaMemcpyToSymbol(MY_CONST(rmass)   , & sdata->atom.rmass.dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(type)    , & sdata->atom.type .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(map_array), & sdata->atom.map_array .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(vatom)   , & sdata->atom.vatom.dev_data, sizeof(ENERGY_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(debugdata), & sdata->debugdata         , sizeof(int*)     );
}

void Cuda_FixShakeCuda_UpdateDomain(cuda_shared_data* sdata)
{
	cudaMemcpyToSymbol(MY_CONST(periodicity), sdata->domain.periodicity		, sizeof(int)*3    );
	cudaMemcpyToSymbol(MY_CONST(prd)		, sdata->domain.prd				, sizeof(X_FLOAT)*3    );
	cudaMemcpyToSymbol(MY_CONST(triclinic)  , &sdata->domain.triclinic		, sizeof(int)      );
	cudaMemcpyToSymbol(MY_CONST(h)			, sdata->domain.h				, sizeof(X_FLOAT)*6    );
}

void Cuda_FixShakeCuda_UpdateBuffer(cuda_shared_data* sdata,int size)
{
		if(sdata->buffersize<size)
		{
			MYDBG(printf("Cuda_FixShakeCuda Resizing Buffer at %p with %i kB to\n",sdata->buffer,sdata->buffersize);)
			CudaWrapper_FreeCudaData(sdata->buffer,sdata->buffersize);
			sdata->buffer = CudaWrapper_AllocCudaData(size);
			sdata->buffersize=size;
			sdata->buffer_new++;
			MYDBG(printf("New buffer at %p with %i kB\n",sdata->buffer,sdata->buffersize);)
			
		}
		cudaMemcpyToSymbol(MY_CONST(buffer) , & sdata->buffer, sizeof(int*)     );
}

void Cuda_FixShakeCuda_Init(cuda_shared_data* sdata,X_FLOAT dtv, F_FLOAT dtfsq,
						void* shake_flag,void* shake_atom,void* shake_type, void* xshake,
						void* bond_distance,void* angle_distance,void* virial,
						int max_iter,X_FLOAT tolerance)
{
	Cuda_FixShakeCuda_UpdateNmax(sdata);
	Cuda_FixShakeCuda_UpdateDomain(sdata);
	cudaMemcpyToSymbol(MY_CONST(shake_atom)        , & shake_atom 	  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(shake_type)        , & shake_type 	  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(shake_flag)        , & shake_flag 	  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(xshake)            , & xshake     	  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(dtv)               , & dtv        	  , sizeof(X_FLOAT));
	cudaMemcpyToSymbol(MY_CONST(dtfsq)             , & dtfsq      	  , sizeof(F_FLOAT));
	cudaMemcpyToSymbol(MY_CONST(bond_distance)     , & bond_distance  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(angle_distance)    , & angle_distance , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(virial)     	   , & virial  		  , sizeof(void*) );
	cudaMemcpyToSymbol(MY_CONST(flag)  			   , &sdata->flag	  , sizeof(int*)  );
	cudaMemcpyToSymbol(MY_CONST(max_iter)  		   , &max_iter  	  , sizeof(int)   );
	cudaMemcpyToSymbol(MY_CONST(tolerance)  	   , &tolerance  	  , sizeof(X_FLOAT));
	
	if(sdata->atom.mass_host)
	cudaMemcpyToSymbol(MY_CONST(mass),& sdata->atom.mass.dev_data , sizeof(V_FLOAT*) );
	cudaMemcpyToSymbol(MY_CONST(rmass_flag), & sdata->atom.rmass_flag       , sizeof(int)      ); //

    cudaMemcpyToSymbol(MY_CONST(flag)  , &sdata->flag, sizeof(int*));
	
}

void Cuda_FixShakeCuda_UnconstrainedUpdate(cuda_shared_data* sdata)
{
	if(sdata->atom.update_nmax) 
		Cuda_FixShakeCuda_UpdateNmax(sdata);
	if(sdata->atom.update_nlocal) 		
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal , sizeof(int)); 
	if(sdata->buffer_new)
		Cuda_FixShakeCuda_UpdateBuffer(sdata,10*sizeof(double));
	int3 layout=getgrid(sdata->atom.nlocal);
	dim3 threads(layout.z, 1, 1);
	dim3 grid(layout.x, layout.y, 1);
	
	FixShakeCuda_UnconstrainedUpdate_Kernel<<<grid, threads>>> ();
	cudaThreadSynchronize();
	CUT_CHECK_ERROR("FixShakeCuda_UnconstrainedUpdate: Kernel execution failed");	
}

void Cuda_FixShakeCuda_Shake(cuda_shared_data* sdata,int vflag,int vflag_atom,int* list,int nlist)
{
	if(sdata->atom.update_nmax) 
		Cuda_FixShakeCuda_UpdateNmax(sdata);
	if(sdata->domain.update) 
		Cuda_FixShakeCuda_UpdateDomain(sdata);
	if(sdata->atom.update_nlocal) 		
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal , sizeof(int)); 
	int3 layout=getgrid(sdata->atom.nlocal,6*sizeof(ENERGY_FLOAT),64);
	dim3 threads(layout.z, 1, 1);
	dim3 grid(layout.x, layout.y, 1);
	if(sdata->buffer_new)
		Cuda_FixShakeCuda_UpdateBuffer(sdata,grid.x*grid.y*6*sizeof(ENERGY_FLOAT));
	
	BindXTypeTexture(sdata);
	
	FixShakeCuda_Shake_Kernel<<<grid, threads,6*threads.x*sizeof(ENERGY_FLOAT)>>> (vflag,vflag_atom,list,nlist);
	cudaThreadSynchronize();
	
	CUT_CHECK_ERROR("FixShakeCuda_Shake: Kernel execution failed");	

    if(vflag) 
	{
		int n=grid.x*grid.y;
		grid.x=6;
		grid.y=1;
		threads.x=256;
		MY_AP(PairVirialCompute_reduce)<<<grid,threads,threads.x*sizeof(ENERGY_FLOAT)>>>(n);
		cudaThreadSynchronize();
		CUT_CHECK_ERROR("Cuda_FixShakeCuda: (no binning) virial compute Kernel execution failed");
	}
	
}

int Cuda_FixShakeCuda_PackComm(cuda_shared_data* sdata,int n,int iswap,void* buf_send,int* pbc,int pbc_flag)
{
	if(sdata->atom.update_nmax) 
		Cuda_FixShakeCuda_UpdateNmax(sdata);
	if(sdata->atom.update_nlocal) 		
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal        , sizeof(int)      );
	int size=n*3*sizeof(X_FLOAT);
	if(sdata->buffer_new or (size>sdata->buffersize))
		Cuda_FixShakeCuda_UpdateBuffer(sdata,size);

	X_FLOAT dx=0.0;
	X_FLOAT dy=0.0;
	X_FLOAT dz=0.0;
 	if (pbc_flag != 0) {
    if (sdata->domain.triclinic == 0) {
      dx = pbc[0]*sdata->domain.prd[0];
      dy = pbc[1]*sdata->domain.prd[1];
      dz = pbc[2]*sdata->domain.prd[2];
    } else {
      dx = pbc[0]*sdata->domain.prd[0] + pbc[5]*sdata->domain.xy + pbc[4]*sdata->domain.xz;
      dy = pbc[1]*sdata->domain.prd[1] + pbc[3]*sdata->domain.yz;
      dz = pbc[2]*sdata->domain.prd[2];
    }}	

	int3 layout=getgrid(n);
	dim3 threads(layout.z, 1, 1);
	dim3 grid(layout.x, layout.y, 1);
	  
	if(sdata->atom.nlocal>0)
	{
	  cudaMemset( sdata->flag,0,sizeof(int));
	  FixShakeCuda_PackComm_Kernel<<<grid, threads,0>>>((int*) sdata->comm.sendlist.dev_data,n,sdata->comm.maxlistlength,iswap,dx,dy,dz);
	  cudaThreadSynchronize();
	  cudaMemcpy(buf_send, sdata->buffer, n*3*sizeof(X_FLOAT), cudaMemcpyDeviceToHost);
	  int aflag;
	  cudaMemcpy(&aflag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);
	  if(aflag!=0) printf("aflag PackComm: %i\n",aflag);
	  CUT_CHECK_ERROR("Cuda_FixShakeCuda_PackComm: Kernel execution failed");
		
	}		
    return 3*n;
}

int Cuda_FixShakeCuda_PackComm_Self(cuda_shared_data* sdata,int n,int iswap,int first,int* pbc,int pbc_flag)
{
	if(sdata->atom.update_nmax) 
		Cuda_FixShakeCuda_UpdateNmax(sdata);
	if(sdata->atom.update_nlocal) 		
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal        , sizeof(int)      );
	int size=n*3*sizeof(X_FLOAT);
	if(sdata->buffer_new or (size>sdata->buffersize))
		Cuda_FixShakeCuda_UpdateBuffer(sdata,size);
	static int count=-1;
	count++;
	X_FLOAT dx=0.0;
	X_FLOAT dy=0.0;
	X_FLOAT dz=0.0;
 	if (pbc_flag != 0) {
    if (sdata->domain.triclinic == 0) {
      dx = pbc[0]*sdata->domain.prd[0];
      dy = pbc[1]*sdata->domain.prd[1];
      dz = pbc[2]*sdata->domain.prd[2];
    } else {
      dx = pbc[0]*sdata->domain.prd[0] + pbc[5]*sdata->domain.xy + pbc[4]*sdata->domain.xz;
      dy = pbc[1]*sdata->domain.prd[1] + pbc[3]*sdata->domain.yz;
      dz = pbc[2]*sdata->domain.prd[2];
    }}	



	int3 layout=getgrid(n);
	dim3 threads(layout.z, 1, 1);
	dim3 grid(layout.x, layout.y, 1);
	if(sdata->atom.nlocal>0)
	{
	  FixShakeCuda_PackComm_Self_Kernel<<<grid, threads,0>>>((int*) sdata->comm.sendlist.dev_data,n,sdata->comm.maxlistlength,iswap,dx,dy,dz,first);
	  cudaThreadSynchronize();
	  CUT_CHECK_ERROR("Cuda_CommCuda_PackComm_Self: Kernel execution failed");
	}	
	
    return 3*n;
}

void Cuda_FixShakeCuda_UnpackComm(cuda_shared_data* sdata,int n,int first,void* buf_recv)
{
	if(sdata->atom.update_nmax) 
		Cuda_FixShakeCuda_UpdateNmax(sdata);
	if(sdata->atom.update_nlocal) 		
		cudaMemcpyToSymbol(MY_CONST(nlocal)  , & sdata->atom.nlocal        , sizeof(int)      );
	int size=n*3*sizeof(X_FLOAT);
	if(sdata->buffer_new or (size>sdata->buffersize))
		Cuda_FixShakeCuda_UpdateBuffer(sdata,size);

	int3 layout=getgrid(n);
	dim3 threads(layout.z, 1, 1);
	dim3 grid(layout.x, layout.y, 1);
	if(sdata->atom.nlocal>0)
	{
	  cudaMemcpy(sdata->buffer,(void*)buf_recv, n*3*sizeof(X_FLOAT), cudaMemcpyHostToDevice);
	  FixShakeCuda_UnpackComm_Kernel<<<grid, threads,0>>>(n,first);
	  cudaThreadSynchronize();
	  CUT_CHECK_ERROR("Cuda_FixShakeCuda_UnpackComm: Kernel execution failed");
		
	}		
}
