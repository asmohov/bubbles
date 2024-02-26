//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rt.cpp
//! \brief Problem generator for RT instabilty.
//!
//! Note the gravitational acceleration is hardwired to be 0.1. Density difference is
//! hardwired to be 2.0 in 2D, and is set by the input parameter `problem/rhoh` in 3D
//! (default value is 3.0). This reproduces 2D results of Liska & Wendroff, 3D results of
//! Dimonte et al.
//!
//! FOR 2D HYDRO:
//! Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to match Liska
//! & Wendroff. Interface is at y=0; perturbation added to Vy. Gravity acts in y-dirn.
//! Special reflecting boundary conditions added in x2 to improve hydrostatic eqm
//! (prevents launching of weak waves) Atwood number A=(d2-d1)/(d2+d1)=1/3. Options:
//!    - iprob = 1  -- Perturb V2 using single mode
//!    - iprob != 1 -- Perturb V2 using multiple mode
//!
//! FOR 3D:
//! Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1, gamma=5/3 to
//! match Dimonte et al.  Interface is at z=0; perturbation added to Vz. Gravity acts in
//! z-dirn. Special reflecting boundary conditions added in x3.  A=1/2.  Options:
//!    - iprob = 1 -- Perturb V3 using single mode
//!    - iprob = 2 -- Perturb V3 using multiple mode
//!    - iprob = 3 -- B rotated by "angle" at interface, multimode perturbation
//!
//! REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// made global to share with BC functions
Real grav_acc;
} // namespace

int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);
  if (mesh_size.nx3 == 1) {  // 2D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, ProjectPressureInnerX2);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, ProjectPressureOuterX2);
  } else { // 3D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, ProjectPressureInnerX3);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, ProjectPressureOuterX3);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Rayleigh-Taylor instability test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::int64_t iseed = -1;
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  Real kx = 2.0*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 2.0*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real kz = 2.0*(PI)/(pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);
  

  std::cout<<"mesh_size_x1max-x1min"<<pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min<<std::endl;
  std::cout<<"kx"<<kx<<std::endl;

  // Read perturbation amplitude, problem switch, density ratio
  Real amp = pin->GetReal("problem","amp");
  int iprob = pin->GetInteger("problem","iprob");
  Real drat = pin->GetOrAddReal("problem","drat",3.0);
  // Read spheromak parameters
  //Real rad_sph = pin->GetReal("problem","rad_sph"); deprecated in favor of alpha paramterization
  Real x0 = pin->GetReal("problem","x0_sph");
  Real y0 = pin->GetReal("problem","y0_sph");
  Real z0 = pin->GetReal("problem","z0_sph");
  Real orient_sph = pin->GetReal("problem","orient_sph");
  Real Alpha = pin->GetReal("problem","alpha_sph");
  Real beta_out = pin->GetReal("problem","beta_out");
  Real beta_in = pin->GetReal("problem","beta_in");
  //set b0 from density at top of box and beta_out
  Real azero= pcoord->x3v(ke);//set exponential coefficient to z length of box
  Real rho_top = azero*std::exp(grav_acc*azero); //define rho_max as the density at the top of the box
  Real b0 = std::sqrt(2*(rho_top/beta_out));
  // 2D PROBLEM ---------------------------------------------------------------

  if (block_size.nx3 == 1) {
    grav_acc = phydro->hsrc.GetG2();
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real den=1.0;
          if (pcoord->x2v(j) > 0.0) den *= drat;

          if (iprob == 1) {
            phydro->u(IM2,k,j,i) = (1.0 + std::cos(kx*pcoord->x1v(i)))*
                                   (1.0 + std::cos(ky*pcoord->x2v(j)))/4.0;
          } else {
            phydro->u(IM2,k,j,i) = (ran2(&iseed) - 0.5)*(1.0+std::cos(ky*pcoord->x2v(j)));
          }

          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) *= (den*amp);
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pcoord->x2v(j)))/gm1;
            phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/den;
          }
        }
      }
    }

    // initialize interface B, same for all iprob
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      Real b0 = pin->GetReal("problem","b0");
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = b0;
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }

    // 3D PROBLEM ----------------------------------------------------------------

  } else {
    //initialize density (cell centered)
    //moved to below b field, leaving original here for safety
    //grav_acc = phydro->hsrc.GetG3();
    //for (int k=ks; k<=ke; k++) {
      //for (int j=js; j<=je; j++) {
        //for (int i=is; i<=ie; i++) {
          //Real azero= pcoord->x3v(ke);//set exponential coefficient to z length of box
          //define x y and z
	  //Real x1 = pcoord->x1v(i);
      	  //Real y1 = pcoord->x2v(j);
      	  //Real z1 = pcoord->x3v(k);
	  //Real den=azero*std::exp(grav_acc*z1);
	  //std::cout<< "attempting to set exponential density: " << den <<std::endl;
          //if (pcoord->x3v(k) > 0.0) den *= drat;

          //if (iprob == 1) {
             //phydro->u(IM3,k,j,i) = (1.0 + std::cos(kx*(pcoord->x1v(i))))/8.0
               //                    *(1.0 + std::cos(ky*pcoord->x2v(j)))
               //                    *(1.0 + std::cos(kz*pcoord->x3v(k)));
          //} else {
            //phydro->u(IM3,k,j,i) = amp*(ran2(&iseed) - 0.5)*(
                //1.0 + std::cos(kz*pcoord->x3v(k)));
          //}

          //phydro->u(IDN,k,j,i) = den;
          //phydro->u(IM1,k,j,i) = 0.0;
          //phydro->u(IM2,k,j,i) = 0.0;
          //phydro->u(IM3,k,j,i) *= (den*amp);
          //if (NON_BAROTROPIC_EOS) {
          //  phydro->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pcoord->x3v(k)))/gm1;
          //  phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/den;
          //}
        //}
      //}
    //}

    // initialize interface B
    if (MAGNETIC_FIELDS_ENABLED) {
      // Read magnetic field strength, angle [in degrees, 0 is along +ve X-axis]
      //determine b0 from plasma beta and box size, moved to global
      //Real azero= pcoord->x3v(ke);//set exponential coefficient to z length of box 
      //Real rho_top = azero*std::exp(grav_acc*azero); //define rho_max as the density at the top of the box
      //Real b0 = std::sqrt(2*(rho_top/beta_out)); 
      //Real b0 = pin->GetReal("problem","b0");
      Real angle = pin->GetReal("problem","angle");
      angle = (angle/180.)*PI;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
		//Bx
	    Real x1 = pcoord->x1f(i);
            Real y1 = pcoord->x2f(j);
            Real z1 = pcoord->x3f(k);
	    //adjust relative to spheromak center
	    Real x=x1-x0;
	    Real y=y1-y0;
	    Real z=z1-z0+1e-15;//catch 0 case
	    Real r0 = 4.493/Alpha; //define r0 in terms of Alpha
	    //check if inside our outside spheromak
	    Real r=std::sqrt(x*x+y*y+z*z);
	    if (r<r0){//internal Bx
		//std::cout<<"inside spheromak radius"<<std::endl;
                Real Bx_int=(1/(Alpha*r*r*r*r*r))*(z*Alpha*(-3*x*r+Alpha*x*x*y*r/(std::sqrt(z*z))+y*Alpha*(y*y+z*z)*r/(std::sqrt(z*z)))*std::cos(r*Alpha) -
                           z*(y*r*Alpha*r/(std::sqrt(z*z))+x*x*x*Alpha*Alpha+x*(-3+y*y*Alpha*Alpha+Alpha*Alpha*z*z))*std::sin(r*Alpha));
		Bx_int = b0*Bx_int;
		//assign according to orientation of spheromak
		if(orient_sph==0){pfield->b.x1f(k,j,i)=Bx_int;}
                if(orient_sph==1){pfield->b.x3f(k,j,i)=Bx_int;}//swap Bx and Bz for this case
		if(orient_sph==2){pfield->b.x1f(k,j,i)=Bx_int;}
		//std::cout<<"inside spheromak radius Bx is"<<Bx_int<<std::endl;
			}
	    else{//external Bx
		Real Bx_ext = -3*r0*r0*r0*x*z/(2*r*r*r*r*r);
		Bx_ext = b0*Bx_ext;
		//assign according to orientation of spheromak
                if(orient_sph==0){pfield->b.x1f(k,j,i)=Bx_ext;}
                if(orient_sph==1){pfield->b.x3f(k,j,i)=Bx_ext;}//swap Bx and Bz for this case
                if(orient_sph==2){pfield->b.x1f(k,j,i)=Bx_ext;}
		//placeholder
		//if(pcoord->x3v(k)>0.0){pfield->b.x1f(k,j,i)=0.0;}
	        //else{pfield->b.x1f(k,j,i)=0.0;}	
		}
          }
        }
      }
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
		  //By
            Real x1 = pcoord->x1f(i);
            Real y1 = pcoord->x2f(j);
            Real z1 = pcoord->x3f(k);
            //adjust relative to spheromak center
            Real x=x1-x0;
            Real y=y1-y0;
            Real z=z1-z0+1e-15;
	    Real r0 = 4.493/Alpha; //define r0 in terms of Alpha
            //check if inside our outside spheromak
            Real r=std::sqrt(x*x+y*y+z*z);
            if (r<r0){//internal By
                Real By_int=(1/(Alpha*r*r*r*r*r))*(z*Alpha*(-3*y*r+Alpha*y*y*x*r/(std::sqrt(z*z))+x*Alpha*(x*x+z*z)*r/(std::sqrt(z*z)))*std::cos(r*Alpha) -
                           z*(x*r*Alpha*r/(std::sqrt(z*z))+y*y*y*Alpha*Alpha+y*(-3+x*x*Alpha*Alpha+Alpha*Alpha*z*z))*std::sin(r*Alpha));
		By_int=b0*By_int;
		//assign according to orientation of spheromak
                if(orient_sph==0){pfield->b.x2f(k,j,i)=By_int;}
                if(orient_sph==1){pfield->b.x2f(k,j,i)=By_int;}
                if(orient_sph==2){pfield->b.x3f(k,j,i)=By_int;}//swap By and Bz in this case
		//std::cout<<"inside spheromak radius By is"<<By_int<<std::endl;
                        }
            else{//external By
                 //placeholder
		 Real By_ext = -3*r0*r0*r0*y*z/(2*r*r*r*r*r);
                 By_ext = b0*By_ext;
		 //assign according to orientation of spheromak
                 if(orient_sph==0){pfield->b.x2f(k,j,i)=By_ext;}
                 if(orient_sph==1){pfield->b.x2f(k,j,i)=By_ext;}
                 if(orient_sph==2){pfield->b.x3f(k,j,i)=By_ext;}//swap By and Bz in this case
		 //std::cout<<"outside spheromak radius By is"<<By_ext<<std::endl;
		 //if(pcoord->x3v(k)>0.0){pfield->b.x2f(k,j,i)=0.0;}
		 //else{pfield->b.x2f(k,j,i)=0;}
                        }
          }
        }
      }
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
		  //Bz
            Real x1 = pcoord->x1v(i);
            Real y1 = pcoord->x2v(j);
            Real z1 = pcoord->x3v(k);
            //adjust relative to spheromak center
            Real x=x1-x0;
            Real y=y1-y0;
            Real z=z1-z0+1e-15;
	    Real r0 = 4.493/Alpha; //define r0 in terms of Alpha
            std::cout<<"spheromake radius is "<<r0<<std::endl;
	    //ORIENTATION 0=z, 1=x, 2=y
            if(orient_sph==1){Real z_temp=z;z=x;x=z_temp;y=z_temp;}
	    //check if inside our outside spheromak
            Real r=std::sqrt(x*x+y*y+z*z);
            if (r<r0){//internal Bz
                 Real Bz_int = b0*((x*x+y*y-2*z*z)*std::cos(r*Alpha)/(r*r*r*r) + 
	              (2*z*z+x*x*x*x*Alpha*Alpha+y*y*y*y*Alpha*Alpha+y*y*(-1+z*z*Alpha*Alpha)+x*x*(-1+2*y*y*Alpha*Alpha+z*z*Alpha*Alpha))*std::sin(r*Alpha)/(Alpha*r*r*r*r*r));
		 //catch 0,0,0 case
		 if(r<1e-5){Bz_int = b0*(2/3)*Alpha*Alpha;}
		 //pfield->b.x3f(k,j,i)=Bz_int;
		 //assign according to orientation of spheromak
                 if(orient_sph==0){pfield->b.x3f(k,j,i)=Bz_int;}
                 if(orient_sph==1){pfield->b.x1f(k,j,i)=Bz_int;}//swap Bx and Bz for this case
                 if(orient_sph==2){pfield->b.x2f(k,j,i)=Bz_int;}//swap By and Bz for this case
		 //std::cout<<"inside spheromak radius Bz is"<<Bz_int<<std::endl;
                        }
            else{//external Bz
		 //placeholder
		 Real Bz_ext = 1+r0*r0*r0*(x*x+y*y+2*z*z)/(2*r*r*r*r*r);
                 Bz_ext = -1*b0*Bz_ext;
		 //assign according to orientation of spheromak
                 if(orient_sph==0){pfield->b.x3f(k,j,i)=Bz_ext;}
                 if(orient_sph==1){pfield->b.x1f(k,j,i)=Bz_ext;}//swap Bx and Bz for this case
                 if(orient_sph==2){pfield->b.x2f(k,j,i)=Bz_ext;}//swap By and Bz for this case
		 //pfield->b.x3f(k,j,i)=0.0;
                        }
          }
        }
      }
      if (NON_BAROTROPIC_EOS) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              phydro->u(IEN,k,j,i) += 0.5*b0*b0;
            }
          }
        }
      }
    }

    //initialize density (cell centered)
    grav_acc = phydro->hsrc.GetG3();
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real azero= pcoord->x3v(ke);//set exponential coefficient to z length of box
          //define x y and z
          Real x1 = pcoord->x1v(i);
          Real y1 = pcoord->x2v(j);
          Real z1 = pcoord->x3v(k);
	  //relative to spheromak center
          Real x=x1-x0;
          Real y=y1-y0;
          Real z=z1-z0+1e-15;//catch 0 case
	  //set density depending on inside or outside bubble
	  Real r0 = 4.493/Alpha;
	  Real den;
	  if (std::sqrt(x*x+y*y+z*z)<r0){//inside bubble
	  den=(1.451*1.4451*b0*b0*.5*beta_in);}//set using B at theta=pi/2,r=r0
	  else{den=azero*std::exp(grav_acc*z1);}//outside bubble
          //std::cout<< "attempting to set exponential density: " << den <<std::endl;
          if (pcoord->x3v(k) > 0.0) den *= drat;

          if (iprob == 1) {
            phydro->u(IM3,k,j,i) = (1.0 + std::cos(kx*(pcoord->x1v(i))))/8.0
                                   *(1.0 + std::cos(ky*pcoord->x2v(j)))
                                   *(1.0 + std::cos(kz*pcoord->x3v(k)));
          } else {
            phydro->u(IM3,k,j,i) = amp*(ran2(&iseed) - 0.5)*(
                1.0 + std::cos(kz*pcoord->x3v(k)));
          }

          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) *= (den*amp);
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = (1.0/gamma + grav_acc*den*(pcoord->x3v(k)))/gm1;
            phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/den;
          }
        }
      }
    }

    
  } // end of 3D initialization

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco,
//!                             AthenaArray<Real> &prim, FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        if (n==(IVY)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVY,k,jl-j,i) = -prim(IVY,k,jl+j-1,i);  // reflect 2-velocity
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,k,jl-j,i) = prim(IPR,k,jl+j-1,i)
                                 - prim(IDN,k,jl+j-1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,k,jl-j,i) = prim(n,k,jl+j-1,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(jl-j),i) =  b.x1f(k,(jl+j-1),i);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = -b.x2f(k,(jl+j  ),i);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) =  b.x3f(k,(jl+j-1),i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco,
//!                             AthenaArray<Real> &prim, FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        if (n==(IVY)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVY,k,ju+j,i) = -prim(IVY,k,ju-j+1,i);  // reflect 2-velocity
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,k,ju+j,i) = prim(IPR,k,ju-j+1,i)
                                 + prim(IDN,k,ju-j+1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,k,ju+j,i) = prim(n,k,ju-j+1,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(ju+j  ),i) =  b.x1f(k,(ju-j+1),i);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j+1),i) = -b.x2f(k,(ju-j+1),i);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j  ),i) =  b.x3f(k,(ju-j+1),i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco,
//!                             AthenaArray<Real> &prim, FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        if (n==(IVZ)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVZ,kl-k,j,i) = -prim(IVZ,kl+k-1,j,i);  // reflect 3-vel
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,kl-k,j,i) = prim(IPR,kl+k-1,j,i)
                                 - prim(IDN,kl+k-1,j,i)*grav_acc*(2*k-1)*pco->dx3f(k);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,kl-k,j,i) = prim(n,kl+k-1,j,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f((kl-k),j,i) =  b.x1f((kl+k-1),j,i);
        }
      }
    }

    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f((kl-k),j,i) =  b.x2f((kl+k-1),j,i);
        }
      }
    }

    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f((kl-k),j,i) = -b.x3f((kl+k  ),j,i);  // reflect 3-field
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco,
//!                             AthenaArray<Real> &prim, FaceField &b, Real time, Real dt,
//!                             int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        if (n==(IVZ)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVZ,ku+k,j,i) = -prim(IVZ,ku-k+1,j,i);  // reflect 3-vel
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,ku+k,j,i) = prim(IPR,ku-k+1,j,i)
                                 + prim(IDN,ku-k+1,j,i)*grav_acc*(2*k-1)*pco->dx3f(k);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,ku+k,j,i) = prim(n,ku-k+1,j,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b3
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f((ku+k  ),j,i) =  b.x1f((ku-k+1),j,i);
        }
      }
    }

    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f((ku+k  ),j,i) =  b.x2f((ku-k+1),j,i);
        }
      }
    }

    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f((ku+k+1),j,i) = -b.x3f((ku-k+1),j,i);  // reflect 3-field
        }
      }
    }
  }

  return;
}


// refinement condition: density jump
int RefinementCondition(MeshBlock *pmb) {
  int f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  AthenaArray<Real> &w = pmb->phydro->w;
  // maximum intercell density ratio
  Real drmax = 1.0;
  for (int k=pmb->ks-f3; k<=pmb->ke+f3; k++) {
    for (int j=pmb->js-f2; j<=pmb->je+f2; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        if (w(IDN,k,j,i-1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i-1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i+1)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j,i+1)/w(IDN,k,j,i);
        if (w(IDN,k,j,i)/w(IDN,k,j,i-1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i-1);
        if (w(IDN,k,j,i)/w(IDN,k,j,i+1) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j,i+1);
        if (f2) {
          if (w(IDN,k,j-1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j-1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j+1,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k,j+1,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k,j-1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j-1,i);
          if (w(IDN,k,j,i)/w(IDN,k,j+1,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k,j+1,i);
        }
        if (f3) {
          if (w(IDN,k-1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k-1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k+1,j,i)/w(IDN,k,j,i) > drmax) drmax = w(IDN,k+1,j,i)/w(IDN,k,j,i);
          if (w(IDN,k,j,i)/w(IDN,k-1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k-1,j,i);
          if (w(IDN,k,j,i)/w(IDN,k+1,j,i) > drmax) drmax = w(IDN,k,j,i)/w(IDN,k+1,j,i);
        }
      }
    }
  }
  if (drmax > 1.5) return 1;
  else if (drmax < 1.2) return -1;
  return 0;
}
