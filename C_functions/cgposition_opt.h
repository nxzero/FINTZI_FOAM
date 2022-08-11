#include "algebra.h"

typedef struct {
  double m,vold,volf;
  pts vel,veld,velf,rhov;
  double dissd,dissf,pf,pfstd;
  int nvof;
  pts Up_f;
}sta;


/** functoin that return 1 if the rank of the current process is 0 or if there is only one process**/


void Uprime_f(sta * st){
  #if dimension == 2 
  double Up_x=0,Up_y=0,pfstd=0;
  foreach(reduction(+:Up_x) reduction(+:Up_y)){
  #elif dimension == 3 
  double Up_z=0,Up_x=0,Up_y=0,pfstd=0;
  foreach(reduction(+:Up_x) reduction(+:Up_y) reduction(+:Up_z)){
  #endif
    foreach_dimension()
      Up_x += sq(u.x[] - st->velf.x) * dv() * (1 - f[]); 
      // Up_x_y += (u.x[] - st->velf.x)*(u.y[] - st->velf.y) * dv() * (1 - f[]); 
      pfstd += sq(p[] - st->pf) * dv() * (1 - f[]); 
  }
  st->pfstd = pfstd/st->volf;
  foreach_dimension()
    st->Up_f.x = Up_x/st->volf; 
    // st->Up_f.x.x = Up_x/st->volf; 
    // st->Up_f.x.x = Up_x/st->volf; 
}

double avg_rho(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv();
      rhv += dv()*(f[]*(rho1-rho2)+rho2);
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}


void dissipation_rate (sta * st)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    #if dimension == 2
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformyx) +sq(SDeformyy));
    #elif dimension == 3 
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    #endif
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  st->dissd = rateWater/st->vold;
  st->dissf = rateAir/st->volf;
}

void calcul_vel_and_rhov(sta * st){
  double vxd,vyd,vxf,vyf,rhovx,rhovy,m,vtd,vtf;
  vxd=vyd=vxf=vyf=rhovx=rhovy=m=vtd=vtf=0;
  // calcul des vol moy , vel moy and momentum moy 
  foreach(reduction(+:vxd) reduction(+:vyd) reduction(+:vxf) reduction(+:vyf) 
          reduction(+:rhovx) reduction(+:rhovy) reduction(+:m) reduction(+:vtd) reduction(+:vtf) ){
    double dm=dv()*(f[]*(rho1-rho2)+rho2);
    //average velocity in the drops
    vxd+=dv()*f[]*u.x[];
    vyd+=dv()*f[]*u.y[];
    //average velocity in the fluid
    vxf+=dv()*(1-f[])*u.x[];
    vyf+=dv()*(1-f[])*u.y[];
    //average momentum
    rhovx+=dm*u.x[];
    rhovy+=dm*u.y[];
    //fluid and drops volume
    vtd+=dv()*f[];
    vtf+=dv()*(1-f[]);
    m+=dm;
  };
  double Vx,Vy,pf;
  Vx=Vy=pf=0;
  foreach(reduction(+:Vx) reduction(+:Vy) reduction(+:pf)){
    Vx+=dv()*u.x[];
    Vy+=dv()*u.y[];
    pf+=dv()*p[]*(1-f[]);
  }
  st->vel.x   = Vx/(Ls*Ls);
  st->vel.y   = Vy/(Ls*Ls);
  st->veld.x  = vxd/vtd;
  st->veld.y  = vyd/vtd;
  st->velf.x  = vxf/vtf;
  st->velf.y  = vyf/vtf;
  st->pf    = pf/vtf;
  st->rhov.x  = rhovx/m;
  st->rhov.y  = rhovy/m;
  st->m       = m;
  st->vold    = vtd;
  st->volf    = vtf;
  #if dimension == 3
  double vzd,vzf,rhovz;
  double Vz;
  vzd=vzf=rhovz=Vz=0;
  foreach(reduction(+:vzd) reduction(+:vzf) reduction(+:rhovz) reduction(+:Vz)){
    double dm=dv()*(f[]*(rho1-rho2)+rho2);
    vzd+=dv()*f[]*u.z[];
    vzf+=dv()*(1-f[])*u.z[];
    rhovz+=dm*u.z[];
    Vz+=dv()*u.z[];
  }  
  st->vel.z   = Vz/(Ls*Ls);
  st->veld.z  = vzd/vtd;
  st->velf.z  = vzf/vtf;
  st->rhov.z  = rhovz/m;
  #endif
}

sta calcul_sta(){
  sta st;
  calcul_vel_and_rhov(&st);
  Uprime_f(&st);
  dissipation_rate(&st);
  int nvof=0;
  for(scalar s in interfaces ){
    reduction(+:nvof);
    nvof+=1;
  }
  st.nvof = nvof;
  return st;
}