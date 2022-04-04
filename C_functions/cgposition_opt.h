#include "algebra.h"
typedef struct {
  double m,vold,volf;
  pts vel,veld,velf,rhov;
  double dissd,dissf;
  int nvof;
  pts Up_f;
}sta;


/** functoin that return 1 if the rank of the current process is 0 or if there is only one process**/


void Uprime_f(sta * st){
  double xp,yp;
  xp=yp=0;
  foreach(reduction(+:xp) reduction(+:yp)){
    xp += (u.x[] - st->velf.x) * dv() * (1 - f[]); 
    yp += (u.y[] - st->velf.y) * dv() * (1 - f[]); 
  }
  st->Up_f.x = xp/st->volf; 
  st->Up_f.y = yp/st->volf; 
  #if dimension == 3
  double zp=0;
  foreach(reduction(+:zp) ){
    zp += (u.z[] - st->velf.z) * dv() * (1 - f[]); 
  }
  st->Up_f.z = zp/st->volf; 
  #endif 
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

double avg_vy(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv()*(f[]*(rho1-rho2)+rho2);
      rhv += dv()*(f[]*(rho1-rho2)+rho2)*u.y[];
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}
double avg_vx(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv()*(f[]*(rho1-rho2)+rho2);
      rhv += dv()*(f[]*(rho1-rho2)+rho2)*u.x[];
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}
#if dimension == 3
double avg_vz(){
  double rhv=0.,volume=0.;
  foreach(reduction(+:rhv) reduction(+:volume)){
    if(f[]!=nodata && dv()>0.){
      volume +=dv()*(f[]*(rho1-rho2)+rho2);
      rhv += dv()*(f[]*(rho1-rho2)+rho2)*u.z[];
    }
  }
  rhv= volume ? rhv/volume : 0.;
  return rhv;
}
#endif


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
  foreach(reduction(+:vxd) reduction(+:vyd) reduction(+:vxf) reduction(+:vyf) reduction(+:rhovx) reduction(+:rhovy) reduction(+:m) reduction(+:vtd) reduction(+:vtf) ){
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
    //fluid and drops violume
    vtd+=dv()*f[];
    vtf+=dv()*(1-f[]);
    m+=dm;
  };
  double Vx,Vy;
  Vx=Vy=0;
  foreach(reduction(+:Vx) reduction(+:Vy)){
    Vx+=dv()*u.x[];
    Vy+=dv()*u.y[];
  }
  st->vel.x   = Vx/(Ls*Ls);
  st->vel.y   = Vy/(Ls*Ls);
  st->veld.x  = vxd/vtd;
  st->veld.y  = vyd/vtd;
  st->velf.x  = vxf/vtf;
  st->velf.y  = vyf/vtf;
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