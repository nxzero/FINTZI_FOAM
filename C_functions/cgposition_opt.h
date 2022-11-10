#include "algebra.h"

typedef struct {
  double m,vold,volf;
  coord V,Vd,Vf,RhoV;
  double dissd,dissf;
  int nvof;
  coord Upf,Upp;
}sta;


/** functoin that return 1 if the rank of the current process is 0 or if there is only one process**/

void copy_file(char * sname,char * tname, FILE * target){
  char ch;
  FILE * Source = fopen(sname,"r");
  target = fopen(tname,"w");
  while( ( ch = fgetc(Source) ) != EOF )
    fputc(ch, target); 
  fclose(Source);
  fclose(target);
}

void Uprime_f(sta * st){
  coord Upf={0},Upp={0};
  foreach(reduction(+:Upf) reduction(+:Upp)){
    foreach_dimension(){
      Upf.x += sq(u.x[] - st->Vf.x) * dv() * (1 - f[]); 
      Upp.x += sq(u.x[] - st->Vd.x) * dv() * f[]; 
    }
  }
  foreach_dimension(){
    st->Upf.x = Upf.x/st->volf;  
    st->Upp.x = Upp.x/st->vold; 
  }
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


coord avg_U(){
  coord U={0};
  double V=0;
  foreach(reduction(+:U) reduction(+:V)){
    foreach_dimension()
      U.x += u.x[]*dv();
    V += dv();
  }
  foreach_dimension() 
    U.x = U.x/V;
  return U;
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
  double m=0,vold=0,volf=0,vol=0;
  coord Vd={0},Vf={0},RhoV={0},V={0};
  // calcul des vol moy , vel moy and momentum moy 
  foreach(reduction(+:Vd) reduction(+:Vf) reduction(+:RhoV) reduction(+:m) reduction(+:vold) reduction(+:volf)){
    double dm=dv()*(f[]*(rho1-rho2)+rho2);
    m +=dm;
    vold+=dv()*f[];
    volf+=dv()*(1-f[]);
    foreach_dimension(){
      Vd.x  +=dv()*f[]*u.x[];//average velocity in the drops
      Vf.x  +=dv()*(1-f[])*u.x[];//average velocity in the fluid
      RhoV.x+=dm*u.x[];//average momentum
    }
  };
  foreach(reduction(+:V) reduction(+:vol)){
    foreach_dimension()
      V.x +=dv()*u.x[];
      vol += dv();
  }
  foreach_dimension(){
    st->V.x   = V.x/vol;
    st->Vd.x  = Vd.x/vold;
    st->Vf.x  = Vf.x/volf;
    st->RhoV.x  = RhoV.x/m;
  }
  st->m       = m;
  st->vold    = vold;
  st->volf    = volf;
}

sta calcul_sta(){
  sta st;
  calcul_vel_and_rhov(&st);
  Uprime_f(&st);
  dissipation_rate(&st);
  st.nvof=0;
  for(scalar s in interfaces ){
    st.nvof+=1;
  }
  // st.nvof = nvof;
  return st;
}