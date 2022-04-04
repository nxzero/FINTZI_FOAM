#include "algebra.h"
typedef struct {
  double omegaz,Jz,vol,diss,distmin,defnorm;
  pts EigVal,pos,vel,omega,J;
  ten EigVec,Sig;
  int tagmin,xper,yper,zper,tag;
} cg;

/** functoin that return 1 if the rank of the current process is 0 or if there is only one process**/

// check if the drop is perio or not 
#if dimension == 2
void isperio(scalar s, scalar tag ,int tg,int tagi, cg n[] ){
  int xper1[tg],xper2[tg],yper1[tg],yper2[tg];
  for (int i = 0; i < tg; i++) yper2[i]=yper1[i]=xper1[i]=xper2[i]=0;
  foreach(reduction(max:xper1[:tg]) reduction(max:xper2[:tg]) reduction(max:yper1[:tg]) reduction(max:yper2[:tg])){
    if (s[] != nodata && dv() > 0.  && s[]>=1e-4) {
      int i=tag[]-1;
      if(x>D) xper1[i]=1;
      else if(x<-D) xper2[i]=1;
      if(y>D) yper1[i]=1;
      else if(y<-D) yper2[i]=1;
    }
  }
  for (int i = 0; i < tg; i++){
    n[tagi+i].xper = (xper1[i]&&xper2[i]);
    n[tagi+i].yper = (yper1[i]&&yper2[i]);
  }
}
#elif dimension == 3
void isperio(scalar s, scalar tag ,int tg,int tagi, cg n[]){
  int xper1[tg],xper2[tg],yper1[tg],yper2[tg],zper1[tg],zper2[tg];
  for (int i = 0; i < tg; i++) yper2[i]=yper1[i]=xper1[i]=xper2[i]=zper1[i]=zper2[i]=0;
  foreach(reduction(max:xper1[:tg]) reduction(max:xper2[:tg]) reduction(max:yper1[:tg]) reduction(max:yper2[:tg]) reduction(max:zper1[:tg]) reduction(max:zper2[:tg])){
    if (s[] != nodata && dv() > 0.  && s[]>=1e-4) {
      int i=tag[]-1;
      if(x>D) xper1[i]=1;
      else if(x<-D) xper2[i]=1;
      if(y>D) yper1[i]=1;
      else if(y<-D) yper2[i]=1;
      if(z>D) zper1[i]=1;
      else if(z<-D) zper2[i]=1;
    }
  }
  for (int i = 0; i < tg; i++){
    n[tagi+i].xper = (xper1[i]&&xper2[i]);
    n[tagi+i].yper = (yper1[i]&&yper2[i]);
    n[tagi+i].zper = (zper1[i]&&zper2[i]);
  }
}
#endif
#if dimension == 2
void pos_vol_and_vel(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double xpos[tg], ypos[tg],volume[tg],vx[tg],vy[tg];
  for (int i = 0; i < tg; i++) xpos[i]= ypos[i]=volume[i]=vx[i]=vy[i]=0;
  foreach(reduction(+:xpos[:tg]) reduction(+:ypos[:tg]) reduction(+:volume[:tg]) reduction(+:vx[:tg]) reduction(+:vy[:tg])){
    if (s[] != nodata && dv() > 0. && s[]>=1e-4){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      if(n[tagi+i].xper)
	      xpos[i]    += dv()*s[]*(x>0?x:x+L0);
      else
	      xpos[i]    += dv()*s[]*x;
      if(n[tagi+i].yper)
	      ypos[i]    += dv()*s[]*(y>0?y:y+L0);
      else
	      ypos[i]    += dv()*s[]*y;
      vx[i]   += dv()*s[]*u.x[];
      vy[i]   += dv()*s[]*u.y[];
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].pos.x = volume[i] ? xpos[i]/volume[i] : 0.;
    n[tagi+i].pos.y = volume[i] ? ypos[i]/volume[i] : 0.;
    n[tagi+i].vel.x= volume[i] ? vx[i]/volume[i] : 0.;
    n[tagi+i].vel.y= volume[i] ? vy[i]/volume[i] : 0.;
    n[tagi+i].vol = volume[i];
  }
}
#elif dimension == 3
void pos_vol_and_vel(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double xpos[tg], ypos[tg], zpos[tg],volume[tg],vx[tg],vy[tg],vz[tg];
  for (int i = 0; i < tg; i++) xpos[i]= ypos[i]= zpos[i]=volume[i]=vx[i]=vy[i]=vz[i]=0;
  foreach(reduction(+:xpos[:tg]) reduction(+:ypos[:tg]) reduction(+:zpos[:tg]) reduction(+:volume[:tg]) reduction(+:vx[:tg]) reduction(+:vy[:tg]) reduction(+:vz[:tg])){
    if (s[] != nodata && dv() > 0. && s[]>=1e-4){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      if(n[tagi+i].xper)
	      xpos[i]    += dv()*s[]*(x>0?x:x+L0);
      else
	      xpos[i]    += dv()*s[]*x;
      if(n[tagi+i].yper)
	      ypos[i]    += dv()*s[]*(y>0?y:y+L0);
      else
	      ypos[i]    += dv()*s[]*y;
      if(n[tagi+i].zper)
	      zpos[i]    += dv()*s[]*(z>0?z:z+L0);
      else
	      zpos[i]    += dv()*s[]*z;
      vx[i]   += dv()*s[]*u.x[];
      vy[i]   += dv()*s[]*u.y[];
      vz[i]   += dv()*s[]*u.z[];
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].pos.x = volume[i] ? xpos[i]/volume[i] : 0.;
    n[tagi+i].pos.y = volume[i] ? ypos[i]/volume[i] : 0.;
    n[tagi+i].pos.z = volume[i] ? zpos[i]/volume[i] : 0.;
    n[tagi+i].vel.x= volume[i] ? vx[i]/volume[i] : 0.;
    n[tagi+i].vel.y= volume[i] ? vy[i]/volume[i] : 0.;
    n[tagi+i].vel.z= volume[i] ? vz[i]/volume[i] : 0.;
    n[tagi+i].vol = volume[i];
  }
}
#endif
#if dimension == 2
void omega(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double omegaz[tg],Jz[tg];
  for (int i = 0; i < tg; i++) omegaz[i]=Jz[i]=0;
  foreach(reduction(+:omegaz[:tg]) reduction(+:Jz[:tg])){
    if (s[] != nodata && dv() > 0. && s[]>=1e-4){
      int i=tag[]-1;
      pts posG = {n[tagi+i].pos.x,n[tagi+i].pos.y};
      pts pos = {x,y};
      if(n[tagi+i].xper)
	      pos.x    = x>0?x:x+L0;
      if(n[tagi+i].yper)
	      pos.y    = y>0?y:y+L0;
      pts lever_arm = diff_pts(pos,posG);
			Jz[i] += dv()*s[]*pow(normL2_pts(lever_arm),2);
      pts vel = {u.x[],u.y[]};
      // projection on the perpendicular axis zz
      omegaz[i] +=  dv()*s[]*cross_pts(lever_arm,vel); 
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].Jz = n[tagi+i].vol ? Jz[i]/n[tagi+i].vol : 0.;
    n[tagi+i].omegaz = Jz[i] ? omegaz[i]/Jz[i] : 0.;
  }
}
#elif dimension == 3
void omega(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double omegax[tg],omegay[tg],omegaz[tg],Jx[tg],Jy[tg],Jz[tg];
  for (int i = 0; i < tg; i++) omegax[i]=omegay[i]=omegaz[i]=Jx[i]=Jy[i]=Jz[i]=0;
  pts omega = {0,0,0};
  foreach(reduction(+:omegaz[:tg])reduction(+:omegay[:tg])reduction(+:omegax[:tg]) reduction(+:Jz[:tg])reduction(+:Jy[:tg])reduction(+:Jx[:tg])){
    if (s[] != nodata && dv() > 0. && s[]>=1e-4){
      int i=tag[]-1;
      pts posG = {n[tagi+i].pos.x,n[tagi+i].pos.y,n[tagi+i].pos.z};
      pts pos = {x,y,z};
      if(n[tagi+i].xper)
	      pos.x    = x>0?x:x+L0;
      if(n[tagi+i].yper)
	      pos.y    = y>0?y:y+L0;
      if(n[tagi+i].zper)
	      pos.z    = z>0?z:z+L0;
      pts r = diff_pts(pos,posG);
			pts lever_armx ={0,r.y,r.z};
			pts lever_army ={r.x,0,r.z};
			pts lever_armz ={r.x,r.y,0};
			Jx[i] += dv()*s[]*pow(normL2_pts(lever_armx),2);
			Jy[i] += dv()*s[]*pow(normL2_pts(lever_army),2);
			Jz[i] += dv()*s[]*pow(normL2_pts(lever_armz),2);
      pts vel = {u.x[],u.y[],u.z[]};
      // projection on the perpendicular axis 
      omega =  mult_pts(cross_pts(r,vel),dv()*s[]); 

      omegaz[i] += omega.z;
      omegay[i] += omega.y;
      omegax[i] += omega.x;
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].omega.z = Jz[i] ? omegaz[i]/Jz[i] : 0.;
    n[tagi+i].J.z = n[tagi+i].vol ? omegaz[i]/n[tagi+i].vol : 0.;
    n[tagi+i].omega.x = Jx[i] ? omegax[i]/Jx[i] : 0.;
    n[tagi+i].J.x = n[tagi+i].vol ? omegax[i]/n[tagi+i].vol : 0.;
    n[tagi+i].omega.y = Jy[i] ? omegay[i]/Jy[i] : 0.;
    n[tagi+i].J.y = n[tagi+i].vol ? omegay[i]/n[tagi+i].vol : 0.;
  }
}
#endif 

#if dimension == 2
void deform(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double xx[tg],yy[tg],xy[tg];
  for (int i = 0; i < tg; i++) xx[i]=yy[i]=xy[i]=0;
  foreach(reduction(+:xx[:tg]) reduction(+:yy[:tg]) reduction(+:xy[:tg])){
    if (s[] != nodata && dv() > 0. &&s[]>=1e-4) {
      int i=tag[]-1;
      if(n[tagi+i].xper) 
        xx[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x));
      else
	      xx[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(x-n[tagi+i].pos.x);
      if(n[tagi+i].yper)
	      yy[i]    += dv()*s[]*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y))*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
      else
	      yy[i]    += dv()*s[]*(y-n[tagi+i].pos.y)*(y-n[tagi+i].pos.y);
      
      if(n[tagi+i].xper&&n[tagi+i].yper)
	      xy[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
	    else if(n[tagi+i].xper)
	      xy[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(y-n[tagi+i].pos.y);
      else if(n[tagi+i].yper)
	      xy[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
      else
	      xy[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(y-n[tagi+i].pos.y);
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].Sig.x.x = n[tagi+i].vol ? xx[i]/n[tagi+i].vol : 0.;
    n[tagi+i].Sig.y.y = n[tagi+i].vol ? yy[i]/n[tagi+i].vol : 0.;
    n[tagi+i].Sig.x.y = n[tagi+i].vol ? xy[i]/n[tagi+i].vol : 0.;
    n[tagi+i].Sig.y.x = n[tagi+i].vol ? xy[i]/n[tagi+i].vol : 0.;
    n[tagi+i].EigVal = EigenValue(n[tagi+i].Sig);
    n[tagi+i].EigVec = EigenVectors(n[tagi+i].Sig);
    n[tagi+i].defnorm = normL2_ten(n[tagi+i].Sig);
  }
}
#elif dimension == 3
void deform(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double xx[tg],yy[tg],zz[tg],xy[tg],zy[tg],zx[tg];
  for (int i = 0; i < tg; i++) xx[tg]=yy[tg]=zz[tg]=xy[tg]=zy[tg]=zx[tg]=0;
  foreach(reduction(+:xx[:tg]) reduction(+:yy[:tg]) reduction(+:xy[:tg]) reduction(+:zz[:tg]) reduction(+:zy[:tg]) reduction(+:zx[:tg]) ){
    if (s[] != nodata && dv() > 0. &&s[]>=1e-4) {
      int i=tag[]-1;

      if(n[tagi+i].xper) 
        xx[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x));
      else
	      xx[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(x-n[tagi+i].pos.x);
      if(n[tagi+i].yper)
	      yy[i]    += dv()*s[]*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y))*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
      else
	      yy[i]    += dv()*s[]*(y-n[tagi+i].pos.y)*(y-n[tagi+i].pos.y);
      if(n[tagi+i].zper)
	      zz[i]    += dv()*s[]*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z))*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z));
      else
	      zz[i]    += dv()*s[]*(z-n[tagi+i].pos.z)*(z-n[tagi+i].pos.z);
      
      if(n[tagi+i].xper&&n[tagi+i].yper)
	      xy[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
	    else if(n[tagi+i].xper)
	      xy[i]    += dv()*s[]*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x))*(y-n[tagi+i].pos.y);
      else if(n[tagi+i].yper)
	      xy[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
      else
	      xy[i]    += dv()*s[]*(x-n[tagi+i].pos.x)*(y-n[tagi+i].pos.y);

      if(n[tagi+i].zper&&n[tagi+i].yper)
	      zy[i]    += dv()*s[]*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z))*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
	    else if(n[tagi+i].zper)
	      zy[i]    += dv()*s[]*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z))*(y-n[tagi+i].pos.y);
      else if(n[tagi+i].yper)
	      zy[i]    += dv()*s[]*(z-n[tagi+i].pos.z)*(y>0?(y-n[tagi+i].pos.y):(y+L0-n[tagi+i].pos.y));
      else
	      zy[i]    += dv()*s[]*(z-n[tagi+i].pos.z)*(y-n[tagi+i].pos.y);

      if(n[tagi+i].zper&&n[tagi+i].xper)
	      zx[i]    += dv()*s[]*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z))*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x));
	    else if(n[tagi+i].zper)
	      zx[i]    += dv()*s[]*(z>0?(z-n[tagi+i].pos.z):(z+L0-n[tagi+i].pos.z))*(x-n[tagi+i].pos.x);
      else if(n[tagi+i].xper)
	      zx[i]    += dv()*s[]*(z-n[tagi+i].pos.z)*(x>0?(x-n[tagi+i].pos.x):(x+L0-n[tagi+i].pos.x));
      else
	      zx[i]    += dv()*s[]*(z-n[tagi+i].pos.z)*(x-n[tagi+i].pos.x);
    }
  }
  for(int i=0;i<tg;i++){    
      n[tagi+i].Sig.x.x = n[tagi+i].vol ? xx[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.y.y = n[tagi+i].vol ? yy[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.x.y = n[tagi+i].vol ? xy[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.y.x = n[tagi+i].vol ? xy[i]/n[tagi+i].vol : 0.;
  
      n[tagi+i].Sig.z.z = n[tagi+i].vol ? zz[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.z.x = n[tagi+i].vol ? zx[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.x.z = n[tagi+i].vol ? zx[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.z.y = n[tagi+i].vol ? zy[i]/n[tagi+i].vol : 0.;
      n[tagi+i].Sig.y.z = n[tagi+i].vol ? zy[i]/n[tagi+i].vol : 0.;
      n[tagi+i].defnorm = normL2_ten(n[tagi+i].Sig);
      // n[tagi+i].EigVal = EigenValue(n[tagi+i].Sig);
      // n[tagi+i].EigVec = EigenVectors(n[tagi+i].Sig);
  }
}
#endif 

void dissipation_rate_on_bubbles (scalar s, scalar tag ,int tg,int tagi,cg n[])
{
  double diss[tg];
  for (int i = 0; i < tg; i++) diss[i]=0;
  foreach (reduction (+:diss[:tg])) {
    if (s[] != nodata && dv() > 0.  && s[]>=1e-4) {
      int i=tag[]-1;
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
      diss[i] += mu1/rho[]*f[]*sqterm; //water
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].diss = n[tagi+i].vol ? diss[i]/n[tagi+i].vol : 0.;
  }
}

// main functions
void cg_bub (scalar s, scalar tag ,int tg,int tagi,cg n[])
{
  isperio(s,tag,tg,tagi,n);
  pos_vol_and_vel(s,tag,tg,tagi,n);
  omega(s,tag,tg,tagi,n);
  deform(s,tag,tg,tagi,n);
  //interfacial NRJ ?? 
  dissipation_rate_on_bubbles(s,tag,tg,tagi,n);
  for(int i=0;i<tg;i++){
    if (n[tagi+i].pos.x>L0/2.) n[tagi+i].pos.x-=L0;
    if (n[tagi+i].pos.y>L0/2.) n[tagi+i].pos.y-=L0;
    #if dimension == 3
    if (n[tagi+i].pos.z>L0/2.) n[tagi+i].pos.z-=L0;
    #endif
  }
}





//dimension or "size" from the center, this have to be the same everywere 
double dist_perio(pts a,pts b){
  if(a.x>0){ 
    b.x = (b.x<a.x-L0/2)*L0+b.x;
  }else{
    b.x = (b.x>a.x+L0/2)*(-L0)+b.x;
  }
  if(a.y>0){ 
    b.y = (b.y<a.y-L0/2)*L0+b.y;
  }else{
    b.y = (b.y>a.y+L0/2)*(-L0)+b.y;
  }
  #if dimension == 3
  if(a.z>0){ 
    b.z = (b.z<a.z-L0/2)*L0+b.z;
  }else{
    b.z = (b.z>a.z+L0/2)*(-L0)+b.z;
  }
  #endif

  return normL2_pts(diff_pts(a,b));
}

// // for the linear algebra op include the lib but for now do it by hands
void assign_tags_from_last_step (cg datat0[],cg datat1[],int tagi){
//   // the first step is to calcul the prediction of the particlue at time n-1 pts of the 
  for (int j = 0; j < tagi; j++)
  {
    pts np = datat0[j].pos;
    pts vp = datat0[j].vel;
    double dt = dtprint;
    np = add_pts(np,mult_pts(vp,dt)); // euleur schem
    double r = sqrt(datat0[j].vol/M_PI)*0.5; // 
    bool probleme_of_assignement = true;
    for (int k = j; k < (tagi+j); k++)
    {
      int idx = k % tagi;
      pts np1 = datat1[idx].pos;
      double dist = dist_perio(np,np1);
      if(dist<r){
        datat1[idx].tag = datat0[j].tag;
        probleme_of_assignement = false;
        break;
      }
    }
    if(probleme_of_assignement){
        fprintf(stdout,"Error of assignement at t = %g for bubble tag = %d\n",t,datat0[j].tag);
    }
  }
}

void print_pair_dist(cg data[],int tagi, FILE * Dist,int step,float t){
  // put the bubbles in order of tags 
  cg data_orderred[tagi];
  for (int j = 0; j < tagi; j++)
  {
    data[j].distmin = Ls*10;
    data_orderred[data[j].tag]=data[j];
  }
  Dist = fopen("Dist.csv","a");
  fprintf(Dist,"%d,%g",step,t);
  for (int k = 0; k < tagi; k++)
  {
    pts npi = data_orderred[k].pos; 
    for (int j = k+1; j<tagi;j++){
      pts npj = data_orderred[j].pos;
      double dist = dist_perio(npi,npj);
      fprintf(Dist,",%g",dist);
      if(dist<data_orderred[k].distmin){
        data_orderred[k].distmin = dist;
        data_orderred[k].tagmin = j;
      } 
      if(dist < data_orderred[j].distmin){
        data_orderred[j].distmin = dist;
        data_orderred[j].tagmin = k;
      }
    }
    for (int j = 0; j < tagi; j++)
    {
      if(data[j].tag == k) {
        data[j].distmin = data_orderred[k].distmin;
        data[j].tagmin = data_orderred[k].tagmin;
      }
    }
  }
  fprintf(Dist,"\n");
  fclose(Dist);
}
