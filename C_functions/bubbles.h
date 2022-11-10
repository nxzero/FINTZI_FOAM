#include "algebra.h"
#include "tag.h"
#define POS {x,y,z} 
#define P x,y,z
FILE * Dist_x;
char name_x[] = "Dist_x.csv";
char name_y[] = "Dist_y.csv";
FILE * Dist_y;
#if dimension == 3
FILE * Dist_z;
char name_z[] = "Dist_z.csv";
#endif
typedef struct {
  double omega_z,J_z,J_y,J_x,J_x_y,vol,diss,defnorm,ct;
  coord EigVal,pos,vel,velc,distmin,per;
  ten EigVec,Sig;
  int tagmin,tag,tagmax;
  double pressure;
  int realtag,si;
  int adj[100];// maximum index of scalar feilds  
  #if dimension == 3
  double omega_y,omega_x,J_y_z,J_z_x;
  int per_z;
  #endif
} cg;

/** functoin that return 1 if the rank of the current process is 0 or if there is only one process**/


coord POS_perio(coord pos,coord per){
  foreach_dimension(){
    if(per.x) pos.x = pos.x>0?pos.x:pos.x+L0;
  }
  return pos;
}
double dist_perio(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-L0/2)*L0+b.x;
    else b.x = (b.x>a.x+L0/2)*(-L0)+b.x;
  }
  return normL2_pts(diff_pts(a,b));
}
coord dist_perioPts(coord a,coord b){
  foreach_dimension(){
    if(a.x>0) b.x = (b.x<a.x-L0/2)*L0+b.x;
    else b.x = (b.x>a.x+L0/2)*(-L0)+b.x;
  }
  return diff_pts(a,b);
}
// check if the drop is perio or not 
void isperio(scalar s, scalar tag ,int tg,int tagi, cg n[] ){
  coord per1[tg],per2[tg];
  for (int i = 0; i < tg; i++) 
    foreach_dimension() {
      per1[i].x = 0;
      per2[i].x = 0;
    }
  foreach(reduction(max:per1[:tg]) reduction(max:per2[:tg]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      coord pos = POS;
      foreach_dimension(){
        if(pos.x>(L0/2-D)) per1[i].x=1;
        else if(pos.x<-(L0/2-D)) per2[i].x=1;
      }
    }
  for (int i = 0; i < tg; i++) 
    foreach_dimension()
      n[tagi+i].per.x = (per1[i].x && per2[i].x);
}

void pos_vol_and_vel(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  double pressure[tg];
  #if dimension == 2
  double pos_x[tg], pos_y[tg],volume[tg],v_x[tg],v_y[tg];
  for (int i = 0; i < tg; i++) pos_x[i]= pos_y[i]=volume[i]=v_x[i]=v_y[i]=0;
  foreach(reduction(+:pos_x[:tg]) reduction(+:pos_y[:tg])
          reduction(+:v_x[:tg])   reduction(+:v_y[:tg])   
          reduction(+:volume[:tg]) reduction(+:pressure[:tg])){
  #elif dimension == 3
  double pos_x[tg], pos_y[tg], pos_z[tg],volume[tg],v_x[tg],v_y[tg],v_z[tg];
  for (int i = 0; i < tg; i++) pos_x[i]= pos_y[i]= pos_z[i]=volume[i]=v_x[i]=v_y[i]=v_z[i]=0;
  foreach(reduction(+:pos_x[:tg])  reduction(+:pos_y[:tg]) reduction(+:pos_z[:tg]) 
          reduction(+:v_x[:tg])    reduction(+:v_y[:tg])   reduction(+:v_z[:tg])
          reduction(+:volume[:tg]) reduction(+:pressure[:tg])){
  #endif
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      volume[i] += dv()*s[];
      pressure[i] += p[]*s[]*dv();
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].per);
      foreach_dimension(){
        pos_x[i]    += dv()*s[]*pos.x;
        v_x[i]   += dv()*s[]*u.x[];
      }
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].pressure = volume[i] ? pressure[i]/volume[i] : 0.;
    foreach_dimension(){
      n[tagi+i].pos.x = volume[i] ? pos_x[i]/volume[i] : 0.;
      n[tagi+i].vel.x = volume[i] ? v_x[i]/volume[i] : 0.;
    }
    n[tagi+i].vol = volume[i];
  }
}



void omega(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  #if dimension == 2
  double omega_z[tg],J_z[tg],J_y[tg],J_x[tg],J_x_y[tg] ;
  for (int i = 0; i < tg; i++) omega_z[i]=J_z[i]=J_y[i]=J_x[i]=J_x_y[tg]=0;
  foreach(reduction(+:omega_z[:tg]) reduction(+:J_z[:tg]) reduction(+:J_x_y[:tg]) 
  reduction(+:J_y[:tg])  reduction(+:J_x[:tg]) ){
  #elif dimension == 3
  double omega_x[tg],omega_y[tg],omega_z[tg],J_z[tg],J_y[tg],J_x[tg],J_z_x[tg],J_y_z[tg],J_x_y[tg];
  for (int i = 0; i < tg; i++) omega_x[i]=omega_y[i]=omega_z[i]=J_z[i]=J_y[i]=J_x[i]=J_z_x[i]=J_y_z[i]=J_x_y[i]=0;
  foreach(reduction(+:omega_z[:tg]) reduction(+:J_z[:tg]) reduction(+:J_x_y[:tg])
          reduction(+:omega_y[:tg]) reduction(+:J_y[:tg]) reduction(+:J_y_z[:tg])
          reduction(+:omega_x[:tg]) reduction(+:J_x[:tg]) reduction(+:J_z_x[:tg])){
  #endif
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].per);
      coord posP = diff_pts(pos,n[tagi+i].pos);
      J_z[i] += dv()*s[]*(pow(posP.x,2)+pow(posP.y,2));
      omega_z[i] +=  dv()*s[]*(posP.x * u.y[] - posP.y * u.x[]); 
      J_x_y[i] -= dv()*s[]*posP.x*posP.y;
      #if dimension == 2
      J_x[i] += dv()*s[]*pow(posP.y,2);
      J_y[i] += dv()*s[]*pow(posP.x,2);
      #elif dimension == 3
      J_x[i] -= dv()*s[]*(pow(posP.y,2)+pow(posP.z,2));
      J_y[i] -= dv()*s[]*(pow(posP.z,2)+pow(posP.x,2));
      J_z_x[i] += dv()*s[]*posP.z*posP.x;
      J_y_z[i] += dv()*s[]*posP.y*posP.z;
      omega_x[i] +=  dv()*s[]*(posP.y * u.z[] - posP.z * u.y[]); 
      omega_y[i] +=  dv()*s[]*(posP.z * u.x[] - posP.x * u.z[]); 
      #endif
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].J_x = n[tagi+i].vol ? J_x[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_y = n[tagi+i].vol ? J_y[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_z = n[tagi+i].vol ? J_z[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_x_y = n[tagi+i].vol ? J_x_y[i]/n[tagi+i].vol : 0.;
    n[tagi+i].omega_z = J_z[i] ? omega_z[i]/J_z[i] : 0.;
    #if dimension == 2
    n[tagi+i].defnorm = sqrt(2*sq(n[tagi+i].J_x_y)+sq(n[tagi+i].J_x)+sq(n[tagi+i].J_y)+sq(n[tagi+i].J_z)); 
    #elif dimension == 3
    n[tagi+i].omega_x = J_x[i] ? omega_x[i]/J_x[i] : 0.;
    n[tagi+i].omega_y = J_y[i] ? omega_y[i]/J_y[i] : 0.;
    n[tagi+i].J_y_z = n[tagi+i].vol ? J_y_z[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_z_x = n[tagi+i].vol ? J_z_x[i]/n[tagi+i].vol : 0.;
    n[tagi+i].defnorm =sqrt(2*sq(n[tagi+i].J_z_x)+2*sq(n[tagi+i].J_y_z)+2*sq(n[tagi+i].J_x_y)+sq(n[tagi+i].J_x)+sq(n[tagi+i].J_y)+sq(n[tagi+i].J_z)); 
    #endif
  }
}
void omega2(scalar s, scalar tag ,int tg,int tagi,cg n[]){
  #if dimension == 2
  double omega_z[tg],J_z[tg],J_y[tg],J_x[tg],J_x_y[tg] ;
  for (int i = 0; i < tg; i++) omega_z[i]=J_z[i]=J_y[i]=J_x[i]=J_x_y[tg]=0;
  foreach(reduction(+:omega_z[:tg]) reduction(+:J_z[:tg]) reduction(+:J_x_y[:tg]) 
  reduction(+:J_y[:tg])  reduction(+:J_x[:tg]) ){
  #elif dimension == 3
  double omega_x[tg],omega_y[tg],omega_z[tg],J_z[tg],J_y[tg],J_x[tg],J_z_x[tg],J_y_z[tg],J_x_y[tg];
  for (int i = 0; i < tg; i++) omega_x[i]=omega_y[i]=omega_z[i]=J_z[i]=J_y[i]=J_x[i]=J_z_x[i]=J_y_z[i]=J_x_y[i]=0;
  foreach(reduction(+:omega_z[:tg]) reduction(+:J_z[:tg]) reduction(+:J_x_y[:tg])
          reduction(+:omega_y[:tg]) reduction(+:J_y[:tg]) reduction(+:J_y_z[:tg])
          reduction(+:omega_x[:tg]) reduction(+:J_x[:tg]) reduction(+:J_z_x[:tg])){
  #endif
    if (s[] != nodata && dv() > 0. && s[]>=EPS){
      int i=tag[]-1;
      coord pos = POS;
      pos = POS_perio(pos,n[tagi+i].per);
      coord posP = diff_pts(pos,n[tagi+i].pos);
      J_z[i] += dv()*s[]*(pow(posP.x,2)+pow(posP.y,2));
      omega_z[i] +=  dv()*s[]*(posP.x * u.y[] - posP.y * u.x[]); 
      J_x_y[i] -= dv()*s[]*posP.x*posP.y;
      #if dimension == 2
      J_x[i] += dv()*s[]*pow(posP.y,2);
      J_y[i] += dv()*s[]*pow(posP.x,2);
      #elif dimension == 3
      J_x[i] -= dv()*s[]*(pow(posP.y,2)+pow(posP.z,2));
      J_y[i] -= dv()*s[]*(pow(posP.z,2)+pow(posP.x,2));
      J_z_x[i] += dv()*s[]*posP.z*posP.x;
      J_y_z[i] += dv()*s[]*posP.y*posP.z;
      omega_x[i] +=  dv()*s[]*(posP.y * u.z[] - posP.z * u.y[]); 
      omega_y[i] +=  dv()*s[]*(posP.z * u.x[] - posP.x * u.z[]); 
      #endif
    }
  }
  for(int i=0;i<tg;i++){
    n[tagi+i].J_x = n[tagi+i].vol ? J_x[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_y = n[tagi+i].vol ? J_y[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_z = n[tagi+i].vol ? J_z[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_x_y = n[tagi+i].vol ? J_x_y[i]/n[tagi+i].vol : 0.;
    n[tagi+i].omega_z = J_z[i] ? omega_z[i]/J_z[i] : 0.;
    #if dimension == 2
    n[tagi+i].defnorm = sqrt(2*sq(n[tagi+i].J_x_y)+sq(n[tagi+i].J_x)+sq(n[tagi+i].J_y)+sq(n[tagi+i].J_z)); 
    #elif dimension == 3
    n[tagi+i].omega_x = J_x[i] ? omega_x[i]/J_x[i] : 0.;
    n[tagi+i].omega_y = J_y[i] ? omega_y[i]/J_y[i] : 0.;
    n[tagi+i].J_y_z = n[tagi+i].vol ? J_y_z[i]/n[tagi+i].vol : 0.;
    n[tagi+i].J_z_x = n[tagi+i].vol ? J_z_x[i]/n[tagi+i].vol : 0.;
    n[tagi+i].defnorm =sqrt(2*sq(n[tagi+i].J_z_x)+2*sq(n[tagi+i].J_y_z)+2*sq(n[tagi+i].J_x_y)+sq(n[tagi+i].J_x)+sq(n[tagi+i].J_y)+sq(n[tagi+i].J_z)); 
    #endif
  }
}


void dissipation_rate_on_bubbles (scalar s, scalar tag ,int tg,int tagi,cg n[])
{
  double diss[tg];
  for (int i = 0; i < tg; i++) diss[i]=0;
  foreach (reduction (+:diss[:tg])) {
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
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


void cg_bub (scalar s, scalar tag ,int tg,int tagi,cg n[])
{
  isperio(s,tag,tg,tagi,n);
  pos_vol_and_vel(s,tag,tg,tagi,n);
  omega(s,tag,tg,tagi,n);
  //interfacial NRJ ?? 
  // Pressure that is directly linked to the interfacial NRJ 
  dissipation_rate_on_bubbles(s,tag,tg,tagi,n);
  for(int i=0;i<tg;i++){
    foreach_dimension()
      if (n[tagi+i].pos.x>L0/2.) n[tagi+i].pos.x-=L0;
  }
}

void points_vel(int tg , cg n[]){
    for(int i=0;i<tg;i++){
      double a = n[i].pos.x ,b=n[i].pos.y;
      #if dimension == 3
        double c=n[i].pos.z;
        double vel_x,vel_y,vel_z,vel2_x,vel2_y,vel2_z;
        foreach_dimension()
          vel_x = interpolate(u.x,a,b,c) < Ls ?  interpolate(u.x,a,b,c) : 0; 
      #elif dimension == 2
        double vel2_x,vel2_y,vel_x,vel_y;
        foreach_dimension()
          vel_x = interpolate(u.x,a,b) < Ls ?  interpolate(u.x,a,b) : 0; 
      #endif
      #if _MPI
          foreach_dimension()  MPI_Allreduce(&vel_x, &vel2_x, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
      #else
        foreach_dimension() vel2_x = vel_x;
      #endif
      foreach_dimension() n[i].velc.x = vel2_x;
    }
}

void assign_tags_from_last_step (cg datat0[],cg datat1[],int tagi){
  for(int k = 0;k<tagi;k++) datat1[k].tag = -1;
  int tags_avail[tagi];
  for(int i=0;i<tagi;i++) tags_avail[i]=0;
  for (int j = 0; j < datat0[0].tagmax; j++)
  {
    coord np = datat0[j].pos;
    coord vp = datat0[j].vel;
    double dt = dtprint;
    np = add_pts(np,mult_pts(vp,dt)); // euleur schem
    #if dimension == 2
    double r = sqrt(datat0[j].vol/M_PI)*0.5; // 
    #else
    double r = datat0[j].vol/M_PI; // 
    #endif
    bool probleme_of_assignement = true;
    for (int k = j; k < (tagi+j); k++)
    {
      int idx = k % tagi;
      coord np1 = datat1[idx].pos;
      double dist = dist_perio(np,np1);
      if(dist<r){
        datat1[idx].tag = datat0[j].tag;
        probleme_of_assignement = false;
        break;
      }
    }
    if(probleme_of_assignement){
      fprintf(stdout,"Pbug : X  %g, Y  %g\n",np.x,np.y);
      fprintf(stdout,"Error of assignement at t = %g for bubble tag = %d tm %d j %d \n",t,datat0[j].tag,datat0[0].tagmax,j);
    }
  }
  
  // if coalescence or breakup
  int ta = 0;
  for (int k = 0; k < tagi; k++){
    bool avail = true;
    for (int j = 0; j < tagi; j++)
      if (k == datat1[j].tag)
        avail =false;
    if(avail) {
      tags_avail[ta] = k;
      ta++;
    }
    if(datat1[k].tag>=tagi) datat1[k].tag=-1;
  }
  // fprintf(stdout,"Tagavail (");
  // for (int i = 0; i < ta; i++)fprintf(stdout,"%d,",tags_avail[i]);
  // fprintf(stdout,")\n");
  // fprintf(stdout,"current tags (");
  // for (int i = 0; i < tagi; i++)fprintf(stdout,"%d,",datat1[i].tag);
  // fprintf(stdout,")\n");

  int new_tag=0;
  for (int k = 0; k < tagi; k++){
    if(datat1[k].tag == -1 && new_tag<=ta){
      if(tags_avail[new_tag] >= tagi) new_tag++;
      datat1[k].tag = tags_avail[new_tag];
      new_tag++;
    }
  }
}

void print_pair_dist(cg data[],cg data0[],int tagi,int step, float t){
  cg data_orderred[tagi];
  cg data_orderred0[data0[0].tagmax];
  for (int j = 0; j < tagi; j++)
  {
    data[j].ct = 0;
    foreach_dimension()
      data[j].distmin.x = Ls*100;
    data_orderred[data[j].tag]=data[j];
  }


  foreach_dimension(){
  Dist_x = fopen(name_x,"a");
  fprintf(Dist_x,"%d,%g",step,t);
  }

  for (int k = 0; k < tagi; k++)
  {
    coord npi = data_orderred[k].pos; 
    for (int j = k+1; j<tagi;j++){
      coord npj = data_orderred[j].pos;
      coord distPts = dist_perioPts(npi,npj);
      double dist = normL2_pts(distPts);
      foreach_dimension() fprintf(Dist_x,",%g",distPts.x);
      if(dist<normL2_pts(data_orderred[k].distmin)){
        data_orderred[k].distmin = distPts;
        data_orderred[k].tagmin = j;
      } 
      if(dist < normL2_pts(data_orderred[j].distmin)){
        data_orderred[j].distmin = mult_pts(distPts,-1);
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
  foreach_dimension(){
  fprintf(Dist_x,"\n");
  fclose(Dist_x);
  }

  for (int j = 0; j < data0[0].tagmax; j++)
    data_orderred0[data0[j].tag]=data0[j];
  for(int k = 0; k<tagi;k++){
    int tag = data[k].tag;
    int tm = data[k].tagmin;
    double Dia  = sqrt(data_orderred[tm].vol/ M_PI )+sqrt(data[k].vol/ M_PI );

    if((normL2_pts(data[k].distmin) <= Dia) && (data[k].adj[data_orderred[tm].si] == true)){
      if(tag >= data0[0].tagmax) data[k].ct = dtprint;
      if(tag <  data0[0].tagmax) data[k].ct = data_orderred0[tag].ct + dtprint;
    }
  }

}

void time_of_contact(scalar s, scalar tag ,int tg,int tagi,cg n[])
{
  int taggg[tg];
  int nvar = datasize/sizeof(double), len = tg*nvar, adj[len],tc[len];
  for (int i = 0; i < len; i++) adj[i] = false;
  for (int i = 0; i < tg; i++) taggg[i]=tc[i]=0;
  foreach(reduction(max:taggg[:tg])  reduction(max:tc[:tg]) reduction(max:adj[:nvar]))
    if (s[] != nodata && dv() > 0.  && s[]>=EPS) {
      int i=tag[]-1;
      taggg[i] = tag[];
      foreach_neighbor()
      for(scalar c in interfaces)
        if(c.i != s.i && c[] > EPS)
          adj[i*c.i + c.i] = true;
  }
  for (int i = 0; i < tg; i++){
    n[tagi+i].realtag = taggg[i];
    for(scalar c in interfaces){
     n[tagi+i].adj[c.i] = adj[i*c.i + c.i];
    }
  } 
}


void coalescence(cg data[],int tagi){
  cg DAT[tagi];
  for (int j = 0; j < tagi; j++) DAT[data[j].tag]=data[j];
  scalar b[];

  for(scalar c in interfaces){
    foreach(){
      b[] = c[]>EPS;}
    tag (b);

    int si[tagi],realtag[tagi]; 
    for (int i = 0; i < tagi; i++) si[i]=realtag[i]=-1;
    for(int j = 0 ;j < tagi ; j++){
      if(DAT[j].ct != 0.){
        fprintf(stdout,"time %g the tag %d %d ct %g\n",t,j,DAT[j].tag,DAT[j].ct);
        cg b1 = DAT[j];
        cg b2 = DAT[DAT[j].tagmin];
        double d_eq = 2*sqrt(b2.vol/M_PI)*sqrt(b1.vol/M_PI) / (sqrt(b2.vol/M_PI) + sqrt(b1.vol/M_PI));
        double V_0 = normL2_pts(diff_pts(b1.vel,b2.vel)); 
        double td = rho_f*sq(d_eq)*V_0/(8*sig);
        double We = rho_f*d_eq*sq(V_0)/(2*sig);
        fprintf(stdout,"coalescence time td = %g We = %g deq %g vel %g\n",td,We,d_eq,V_0);
        if (DAT[j].ct > td && DAT[j].si  != DAT[DAT[j].tagmin].si && c.i == DAT[j].si) {
          // for(int l = 0 ;l < tagi ; l++){
            // int tm = DAT[j].tag;
            // if(DAT[l].tag != DAT[j].tag && DAT[l].tag != DAT[tm].tag){
            //   if(DAT[l].adj[DAT[j].si] == 1 && DAT[l].si == DAT[tm].si)
            //   if(DAT[l].tagmin == DAT[j].tag){
            //   fprintf(stdout,"l too close too [");
            //   for(scalar s in interfaces) fprintf(stdout,"(%d,%d),",DAT[l].adj[s.i],s.i);
            //   fprintf(stdout,"]\n");
            //   }
            // }else{
              realtag[j] = DAT[j].realtag;
              si[j] = DAT[DAT[j].tagmin].si;
              DAT[j].si  = DAT[DAT[j].tagmin].si;
              fprintf(stdout,"add to caaol realt %d si objc %d\n",si[j],realtag[j]);
          //   }
          // }
        }
      }
    }
    #if _MPI
        MPI_Bcast(si,tagi,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(realtag,tagi,MPI_INT,0,MPI_COMM_WORLD);
    #endif


    for (int i = 0; i < tagi; i++)
      if(si[i] != -1)
      foreach()
        if(b[] == realtag[i]){
        scalar t = {si[i]};
        t[] = c[]; 
        c[] = 0.;
        }
    scalar * list = list_copy ({c});
    for (int i = 0; i < tagi; i++)
      if(si[i] != -1)
        list = list_add(list,(scalar){si[i]});
    boundary (list);
    free (list);
  }
}
