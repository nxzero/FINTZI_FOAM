/**
# Simulation of a rising suspension of drops

*/
/** Necessary header files. no-coalescence.h can be found [here](http://basilisk.fr/sandbox/popinet/no-coalescence.h)*/

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "no-coalescence.h"
#include "view.h"
#include "output.h"
// #include "navier-stokes/conserving.h"
#include "parameters.h"
#include "algebra.h"
#include "cgposition_opt.h"
#include "bubbles.h"
// FILE * Pos;
FILE * Pos;
FILE * Dist;
FILE * Stats;
FILE * PT_tensor;
FILE * Para;

/** Function used to initialise the drops on a square grid so the maximum possible volume fraction $\phi=\frac{\pi}{4}\approx0.785$*/
#if dimension == 2
double geometry( double x, double y, double z){
  double inout=-1.,circ=0.;
  srand(3458829);
  for(double i=(-Ls/2.+DP/2.); i<= Ls/2; i+=DP){
    for( double j=(-Ls/2.+DP/2.); j<= Ls/2.; j+=DP){
      int A=rand();
      int B=rand();
      circ= -sq(x-i+RND*(A*2./RAND_MAX-1))-sq(y-j+RND*(B*2./RAND_MAX-1))+sq(D/2);
      if (circ>=0){
	      inout=circ;
      }
    }
  }
  return(inout);
}
#elif dimension == 3
double geometry( double x, double y, double z){
  double inout=-1.,circ=0.;
  srand(3458829);
  for(double i=(-Ls/2.+DP/2.); i<= Ls/2; i+=DP){
    for( double j=(-Ls/2.+DP/2.); j<= Ls/2.; j+=DP){
      for( double k=(-Ls/2.+DP/2.); k<= Ls/2.; k+=DP){
        int A=rand();
        int B=rand();
        int C=rand();
        circ= -sq(x-i+RND*(A*2./RAND_MAX-1))-sq(y-j+RND*(B*2./RAND_MAX-1))-sq(z-k+RND*(C*2./RAND_MAX-1))+sq(D/2);
        if (circ>=0){
	        inout=circ;
        }
      }
    }
  }
  return(inout);
}
#endif

/** Event to recover the global parameters of each drop works by tagging each drop and then outputs to stderr the iteration, time, position of the CG, velocity, and tag value for the drop*/
scalar tag_level[],tg_l[];
cg BubblesPos[Nb+2]; // not sure the +2 is really usefulll but in case there is fragment we need more space  
/** tmp variable to store the previous step in order to track bubbles. 
 * the hole process add less than 1% simulation time */ 
cg BubblesPosTmp[Nb+2]; 

event track_bub(t=FIRSTSTEP; t<= TMAX;t+=dtprint){
  foreach()
    tg_l[]=0;
  int tagi = 0;
  int nbubble;
  for(scalar s in interfaces){
    foreach()
      tag_level[]=s[]>1e-4;
    nbubble=tag(tag_level);
    cg_bub(s,tag_level,nbubble,tagi,BubblesPos);
    for(int j=0;j<nbubble;j++){
      if(BubblesPos[tagi+j].vol<1e-2){//Mettre relativement a la taille d'une maille 
        foreach(){ //Fragment suppression
          if(tag_level[]==j+1){
            s[]=0;
          }
        }
        // remove the fragment from the list by replacing the other 
        for (int k = j; k < nbubble-1; k++) { 
          BubblesPos[tagi+k] = BubblesPos[tagi+k + 1];
        }
        nbubble -= 1;
        j -=1;
      }
    }
    foreach(){
      tag_level[]*=s[];
      tg_l[]+=tag_level[];
    }
    for(int k=0;k<nbubble;k++)
    {
      BubblesPos[tagi].tag = tagi;
      tagi++;
    }
  }
  if(t==FIRSTSTEP){
    memcpy(BubblesPosTmp, BubblesPos, sizeof(cg)*tagi);// copy the arry
  }

  assign_tags_from_last_step(BubblesPosTmp,BubblesPos,tagi);
  memcpy(BubblesPosTmp, BubblesPos, sizeof(cg)*tagi);// copy the arry

  // print the datas with the main process only
  if(main_process()){
    print_pair_dist(BubblesPos,tagi,Dist,i,t);
    Pos = fopen("Pos.csv","a");
    for(int j=0;j<tagi;j++){
      #if dimension == 2
  	  fprintf(Pos,"%d,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d\n",i,t,tagi,BubblesPos[j].pos.x,BubblesPos[j].pos.y,BubblesPos[j].vel.x,BubblesPos[j].vel.y,BubblesPos[j].omegaz,BubblesPos[j].vol,BubblesPos[j].EigVec.x.x,BubblesPos[j].EigVec.x.y,BubblesPos[j].defnorm,BubblesPos[j].diss,BubblesPos[j].distmin,BubblesPos[j].tagmin,BubblesPos[j].tag);
      #elif dimension == 3
  	  fprintf(Pos,"%d,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d\n",i,t,tagi,BubblesPos[j].pos.x,BubblesPos[j].pos.y,BubblesPos[j].pos.z,BubblesPos[j].vel.x,BubblesPos[j].vel.y,BubblesPos[j].vel.z,BubblesPos[j].omega.x,BubblesPos[j].omega.y,BubblesPos[j].omega.z,BubblesPos[j].vol,BubblesPos[j].defnorm,BubblesPos[j].diss,BubblesPos[j].distmin,BubblesPos[j].tagmin,BubblesPos[j].tag);
      #endif
    }
    fclose(Pos);
  }
}


event average_stats(t+=dtprint){
  sta st = calcul_sta();
  if(main_process()){
    Stats = fopen("Stats.csv","a");
    #if dimension == 2
    fprintf(Stats,"%d,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",i,t,st.nvof,st.veld.x,st.veld.y,st.velf.x,st.velf.y,st.rhov.x,st.rhov.y,st.vel.x,st.vel.y,st.Up_f.x,st.Up_f.y,st.dissd,st.dissf);
    #elif dimension == 3
    fprintf(Stats,"%d,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",i,t,st.nvof,st.veld.x,st.veld.y,st.veld.z,st.velf.x,st.velf.y,st.velf.z,st.rhov.x,st.rhov.y,st.rhov.z,st.vel.x,st.vel.y,st.vel.z,st.Up_f.x,st.Up_f.y,st.Up_f.z,st.dissd,st.dissf);
    #endif
    fclose(Stats);
  }
}


/** The flow parameters are calculated based on the previously defined non-dimensional parameters. The tolerance is reduced to 1e-4 and the boundaries are set to be periodic*/
int main(){
  size(Ls);
  init_grid(1<<LEVEL);
  origin(-Ls/2.,-Ls/2.,-Ls/2.);
  rho1=r1;
  mu1=mu_d;//sqrt(fabs(1-rho_r)*rho_r*g*(D*D*D))*rho1/(Ga*mu_r);
  rho2=rho_f;//rho1*rho_r;
  mu2=mu_f;//mu1*mu_r;
  f.sigma=sig;//fabs(1-rho_r)*rho1*g*D*D/Bo;
  TOLERANCE=1e-4;
  foreach_dimension()
    periodic(right);
  run();
}
/**Initializes the domain with zero velocity and outputs the flow parameters*/
event init (t=0){
  fraction(f, geometry(x,y,z));
  foreach()
    u.x[]=u.y[]=u.z[]=0.;
  if(main_process()){
    // files where the particules stats are stored 
    Pos = fopen("Pos.csv","w");
    #if dimension ==2
    fprintf(Pos,"i,t,tagi,x,y,vx,vy,omegaz,vol,VIx,VIy,defnorm,diss,distmin,tagmin,tag\n");
    #elif dimension ==3
    fprintf(Pos,"i,t,tagi,x,y,z,vx,vy,vz,omegax,omegay,omegaz,vol,defnorm,diss,distmin,tagmin,tag\n");
    #endif
    fclose(Pos);
    // files where the dist between all particules are stored 
    Dist = fopen("Dist.csv","w");
    fprintf(Dist,"i,t");
    for (int j = 0; j < Nb; j++)
    {
      for (int k=j+1; k < Nb; k++){
        fprintf(Dist,",b%db%d",j,k);
      }
    }
    // files where the fluid and solid phase stats are stored 
    fprintf(Dist,"\n");
    fclose(Dist);
    Stats = fopen("Stats.csv","w");
    #if dimension == 2
    fprintf(Stats,"i,t,nvof,vxdrops,vydrops,vxfluid,vyfluid,rhovx,rhovy,Vx,Vy,xp,yp,dissD,dissF\n");
    #elif dimension == 3
    fprintf(Stats,"i,t,nvof,vxdrops,vydrops,vzdrops,vxfluid,vyfluid,vzfluid,rhovx,rhovy,rhovz,Vx,Vy,Vz,xp,yp,zp,dissF,dissD\n");
    #endif
    fclose(Stats);
  }
}


/** Overloads the acceleration event to apply the effect of gravity. Due to the periodic boundary conditions the acceleration needs to be reduced by $\frac{\rho_{av}}{\rho}g$ rho is not available directly as a face centered vector so the stagger of f[] is adjusted to calculate it */
event acceleration (i++) {
  double rhoav=avg_rho();
  face vector av = a;
  foreach_face(y)
    av.y[] -= (1-rhoav/((f[]+f[-1])/2*(rho1-rho2)+rho2))*g;
}
/** Outputs videos of the velocity, Pressure, vorticity, and tag fields*/
// event movies(t=FIRSTSTEP; t<=TMAX; t+=MOVIES)
event movies(i++)
{
  #if dimension == 3
    view (fov = 30.86528,
	quat = {0.515965,0.140691,0.245247,0.808605},
	tx = -0.07438, ty = -0.0612925,
	width = 1024, height = 768);
  #endif

  box();
  draw_vof("f", lw = 2.);
  squares("u.y");
  save("uy.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("u.x");
  save("ux.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares("p");
  save("P.mp4");

  clear();

  scalar omega[];
  vorticity (u, omega);

  box();
  draw_vof("f", lw = 2.);
  squares("omega");
  save("vort.mp4");

  clear();

  box();
  draw_vof("f", lw = 2.);
  squares(color="tg_l",map=randomap,min=0,max=25);
  save("tag.mp4");
  clear();
  dump (file = "dump-last");

}

event SaveDumps(t+=SAVESTEP){
  char named[80];
  sprintf (named, "dump-%d", i);
  dump (file = named);
}
event stop(t=TMAX){
  return 1;
}
