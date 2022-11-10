// definition of point structure and simple algebra operators
int main_process(){
  #if _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank==0;
  #else
  return 1;
  #endif
}

typedef struct ten{
  coord x;
  coord y;
  #if dimension == 3
  coord z;
  #endif
}ten;

coord add_pts(coord a,coord b){
  #if dimension == 2
  coord c = {a.x+b.x, a.y+b.y};
  #elif dimension == 3
  coord c = {a.x+b.x, a.y+b.y, a.z+b.z};
  #endif
  return c;
}
coord diff_pts(coord a,coord b){
  #if dimension == 2
  coord c = {a.x-b.x,a.y-b.y};
  #elif dimension == 3
  coord c = {a.x-b.x, a.y-b.y, a.z-b.z};
  #endif
  return c;
}
coord mult_pts(coord a,double b){
  #if dimension == 2
  coord c = {a.x*b,  a.y*b};
  #elif dimension == 3
  coord c = {a.x*b,  a.y*b,  a.z*b};
  #endif
  return c;
}
coord div_pts(coord a,double b){
  #if dimension == 2
  coord c = {a.x/b,  a.y/b};
  #elif dimension == 3
  coord c = {a.x/b,  a.y/b,  a.z/b};
  #endif
  return c;
}
double normL2_pts(coord a){
  #if dimension == 2
  double c = sqrt(pow(a.x,2)+pow(a.y,2));
  #elif dimension == 3
  double c = sqrt(pow(a.x,2)+pow(a.y,2)+pow(a.z,2));
  #endif
  return c;
}
double normL2_pts_sq(coord a){
  #if dimension == 2
  double c = sq(a.x)+sq(a.y);
  #elif dimension == 3
  double c = sq(a.x)+sq(a.y)+sq(a.z);
  #endif
  return c;
}
double normL2_ten(ten a){
  #if dimension == 2
  double c = sqrt(pow(a.x.x,2)+pow(a.y.y,2)+pow(a.x.y,2)+pow(a.y.x,2));
  #elif dimension == 3
  double c = sqrt(pow(a.x.x,2)+pow(a.y.y,2)+pow(a.z.z,2)+pow(a.x.y,2)+pow(a.y.x,2)+pow(a.x.z,2)+pow(a.z.x,2)+pow(a.z.y,2)+pow(a.y.z,2));
  #endif
  return c;
}
//only in 3D for the cross product 
#if dimension == 3
coord cross_pts(coord a,coord b){
  coord c = {a.y*b.z-a.z*b.y, 
           a.z*b.x-a.x*b.z, 
           a.x*b.y-a.y*b.x};
#elif dimension == 2
double cross_pts(coord a,coord b){
  double c = a.x*b.y-a.y*b.x;
#endif
  return c;
}

coord EigenValue(ten A){
  double Delta = pow(A.x.x,2.) - 2.*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x;
  coord Eig = {0.,0.};
  if(Delta > 0){
    Eig.x = 1./2.*(  A.x.x + A.y.y - sqrt( Delta ));
    Eig.y = 1./2.*(  A.x.x + A.y.y + sqrt( Delta ));
  }
  return Eig;
}

ten EigenVectors(ten A){
  ten V;
  if(A.y.x != 0){
    V.x.x = -1./(2.*A.y.x)*( - A.x.x + A.y.y + sqrt( pow(A.x.x,2.) - 2.*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x ) ); 
    V.x.y = 1.;
    V.y.x = -1./(2.*A.y.x)*( - A.x.x + A.y.y - sqrt( pow(A.x.x,2.) - 2*A.x.x*A.y.y+pow(A.y.y,2.) +4.*A.x.y*A.y.x ) );
    V.y.y = 1.;
  }else if((A.x.x-A.y.y)!=0){
    V.x.x = 1.;
    V.x.y = 0.;
    V.y.x = -A.x.y/(A.x.x-A.y.y);
    V.y.y = 1.;
  }else{
    V.x.x = 1.;
    V.x.y = 0.;
    V.y.x = 0.;
    V.y.y = 1.;
  }
  V.x = div_pts(V.x,normL2_pts(V.x));
  V.y = div_pts(V.y,normL2_pts(V.y));
  return V;
}


