#include "dca.h"

int test() {
  double rec_x[4], rec_y[4], rec_z[4];
  
  rec_x[0]=1;rec_y[0]=0;rec_z[0]=0;
  rec_x[1]=0;rec_y[1]=1;rec_z[1]=0;
  rec_x[2]=0;rec_y[2]=0;rec_z[2]=2;
  rec_x[3]=0;rec_y[3]=2;rec_z[3]=0;

  double dca_x, dca_y, dca_z, dca;

  struct Line l1;
  l1.p1.x=rec_x[0];
  l1.p1.y=rec_y[0];
  l1.p1.z=rec_z[0];
  l1.p2.x=rec_x[1];
  l1.p2.y=rec_y[1];
  l1.p2.z=rec_z[1];
  printf("l1: (%f %f %f)\n",l1.p2.x-l1.p1.x,l1.p2.y-l1.p1.y,l1.p2.z-l1.p1.z);
  struct Line l2;
  l2.p1.x=rec_x[2];
  l2.p1.y=rec_y[2];
  l2.p1.z=rec_z[2];
  l2.p2.x=rec_x[3];
  l2.p2.y=rec_y[3];
  l2.p2.z=rec_z[3];
  printf("l2: (%f %f %f)\n",l2.p2.x-l2.p1.x,l2.p2.y-l2.p1.y,l2.p2.z-l2.p1.z);

  struct Line lper = line_perpendicular(l1,l2);
  printf("dca=%f p(%f %f %f)\n",pub_dca,pub_dca_x,pub_dca_y,pub_dca_z);
  cout<<"status = "<<lper.status<<endl;

  double dca_x2, dca_y2, dca_z2, dca2;
  struct point PoCA = cal_PoCA(l1,l2);
  dca_x2 = PoCA.x;
  dca_y2 = PoCA.y;
  dca_z2 = PoCA.z;
  dca2 = cal_DCA(l1,l2);
  printf("dca=%f p(%f %f %f)\n",dca2,dca_x2,dca_y2,dca_z2);

    return 0;
}

