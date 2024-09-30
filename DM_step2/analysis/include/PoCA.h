#ifndef POCA_h
#define POCA_h
using namespace std;

#include <iostream>
#include <cmath>
#include <vector>
#include "vector.h"
#include "draw.h"

V3 GetPoCAPoint(V3 const& p1, V3 const& p2,
                V3 const& p3, V3 const& p4);

double GetDCA(V3 const& p1, V3 const& p2,
              V3 const& p3, V3 const& p4);

bool CheckPoCAStatus(V3 const& p1, V3 const& p2,
                     V3 const& p3, V3 const& p4);

bool CollectPath(V3 const& p1, V3 const& p2,
                 V3 const& p3, V3 const& p4,
                 double E,
                 std::vector< int >& count,
                 std::vector< double >& sig);

void median_filter(const vector<double> &lambda,
                   vector<double> &next_lambda,
                   const int size = 1);

double pub_dca_x, pub_dca_y, pub_dca_z, pub_dca;
const double delta=1e-18;
struct point {
  double x;
  double y;
  double z;
  bool status;
};

struct Line {
  struct point p1;
  struct point p2;
  bool status;
};

double cal_ang(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4);
struct Line line_perpendicular(struct Line l1, struct Line l2);
struct point cal_normal_vector(struct Line l1, struct Line l2);
struct point normalization(struct point p1);

V3 GetPoCAPoint(V3 const& p1, V3 const& p2,
                V3 const& p3, V3 const& p4) {
    V3 v_in = p2 - p1; // 入射矢量，从 p1 指向 p2
    V3 v_out = p4 - p3; // 出射矢量，从 p3 指向 p4
    V3 v_n = v_in.cross(v_out);
    v_n = v_n.normalize(); // 法向量
    double d = (p3 - p2).dot(v_n); // 距离
    double t_i = (v_out.x * (d * v_n.y + p2.y - p3.y) - v_out.y * (d * v_n.x + p2.x - p3.x)) / (v_out.x * v_in.y - v_in.x * v_out.y);
    return p2 - (t_i * v_in) + (0.5 * d * v_n);
}

double GetDCA(V3 const& p1, V3 const& p2,
              V3 const& p3, V3 const& p4) {
    V3 v_in = p2 - p1; // 入射矢量，从 p1 指向 p2
    V3 v_out = p4 - p3; // 出射矢量，从 p3 指向 p4
    V3 v_n = v_in.cross(v_out);
    v_n = v_n.normalize(); // 法向量
    double d = (p3 - p2).dot(v_n); // 距离
    return fabs(d);
}

bool CheckPoCAStatus(V3 const& p1, V3 const& p2,
                     V3 const& p3, V3 const& p4) {
    V3 v_in = p2 - p1; // 入射矢量，从 p1 指向 p2
    V3 v_out = p4 - p3; // 出射矢量，从 p3 指向 p4
    V3 v_n = v_in.cross(v_out);
    double length = v_n.length();
    if(isnan(length)) return false;
    if(length < delta) return false;
    else return true;
}

extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;

// 在 ../utils/draw.h 中定义 Nx、Ny、Nz、Lx、Ly、Lz
// 在 ../utils/draw.cpp 中定义 Z1、Z2、Z3、Z4

bool CollectPath(V3 const &p1, V3 const &p2,
                 V3 const &p3, V3 const &p4,
                 double E,
                 vector<int> &count, // M_j
                 vector<double> &sig // k
)
{
    const double dx = Lx / Nx, // 2mm
                 dy = Ly / Ny, // 2mm
                 dz = Lz / Nz; // 5mm
                 //dz = (Z3 - Z2) / Nz; // 6mm

    // 入射矢量、出射矢量、PoCA 点
    V3 v_in = p2 - p1;                   // p2 = P_in
    V3 v_out = p4 - p3;                  // p3 = P_out
    V3 m = GetPoCAPoint(p1, p2, p3, p4); // m = PoCA

    // PoCA 点所在的格子
    int x_index = (m.x + Lx / 2) / dx,
        y_index = (m.y + Ly / 2) / dy,
        z_index = (m.z + Lz / 2) / dz;
        //z_index = (m.z - Z2) / dz;

    if (x_index >= Nx || x_index < 0 ||
        y_index >= Ny || y_index < 0 ||
        z_index >= Nz || z_index < 0)
        return false;

    // lambda_j = 1/M_j \sum_i{1/L_ij (p_i/p_0)^2 (\Delta\theta_{x,ij}^2 + \Delta\theta_{y,ij}^2)/2}

    double d_theta_x = atan2(v_out.x, v_out.z) - atan2(v_in.x, v_in.z);
    double d_theta_y = atan2(v_out.y, v_out.z) - atan2(v_in.y, v_in.z);

    // pr_pow2 = (p/p_0)^2 = (p^2 c^2)/(p_0^2 c^2) = (E^2 - (m_0 c^2)^2)/(p_0^2 c^2)
    // m_0 c^2 = 105.66 MeV
    double pr_pow2 = (E * E - 105.65836668 * 105.65836668) / 9000000;

    // 粒子路径在 PoCA 点所在的格子里的长度 L_ij
    double path_length =
        (p3 - m).length() * ((z_index + 1) * dz + Z2 - m.z) / (Z3 - m.z) +
        (m - p2).length() * (m.z - Z2 - z_index * dz) / (m.z - Z2);
    path_length /= 1000; // 毫米换成米

    // S = 1/L_ij (p_i/p_0)^2 (\Delta\theta_{x,ij}^2 + \Delta\theta_{y,ij}^2)/2
    double S = pr_pow2 *
               (d_theta_x * d_theta_x + d_theta_y * d_theta_y) / 2 / path_length;

    // lambda_j = 1/M_j \sum_i{S}
    auto index = INDEX(x_index, y_index, z_index);
    count[index] += 1;
    // sig[index] += S; // PoCA 点处 \lambda = \lambda + \Delta\lambda

    double cos_theta = (m - p2).dot(p3 - m) / ((m - p2).length() * (p3 - m).length());
    double sigma = 100 * acos(abs(cos_theta));
    double delta_1 = 0.5 * erf(1 / (sigma * sqrt(2.)));
    double delta_2 = 0.5 * erf(2 / (sigma * sqrt(2.))) - delta_1;
    double delta_3 = 0.5 - delta_1 - delta_2;

    if (z_index > 1 && z_index < (Nz - 2))
    {
        sig[index] += 2 * delta_1 * S;
        sig[INDEX(x_index, y_index, z_index - 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index + 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index - 2)] += delta_3 * S;
        sig[INDEX(x_index, y_index, z_index + 2)] += delta_3 * S;
    }
    else if (z_index == 1)
    {
        sig[index] += 2 * delta_1 * S;
        sig[INDEX(x_index, y_index, z_index - 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index + 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index + 2)] += delta_3 * S;
    }
    else if (z_index == 0)
    {
        sig[index] += 2 * delta_1 * S;
        sig[INDEX(x_index, y_index, z_index + 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index + 2)] += delta_3 * S;
    }
    else if (z_index == Nz - 2)
    {
        sig[index] += 2 * delta_1 * S;
        sig[INDEX(x_index, y_index, z_index - 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index + 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index - 2)] += delta_3 * S;
    }
    else if (z_index == Nz - 1)
    {
        sig[index] += 2 * delta_1 * S;
        sig[INDEX(x_index, y_index, z_index - 1)] += delta_2 * S;
        sig[INDEX(x_index, y_index, z_index - 2)] += delta_3 * S;
    }

    // 粒子经过哪些格子
    for (int i = 0; i < Nz; ++i)
    {
        //double k, x, y, z = (i + 0.5) * dz + Z2;
        double k, x, y, z = (i + 0.5) * dz - Lz / 2;
        if (i < z_index)
        { // 粒子的入射轨迹：m - p2
            k = (z - Z2) / (m.z - Z2);
            x = k * (m.x - p2.x) + p2.x;
            y = k * (m.y - p2.y) + p2.y;
        }
        else if (i > z_index)
        { // 粒子的出射轨迹：p3 - m
            k = (z - m.z) / (Z3 - m.z);
            x = k * (p3.x - m.x) + m.x;
            y = k * (p3.y - m.y) + m.y;
        }
        else
            continue; // continue 会跳过当前循环中的代码，强迫开始下一次循环

        x_index = (x + Lx / 2) / dx;
        y_index = (y + Ly / 2) / dy;

        // out of box
        if (x_index >= Nx || x_index < 0 ||
            y_index >= Ny || y_index < 0)
            continue;

        count[INDEX(x_index, y_index, i)] += 1;
    }

    return true;
}

// 中值滤波算法
void median_filter(const vector<double> &lambda,
                   vector<double> &next_lambda,
                   const int size = 1)
{
    double values[1000];
    for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
            {
                int n = 0;
                for (int k_ = k - size; k_ <= k + size; ++k_)
                    for (int j_ = j - size; j_ <= j + size; ++j_)
                        for (int i_ = i - size; i_ <= i + size; ++i_)
                            if (i_ >= 0 && i_ < Nx &&
                                j_ >= 0 && j_ < Ny &&
                                k_ >= 0 && k_ < Nz)
                                values[n++] = lambda[INDEX(i_, j_, k_)];
                sort(values, values + n);
                int delnum = 3; // numbers to discard, from two sides
                if (n - 2 * delnum > 0)
                {
                    double sum = 0;
                    for (int index = delnum; index < n - delnum; index++)
                    {
                        sum = sum + values[index];
                    }
                    next_lambda[INDEX(i, j, k)] = sum / (n - 2 * delnum);
                }
                else
                {
                    // if we discard too many points, use middle
                    next_lambda[INDEX(i, j, k)] = (values[n / 2] +
                                                   values[(n - 1) / 2]) /
                                                  2;
                }
            }
}

//old functions
double cal_ang(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4){
  // 计算两条直线的向量
  double v1x = x2 - x1;
  double v1y = y2 - y1;
  double v1z = z2 - z1;
  double v2x = x4 - x3;
  double v2y = y4 - y3;
  double v2z = z4 - z3;

  // 计算两条直线的点积和向量长度
  double dot_product = v1x * v2x + v1y * v2y + v1z * v2z;
  double length1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
  double length2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

  // 计算夹角
  double angle = acos(dot_product / (length1 * length2)) * 180 / M_PI;
  return angle;
};

struct Line line_perpendicular(struct Line l1, struct Line l2) {
  //printf("\nYuxudong =>\n");
  double x1 = l1.p1.x;
  double y1 = l1.p1.y;
  double z1 = l1.p1.z;
  double x2 = l1.p2.x;
  double y2 = l1.p2.y;
  double z2 = l1.p2.z;
  double x3 = l2.p1.x;
  double y3 = l2.p1.y;
  double z3 = l2.p1.z;
  double x4 = l2.p2.x;
  double y4 = l2.p2.y;
  double z4 = l2.p2.z;

  // step 1: the normal vector of l1 and l2
  // n12 = l1 x l2
  struct point n12_vec = cal_normal_vector(l1, l2);
  //cout<<"n12_vec = ("<<n12_vec.x<<","<<n12_vec.y<<","<<n12_vec.z<<")"<<endl;
  struct point n12_unit = normalization(n12_vec);
  //cout<<"n12_unit = ("<<n12_unit.x<<","<<n12_unit.y<<","<<n12_unit.z<<")"<<endl;
  if(n12_unit.status==false)
  {
    struct point p1{0,0,0};
    struct point p2{1,1,1};
    struct Line lper;
    lper.p1 = p1;
    lper.p2 = p2;
    lper.status=false;
    return lper;
  };

  // the direction vector of l1 and l2
  struct point l1_vec;
  l1_vec.x = x2-x1;
  l1_vec.y = y2-y1;
  l1_vec.z = z2-z1;
  struct point l1_unit = normalization(l1_vec);
  //cout<<"l1_unit = ("<<l1_unit.x<<","<<l1_unit.y<<","<<l1_unit.z<<")"<<endl;
  struct point l2_vec;
  l2_vec.x = x4-x3;
  l2_vec.y = y4-y3;
  l2_vec.z = z4-z3;
  //cout<<"l2_vec = ("<<l2_vec.x<<","<<l2_vec.y<<","<<l2_vec.z<<")"<<endl;
  struct point l2_unit = normalization(l2_vec);
  //cout<<"l2_unit = ("<<l2_unit.x<<","<<l2_unit.y<<","<<l2_unit.z<<")"<<endl;

  // the normal vector of l1 and n12
  struct Line l12 = {{0,0,0}, {n12_unit.x,n12_unit.y,n12_unit.z}};
  struct point n12prime_vec = cal_normal_vector(l1, l12);
  struct point n12prime_unit = normalization(n12prime_vec);
  //cout<<"n12prime_unit = ("<<n12prime_unit.x<<","<<n12prime_unit.y<<","<<n12prime_unit.z<<")"<<endl;

  // step 2: find point in l2 in the plane (n12 x l1)
  // define intersection point is p1
  // xp1 = x3 + t1*l2_unit.x
  // yp1 = y3 + t1*l2_unit.y
  // zp1 = z3 + t1*l2_unit.z
  // line (p1,x1) should be perpendicular to n12prime_unit
  // line (p1,x1) dot n12prime_unit = 0

  double t1 = 0;
  //if(fabs((l2_unit.x*n12prime_unit.x+l2_unit.y*n12prime_unit.y+l2_unit.z*n12prime_unit.z))>delta) t1=-(x3*n12prime_unit.x+y3*n12prime_unit.y+z3*n12prime_unit.z) / (l2_unit.x*n12prime_unit.x+l2_unit.y*n12prime_unit.y+l2_unit.z*n12prime_unit.z);
  if(fabs((l2_unit.x*n12prime_unit.x+l2_unit.y*n12prime_unit.y+l2_unit.z*n12prime_unit.z))>delta) t1=-(x3*n12prime_unit.x+y3*n12prime_unit.y+z3*n12prime_unit.z - x1*n12prime_unit.x-y1*n12prime_unit.y-z1*n12prime_unit.z) / (l2_unit.x*n12prime_unit.x+l2_unit.y*n12prime_unit.y+l2_unit.z*n12prime_unit.z);
  //cout<<"t1 = "<<t1<<endl;
  struct point p1;
  p1.x = x3+t1*l2_unit.x;
  p1.y = y3+t1*l2_unit.y;
  p1.z = z3+t1*l2_unit.z;
  //cout<<"p1 = ("<<p1.x<<","<<p1.y<<","<<p1.z<<")"<<endl;

  // step 3: find point p2 in l1
  // line (p2, p1) is parrallel to n12_unit
  // xp2 = xp1 + t2*n12_unit.x
  // yp2 = yp1 + t2*n12_unit.y
  // zp2 = zp1 + t2*n12_unit.z
  // line (p2,x1) should parrallel to l1
  // note: 
  // line (p2,x1) cross l1 = 0
  double t2 = 0;
 
  if( fabs((n12_unit.y*l1_unit.z-n12_unit.z*l1_unit.y))>delta ){
    t2 = ((p1.z-z1)*l1_unit.y-(p1.y-y1)*l1_unit.z) / (n12_unit.y*l1_unit.z-n12_unit.z*l1_unit.y);
  }  else if ( fabs((n12_unit.x*l1_unit.y-n12_unit.y*l1_unit.x))>delta ) {
    t2 = ((p1.y-y1)*l1_unit.x-(p1.x-x1)*l1_unit.y) / (n12_unit.x*l1_unit.y-n12_unit.y*l1_unit.x);
}  else if ( fabs((n12_unit.z*l1_unit.x-n12_unit.x*l1_unit.z))>delta ) {
    t2 = ((p1.x-x1)*l1_unit.z-(p1.z-z1)*l1_unit.x) / (n12_unit.z*l1_unit.x-n12_unit.x*l1_unit.z);
  }
  //cout<<"t2 = "<<t2<<endl;
  struct point p2;
  p2.x = p1.x+t2*n12_unit.x;
  p2.y = p1.y+t2*n12_unit.y;
  p2.z = p1.z+t2*n12_unit.z;
  //cout<<"p2 = ("<<p2.x<<","<<p2.y<<","<<p2.z<<")"<<endl;

  struct Line lper;
  lper.p1 = p1;
  lper.p2 = p2;
  lper.status=true;

  pub_dca_x = (lper.p1.x + lper.p2.x) / 2;
  pub_dca_y = (lper.p1.y + lper.p2.y) / 2;
  pub_dca_z = (lper.p1.z + lper.p2.z) / 2;
  double l12x=lper.p1.x-lper.p2.x;
  double l12y=lper.p1.y-lper.p2.y;
  double l12z=lper.p1.z-lper.p2.z;
 
  pub_dca = sqrt(l12x*l12x+l12y*l12y+l12z*l12z);

  return lper;
};

struct point cal_PoCA(struct Line l1, struct Line l2) {
  //printf("\nZhangzhijun =>\n");
  double x1 = l1.p1.x;
  double y1 = l1.p1.y;
  double z1 = l1.p1.z;
  double x2 = l1.p2.x;
  double y2 = l1.p2.y;
  double z2 = l1.p2.z;
  double x3 = l2.p1.x;
  double y3 = l2.p1.y;
  double z3 = l2.p1.z;
  double x4 = l2.p2.x;
  double y4 = l2.p2.y;
  double z4 = l2.p2.z;

  struct point v_in, v_out;
  v_in.x = x2-x1;
  v_in.y = y2-y1;
  v_in.z = z2-z1;
  v_out.x = x4-x3;
  v_out.y = y4-y3;
  v_out.z = z4-z3;
  //printf("v_in = (%f %f %f)\n",v_in.x,v_in.y,v_in.z);
  //printf("v_out = (%f %f %f)\n",v_out.x,v_out.y,v_out.z);

  struct point v_n = cal_normal_vector(l1, l2);
  v_n = normalization(v_n);
  //printf("v_n = (%f %f %f)\n",v_n.x,v_n.y,v_n.z);

  double d = (x3-x2)*v_n.x + (y3-y2)*v_n.y + (z3-z2)*v_n.z;
  //printf("d = %f\n",d);

  double t_i = (v_out.x * (d * v_n.y + y2 - y3) - v_out.y * (d * v_n.x + x2 - x3)) / (v_out.x * v_in.y - v_in.x * v_out.y);
  //printf("t_i = %f\n",t_i);

  struct point PoCA;
  PoCA.x = x2 - (t_i * v_in.x) + (0.5 * d * v_n.x);
  PoCA.y = y2 - (t_i * v_in.y) + (0.5 * d * v_n.y);
  PoCA.z = z2 - (t_i * v_in.z) + (0.5 * d * v_n.z);
  //printf("PoCA = (%f %f %f)\n",PoCA.x,PoCA.y,PoCA.z);

  return PoCA;
};

double cal_DCA(struct Line l1, struct Line l2) {
  double x1 = l1.p1.x;
  double y1 = l1.p1.y;
  double z1 = l1.p1.z;
  double x2 = l1.p2.x;
  double y2 = l1.p2.y;
  double z2 = l1.p2.z;
  double x3 = l2.p1.x;
  double y3 = l2.p1.y;
  double z3 = l2.p1.z;
  double x4 = l2.p2.x;
  double y4 = l2.p2.y;
  double z4 = l2.p2.z;

  struct point v_in, v_out;
  v_in.x = x2-x1;
  v_in.y = y2-y1;
  v_in.z = z2-z1;
  v_out.x = x4-x3;
  v_out.y = y4-y3;
  v_out.z = z4-z3;

  struct point v_n = cal_normal_vector(l1, l2);
  v_n = normalization(v_n);

  double d = (x3-x2)*v_n.x + (y3-y2)*v_n.y + (z3-z2)*v_n.z;

  return d;
};

struct point cal_normal_vector(struct Line l1, struct Line l2){
  double x1 = l1.p1.x;
  double y1 = l1.p1.y;
  double z1 = l1.p1.z;
  double x2 = l1.p2.x;
  double y2 = l1.p2.y;
  double z2 = l1.p2.z;
  double x3 = l2.p1.x;
  double y3 = l2.p1.y;
  double z3 = l2.p1.z;
  double x4 = l2.p2.x;
  double y4 = l2.p2.y;
  double z4 = l2.p2.z;

  double x21 = x2-x1;
  double y21 = y2-y1;
  double z21 = z2-z1;
  double x43 = x4-x3;
  double y43 = y4-y3;
  double z43 = z4-z3;

  double xn = (y21*z43-y43*z21);
  double yn = -(x21*z43-x43*z21);
  double zn = (x21*y43-x43*y21);

  struct point normal_direction;
  normal_direction.x = xn;
  normal_direction.y = yn;
  normal_direction.z = zn;

  return normal_direction;
};

struct point normalization(struct point p1){
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  //cout<<"x1 = "<<x1<<" , y1 = "<<y1<<" , z1 = "<<z1<<endl;

  double length = sqrt(x1*x1+y1*y1+z1*z1);
  //cout<<"length = "<<length<<endl;
  struct point unit;
  if(length>delta){
    unit.x = x1/length;
    unit.y = y1/length;
    unit.z = z1/length;
    unit.status = true;
  }else{
    unit.x = 1;
    unit.y = 0;
    unit.z = 0;
    unit.status = false;
  }

  return unit;
};

#endif

