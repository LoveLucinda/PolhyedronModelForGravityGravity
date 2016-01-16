

#ifndef COORDINATETRANS
#define COORDINATETRANS
#include <iostream>
using namespace std;

#define PI 3.141592653

#define ELL_WGS84 1
#define ELL_WGS72 2

struct EllipsolidPar
{
	//% a - major semi - axis(m)
	//% b - minor semi - axis(m)
	//% e2 - eccentricity squared
	//% finv - inverse of flattening
	double a, b, e2, finv;
};
struct POINT3D_Descart
{
	double x, y, z;//��λ��m
	double x2, y2, z2;//��h2��Ӧ
};
struct POINT3D_Spherical
{
	double lat, lon, h;//��λ���ȣ�m
	double h2;//m//hΪ�ϴ�İ뾶��hΪ��С�İ뾶���̣߳�
};

EllipsolidPar refell(int type);
POINT3D_Descart Ell2XYZ(POINT3D_Spherical point3d, EllipsolidPar ellpar);
void ELL2XYZ(const int number_point,POINT3D_Spherical *point3d_sphere, EllipsolidPar ellpar, POINT3D_Descart* point3d_descart);
double deg2rad(double deg);
double rad2deg(double rad);

#endif