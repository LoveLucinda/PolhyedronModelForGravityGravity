

#include "CoordinateTrans.h"

EllipsolidPar refell(int type = ELL_WGS84)
{
	EllipsolidPar temp_ellpsolidPar;
	switch (type)
	{
	case ELL_WGS84:
		temp_ellpsolidPar.a = 6378137.0;
		temp_ellpsolidPar.finv = 298.257223563;
		break;
	case ELL_WGS72:
		temp_ellpsolidPar.a = 6378135.0;
		temp_ellpsolidPar.finv = 298.26;
		break;
	default:
		break;
	}
	double f = 1 / temp_ellpsolidPar.finv;
	temp_ellpsolidPar.b = temp_ellpsolidPar.a*(1 - f);
	temp_ellpsolidPar.e2 = 1 - (1 - f) *(1 - f);

	return temp_ellpsolidPar;
}

POINT3D_Descart Ell2XYZ(POINT3D_Spherical point3d , EllipsolidPar ellpar)
{
	POINT3D_Descart temp_Point3d_Descart;
	point3d.lat = deg2rad(point3d.lat);//×ª»»Îª»¡¶È
	point3d.lon = deg2rad(point3d.lon);

	double v = ellpar.a / sqrt(1 - ellpar.e2*sin(point3d.lat)*sin(point3d.lat));
	temp_Point3d_Descart.x = (v + point3d.h)*cos(point3d.lat)*cos(point3d.lon);
	temp_Point3d_Descart.y = (v + point3d.h)*cos(point3d.lat)*sin(point3d.lon);
	temp_Point3d_Descart.z = (v*(1 - ellpar.e2) + point3d.h)*sin(point3d.lat);

	temp_Point3d_Descart.x2 = (v + point3d.h2)*cos(point3d.lat)*cos(point3d.lon);
	temp_Point3d_Descart.y2 = (v + point3d.h2)*cos(point3d.lat)*sin(point3d.lon);
	temp_Point3d_Descart.z2 = (v*(1 - ellpar.e2) + point3d.h2)*sin(point3d.lat);

	return temp_Point3d_Descart;
}

void ELL2XYZ(const int number_point, POINT3D_Spherical *point3d_sphere, EllipsolidPar ellpar, POINT3D_Descart* point3d_descart)
{
	for (int i = 0; i < number_point; i++)
	{
		point3d_descart[i] = Ell2XYZ(point3d_sphere[i], ellpar);
	}
}

double deg2rad(double deg)
{
	return deg / 180.0*PI;
}
double rad2deg(double rad)
{
	return rad / PI * 180;
}