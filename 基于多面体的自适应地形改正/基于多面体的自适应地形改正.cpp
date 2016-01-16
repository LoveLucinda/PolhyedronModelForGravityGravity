// 基于多面体的自适应地形改正.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

#include "Forward_Hexahedron.h"
void TestForward_Hexahedron();
int _tmain(int argc, _TCHAR* argv[])
{
	TestForward_Hexahedron();
	return 0;
}

void TestForward_Hexahedron()
{
	//设置模型体的顶点坐标-----------------------------------------------------------------
	vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
	EllipsolidPar ellpar = refell(ELL_WGS84);
	const int number_point = 4;
	POINT3D_Spherical point3d_sphere[number_point];
	POINT3D_Descart point3d_descart[number_point];
	//vtkSmartPointer<vtkUnstructuredGrid> aHexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//aHexahedronGrid->Allocate(1, 1);

	point3d_sphere[0].lon = 10;	point3d_sphere[0].lat = -40;	point3d_sphere[0].h = 0;	point3d_sphere[0].h2 = -1500000;
	point3d_sphere[1].lon = 23;	point3d_sphere[1].lat = -40;	point3d_sphere[1].h = 0;	point3d_sphere[1].h2 = -1500000;
	point3d_sphere[2].lon = 23;	point3d_sphere[2].lat = -30;	point3d_sphere[2].h = 0;	point3d_sphere[2].h2 = -1500000;
	point3d_sphere[3].lon = 10;	point3d_sphere[3].lat = -30;	point3d_sphere[3].h = 0;	point3d_sphere[3].h2 = -1500000;
	ELL2XYZ(number_point, point3d_sphere, ellpar, point3d_descart);

	int number_points = hexahedronPoints->GetNumberOfPoints();

	for (size_t i = 0; i < number_point; i++)
	{
		hexahedronPoints->InsertPoint(i + number_points, point3d_descart[i].x, point3d_descart[i].y, point3d_descart[i].z);
	}
	for (size_t i = 0; i < number_point; i++)
	{
		hexahedronPoints->InsertPoint(i + number_points + number_point, point3d_descart[i].x2, point3d_descart[i].y2, point3d_descart[i].z2);
	}
	vtkSmartPointer<vtkHexahedron>aHexahedron = vtkSmartPointer<vtkHexahedron>::New();
	for (int i = 0; i < number_point * 2; i++)
	{
		aHexahedron->GetPointIds()->SetId(i, i + number_points);
	}
	//aHexahedronGrid->InsertNextCell(aHexahedron->GetCellType(), aHexahedron->GetPointIds());

	//调用正演函数
	Forward_V(aHexahedron.GetPointer());
}

