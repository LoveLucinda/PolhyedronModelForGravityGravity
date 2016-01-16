// VTK模型显示.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "vtkNew.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkMath.h"
//#include "vtkCubeAxesActor.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkSmartPointer.h"
#include "vtkTextProperty.h"
#include "vtkPoints.h"
#include "vtkHexahedron.h"
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkBoundingBox.h>
#include <vtkEarthSource.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include<iostream>
using namespace std;

//定义相机不同视觉
#define ID_CAMERA_GENERAL	1
#define ID_CAMERA_FRONT	2
#define ID_CAMERA_BACK	3
#define ID_CAMERA_LEFT	4
#define ID_CAMERA_RIGHT	5
#define ID_CAMERA_UP	6
#define ID_CAMERA_DOWN	7

#include "CoordinateTrans.h"
void DisplayModel();
void DisplayEarthModel(vtkSmartPointer<vtkActor> earthactor,double* bounds);
void DisplayPolyhedronModel(vtkSmartPointer<vtkPoints> hexahedronPoints, vtkSmartPointer<vtkUnstructuredGrid> aHexahedronGrid, vtkSmartPointer<vtkActor> polyhedronActor, POINT3D_Descart* point3d_descart, const int number_point);
int SetCamera(vtkSmartPointer<vtkRenderer> renderer, vtkBoundingBox boundingbox, int type);
void TestPolyhedronDisplay();
int _tmain(int argc, _TCHAR* argv[])
{
	cout << "hello model\n";

	DisplayModel();
	//TestPolyhedronDisplay();
	return 0;
}


#include <vtkCubeSource.h>
#include <vtkElevationFilter.h>
#include <vtkCellArray.h>
#include "vtkPolyhedron.h"
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkPointLocator.h>
#include <vtkShrinkFilter.h>

#define compare_doublevec(x, y, e) \
	(((x[0] - y[0])<e) && ((x[0] - y[0])>-e) && \
	((x[1] - y[1])<e) && ((x[1] - y[1])>-e) && \
	((x[2] - y[2])<e) && ((x[2] - y[2])>-e))
#define compare_double(x, y, e) ((x)-(y)<e && (x)-(y)>-e)

void TestPolyhedronDisplay()
{
	vtkNew<vtkCamera> camera;
	camera->SetClippingRange(1.0, 100.0);
	camera->SetFocalPoint(1.26612, -0.81045, 1.24353);
	camera->SetPosition(-5.66214, -2.58773, 11.243);
	vtkNew<vtkRenderer> renderer;
	renderer->SetActiveCamera(camera.GetPointer());
	vtkNew<vtkRenderWindow> renWin;
	renWin->SetMultiSamples(0);
	renWin->AddRenderer(renderer.GetPointer());
	renWin->SetWindowName("Digital Earth");
	renWin->SetSize(600, 600);
	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetRenderWindow(renWin.GetPointer());
	renderer->SetBackground(0.2, 0.4, 0.1);
	//--------------------------------------------------------
	// create the a cube
	vtkSmartPointer<vtkCubeSource> cube =
		vtkSmartPointer<vtkCubeSource>::New();
	cube->SetXLength(10);
	cube->SetYLength(10);
	cube->SetZLength(20);
	cube->SetCenter(0, 0, 0);
	cube->Update();

	// add scaler
	vtkSmartPointer<vtkElevationFilter> ele =
		vtkSmartPointer<vtkElevationFilter>::New();
	ele->SetInputConnection(cube->GetOutputPort());
	ele->SetLowPoint(0, 0, -10);
	ele->SetHighPoint(0, 0, 10);
	ele->Update();
	vtkPolyData* poly = vtkPolyData::SafeDownCast(ele->GetOutput());

	// create a test polyhedron
	vtkIdType pointIds[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType face0[4] = { 0, 2, 6, 4 };
	vtkIdType face1[4] = { 1, 3, 7, 5 };
	vtkIdType face2[4] = { 0, 1, 3, 2 };
	vtkIdType face3[4] = { 4, 5, 7, 6 };
	vtkIdType face4[4] = { 0, 1, 5, 4 };
	vtkIdType face5[4] = { 2, 3, 7, 6 };
	faces->InsertNextCell(4, face0);
	faces->InsertNextCell(4, face1);
	faces->InsertNextCell(4, face2);
	faces->InsertNextCell(4, face3);
	faces->InsertNextCell(4, face4);
	faces->InsertNextCell(4, face5);
	vtkSmartPointer<vtkUnstructuredGrid> ugrid0 =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	ugrid0->SetPoints(poly->GetPoints());
	ugrid0->GetPointData()->DeepCopy(poly->GetPointData());

	ugrid0->InsertNextCell(VTK_POLYHEDRON, 8, pointIds,
		6, faces->GetPointer());

	vtkPolyhedron *polyhedron = static_cast<vtkPolyhedron*>(ugrid0->GetCell(0));

	vtkCellArray * cell = ugrid0->GetCells();
	vtkIdTypeArray * pids = cell->GetData();
	cout << "num of cells: " << cell->GetNumberOfCells() << endl;
	cout << "num of tuples: " << pids->GetNumberOfTuples() << endl;
	for (int i = 0; i < pids->GetNumberOfTuples(); i++)
	{
		cout << pids->GetValue(i) << " ";
	}
	cout << endl;
	cell->Print(cout);
	// Print out basic information
	//cout << "Testing polyhedron is a cube of with bounds "
	//	<< "[-5, 5, -5, 5, -10, 10]. It has "
	//	<< polyhedron->GetNumberOfEdges() << " edges and "
	//	<< polyhedron->GetNumberOfFaces() << " faces." << endl;

	//double p1[3] = { -100, 0, 0 };
	//double p2[3] = { 100, 0, 0 };
	//double tol = 0.001;
	//double t, x[3], pc[3];
	//int subId = 0;
	//// test writer
	//vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
	//	vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	//writer->SetInputData(ugrid0);
	//writer->SetFileName("test.vtu");
	//writer->SetDataModeToAscii();
	//writer->Update();
	//cout << "finished writing the polyhedron mesh to test.vth " << endl;

	//// test reader
	//vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
	//	vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	//reader->SetFileName("test.vtu");
	//reader->Update();
	//cout << "finished reading the polyhedron mesh from test.vth " << endl;

	//vtkUnstructuredGrid * ugrid = reader->GetOutput();
	//polyhedron = vtkPolyhedron::SafeDownCast(ugrid->GetCell(0));

	//// write again to help compare
	//writer->SetInputData(ugrid);
	//writer->SetFileName("test1.vtu");
	//writer->SetDataModeToAscii();
	//writer->Update();
	//// test the polyhedron functions
	//// test intersection
	//int numInts = polyhedron->IntersectWithLine(p1, p2, tol, t, x, pc, subId); //should be 2
	//if (numInts != 2)
	//{
	//	cerr << "Expect 2 intersections, but get " << numInts << endl;
	//	return ;
	//}

	//// test inside
	//int inside = polyhedron->IsInside(p1, tol); //should be out
	//if (inside)
	//{
	//	cerr << "Expect point [" << p1[0] << ", " << p1[1] << ", " << p1[2]
	//		<< "] to be outside the polyhedral, but it's inside." << endl;
	//	return ;
	//}

	//p2[0] = 0.0; p2[1] = 0.0; p2[2] = 0.0;
	//inside = polyhedron->IsInside(p2, tol); //should be in
	//if (!inside)
	//{
	//	cerr << "Expect point [" << p2[0] << ", " << p2[1] << ", " << p2[2]
	//		<< "] to be inside the polyhedral, but it's outside." << endl;
	//	return ;
	//}
	//// test EvaluatePosition and interpolation function
	//double weights[8], closestPoint[3], dist2;

	//for (int i = 0; i < 8; i++)
	//{
	//	double v;
	//	poly->GetPointData()->GetScalars()->GetTuple(i, &v);
	//	cout << v << " ";
	//}
	//cout << endl;

	//// case 0: point on the polyhedron
	//x[0] = 5.0; x[1] = 0.0; x[2] = 0.0;
	//polyhedron->EvaluatePosition(x, closestPoint, subId, pc, dist2, weights);

	//cout << "weights for point ["
	//	<< x[0] << ", " << x[1] << ", " << x[2] << "]:" << endl;
	//for (int i = 0; i < 8; i++)
	//{
	//	cout << weights[i] << " ";
	//}
	//cout << endl;

	//double refWeights[8] = { 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25 };
	//for (int i = 0; i < 8; i++)
	//{
	//	if (!compare_double(refWeights[i], weights[i], 0.00001))
	//	{
	//		cout << "Error computing the weights for a point on the polyhedron."
	//			<< endl;
	//		return ;
	//	}
	//}

	//double refClosestPoint[3] = { 5.0, 0.0, 0.0 };
	//if (!compare_doublevec(closestPoint, refClosestPoint, 0.00001))
	//{
	//	cout << "Error finding the closet point of a point on the polyhedron."
	//		<< endl;
	//	return ;
	//}

	//double refDist2 = 0.0;
	//if (!compare_double(dist2, refDist2, 0.000001))
	//{
	//	cout << "Error computing the distance for a point on the polyhedron."
	//		<< endl;
	//	return ;
	//}
	//// case 1: point inside the polyhedron
	//x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
	//polyhedron->EvaluatePosition(x, closestPoint, subId, pc, dist2, weights);

	//cout << "weights for point ["
	//	<< x[0] << ", " << x[1] << ", " << x[2] << "]:" << endl;
	//for (int i = 0; i < 8; i++)
	//{
	//	cout << weights[i] << " ";
	//}
	//cout << endl;

	//double refWeights1[8] = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 };
	//for (int i = 0; i < 8; i++)
	//{
	//	if (!compare_double(refWeights1[i], weights[i], 0.00001))
	//	{
	//		cout << "Error computing the weights for a point inside the polyhedron."
	//			<< endl;
	//		return ;
	//	}
	//}

	//if (!compare_double(dist2, refDist2, 0.000001))
	//{
	//	cout << "Error computing the distance for a point inside the polyhedron."
	//		<< endl;
	//	return ;
	//}

	//// case 2: point outside the polyhedron
	//x[0] = 8.0; x[1] = 0.0; x[2] = 0.0;
	//polyhedron->EvaluatePosition(x, closestPoint, subId, pc, dist2, weights);

	//cout << "weights for point ["
	//	<< x[0] << ", " << x[1] << ", " << x[2] << "]:" << endl;
	//for (int i = 0; i < 8; i++)
	//{
	//	cout << weights[i] << " ";
	//}
	//cout << endl;

	//double refWeights2[8] = { 0.0307, 0.0307, 0.0307, 0.0307,
	//	0.2193, 0.2193, 0.2193, 0.2193 };
	//for (int i = 0; i < 8; i++)
	//{
	//	if (!compare_double(refWeights2[i], weights[i], 0.0001))
	//	{
	//		cout << "Error computing the weights for a point outside the polyhedron."
	//			<< endl;
	//		return ;
	//	}
	//}

	//if (!compare_doublevec(closestPoint, refClosestPoint, 0.00001))
	//{
	//	cout << "Error finding the closet point of a point outside the polyhedron."
	//		<< endl;
	//	return ;
	//}

	//refDist2 = 9.0;
	//if (!compare_double(dist2, refDist2, 0.000001))
	//{
	//	cout << "Error computing the distance for a point outside the polyhedron."
	//		<< endl;
	//	return ;
	//}

	//// test evaluation location
	//double weights1[8];
	//polyhedron->EvaluateLocation(subId, pc, x, weights1);

	//double refPoint[3] = { 8.0, 0.0, 0.0 };
	//if (!compare_doublevec(refPoint, x, 0.00001))
	//{
	//	cout << "Error evaluate the point location for its parameter coordinate."
	//		<< endl;
	//	return ;
	//}

	//for (int i = 0; i < 8; i++)
	//{
	//	if (!compare_double(refWeights2[i], weights1[i], 0.0001))
	//	{
	//		cout << "Error computing the weights based on parameter coordinates."
	//			<< endl;
	//		return ;
	//	}
	//}

	//// test derivative
	//pc[0] = 0;  pc[1] = 0.5;  pc[2] = 0.5;
	//polyhedron->EvaluateLocation(subId, pc, x, weights1);

	//double deriv[3], values[8];
	//vtkDataArray * dataArray = poly->GetPointData()->GetScalars();
	//for (int i = 0; i < 8; i++)
	//{
	//	dataArray->GetTuple(i, values + i);
	//}
	//polyhedron->Derivatives(subId, pc, values, 1, deriv);

	//cout << "derivative for point ["
	//	<< x[0] << ", " << x[1] << ", " << x[2] << "]:" << endl;
	//for (int i = 0; i < 3; i++)
	//{
	//	cout << deriv[i] << " ";
	//}
	//cout << endl;

	//double refDeriv[3] = { 0.0, 0.0, 0.05 };
	//if (!compare_doublevec(refDeriv, deriv, 0.00001))
	//{
	//	cout << "Error computing derivative for a point inside the polyhedron."
	//		<< endl;
	//	return ;
	//}
	//// test triangulation
	//vtkSmartPointer<vtkPoints> tetraPoints = vtkSmartPointer<vtkPoints>::New();
	//vtkSmartPointer<vtkIdList> tetraIdList = vtkSmartPointer<vtkIdList>::New();
	//polyhedron->Triangulate(0, tetraIdList, tetraPoints);

	//cout << endl << "Triangulation result:" << endl;

	//for (int i = 0; i < tetraPoints->GetNumberOfPoints(); i++)
	//{
	//	double *pt = tetraPoints->GetPoint(i);
	//	cout << "point #" << i << ": [" << pt[0] << ", "
	//		<< pt[1] << ", " << pt[2] << "]" << endl;
	//}

	//vtkIdType * ids = tetraIdList->GetPointer(0);
	//for (int i = 0; i < tetraIdList->GetNumberOfIds(); i += 4)
	//{
	//	cout << "tetra #" << i / 4 << ":" << ids[i] << " "
	//		<< ids[i + 1] << " " << ids[i + 2] << " " << ids[i + 3] << endl;
	//}

	//vtkSmartPointer<vtkUnstructuredGrid> tetraGrid =
	//	vtkSmartPointer<vtkUnstructuredGrid>::New();
	//for (int i = 0; i < tetraIdList->GetNumberOfIds(); i += 4)
	//{
	//	tetraGrid->InsertNextCell(VTK_TETRA, 4, ids + i);
	//}
	//tetraGrid->SetPoints(poly->GetPoints());
	//tetraGrid->GetPointData()->DeepCopy(poly->GetPointData());

	//// test contour
	//vtkSmartPointer<vtkPointLocator> locator =
	//	vtkSmartPointer<vtkPointLocator>::New();
	//vtkSmartPointer<vtkCellArray> resultPolys =
	//	vtkSmartPointer<vtkCellArray>::New();
	//vtkSmartPointer<vtkPointData> resultPd =
	//	vtkSmartPointer<vtkPointData>::New();
	//vtkSmartPointer<vtkCellData> resultCd =
	//	vtkSmartPointer<vtkCellData>::New();
	//vtkSmartPointer<vtkPoints> resultPoints =
	//	vtkSmartPointer<vtkPoints>::New();
	//resultPoints->DeepCopy(ugrid0->GetPoints());
	//locator->InitPointInsertion(resultPoints, ugrid0->GetBounds());

	//polyhedron->Contour(0.5, tetraGrid->GetPointData()->GetScalars(), locator,
	//	NULL, NULL, resultPolys,
	//	tetraGrid->GetPointData(), resultPd,
	//	tetraGrid->GetCellData(), 0, resultCd);

	//// output the contour
	//vtkSmartPointer<vtkUnstructuredGrid> contourResult =
	//	vtkSmartPointer<vtkUnstructuredGrid>::New();
	//contourResult->SetPoints(locator->GetPoints());
	//contourResult->SetCells(VTK_POLYGON, resultPolys);
	//contourResult->GetPointData()->DeepCopy(resultPd);

	//// test clip
	//vtkSmartPointer<vtkPointLocator> locator1 =
	//	vtkSmartPointer<vtkPointLocator>::New();
	//vtkSmartPointer<vtkCellArray> resultPolys1 =
	//	vtkSmartPointer<vtkCellArray>::New();
	//vtkSmartPointer<vtkPointData> resultPd1 =
	//	vtkSmartPointer<vtkPointData>::New();
	//vtkSmartPointer<vtkCellData> resultCd1 =
	//	vtkSmartPointer<vtkCellData>::New();
	//vtkSmartPointer<vtkPoints> resultPoints1 =
	//	vtkSmartPointer<vtkPoints>::New();
	//resultPoints1->DeepCopy(ugrid0->GetPoints());
	//locator1->InitPointInsertion(resultPoints1, ugrid0->GetBounds());

	//polyhedron->Clip(0.5, tetraGrid->GetPointData()->GetScalars(), locator1,
	//	resultPolys1, tetraGrid->GetPointData(), resultPd1,
	//	tetraGrid->GetCellData(), 0, resultCd1, 0);

	//// output the clipped polyhedron
	//vtkSmartPointer<vtkUnstructuredGrid> clipResult =
	//	vtkSmartPointer<vtkUnstructuredGrid>::New();
	//clipResult->SetPoints(locator1->GetPoints());
	//clipResult->SetCells(VTK_POLYHEDRON, resultPolys1);
	//clipResult->GetPointData()->DeepCopy(resultPd1);

	//// shrink to show the gaps between tetrahedrons.
	//vtkSmartPointer<vtkShrinkFilter> shrink =
	//	vtkSmartPointer<vtkShrinkFilter>::New();
	//shrink->SetInputData(tetraGrid);
	//shrink->SetShrinkFactor(0.7);
	
	// create actors
	vtkSmartPointer<vtkDataSetMapper> mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(poly);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	/*vtkSmartPointer<vtkDataSetMapper> contourMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	contourMapper->SetInputData(contourResult);

	vtkSmartPointer<vtkActor> contourActor =
		vtkSmartPointer<vtkActor>::New();
	contourActor->SetMapper(contourMapper);

	vtkSmartPointer<vtkDataSetMapper> clipPolyhedronMapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	clipPolyhedronMapper->SetInputData(clipResult);

	vtkSmartPointer<vtkActor> clipPolyhedronActor =
		vtkSmartPointer<vtkActor>::New();
	clipPolyhedronActor->SetMapper(clipPolyhedronMapper);*/

	// Create rendering infrastructure
	vtkSmartPointer<vtkProperty> prop = vtkSmartPointer<vtkProperty>::New();
	prop->LightingOff();
	prop->SetRepresentationToSurface();
	prop->EdgeVisibilityOn();
	prop->SetLineWidth(3.0);
	prop->SetOpacity(0.8);

	// set property
	//actor->SetProperty(prop);
	//contourActor->SetProperty(prop);
	//clipPolyhedronActor->SetProperty(prop);

	renderer->AddActor(actor);
	//renderer->AddActor(contourActor);
	//renderer->AddActor(clipPolyhedronActor);
	renderer->SetBackground(.5, .5, .5);
	//--------------------------------------------
	vtkBoundingBox boundbox(actor->GetBounds());
	SetCamera(renderer.GetPointer(), boundbox, ID_CAMERA_FRONT);
	renWin->Render();

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);
	iren->Start();
}
void DisplayModel()
{
	vtkNew<vtkCamera> camera;
	camera->SetClippingRange(1.0, 100.0);
	camera->SetFocalPoint(1.26612, -0.81045, 1.24353);
	camera->SetPosition(-5.66214, -2.58773, 11.243);
	vtkNew<vtkRenderer> renderer;
	renderer->SetActiveCamera(camera.GetPointer());
	vtkNew<vtkRenderWindow> renWin;
	renWin->SetMultiSamples(0);
	renWin->AddRenderer(renderer.GetPointer());
	renWin->SetWindowName("Digital Earth");
	renWin->SetSize(600, 600);
	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetRenderWindow(renWin.GetPointer());
	renderer->SetBackground(0,0,0);
	//地球模型
	vtkSmartPointer<vtkActor> earthActor = vtkSmartPointer<vtkActor>::New();
	renderer->AddActor(earthActor.GetPointer());
	//double bounds[4] = {340,80,20,80};
	double bounds[4] = { 0, 360, 0, 180 };
	DisplayEarthModel(earthActor.GetPointer(),bounds);
	//设置模型体的顶点坐标-----------------------------------------------------------------
	vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
	EllipsolidPar ellpar = refell(ELL_WGS84);
	const int number_point = 4;
	POINT3D_Spherical point3d_sphere[number_point];
	POINT3D_Descart point3d_descart[number_point];
	vtkSmartPointer<vtkActor> polyhedronActor = vtkSmartPointer<vtkActor>::New();
	renderer->AddActor(polyhedronActor.GetPointer());
	polyhedronActor->GetProperty()->SetDiffuseColor(0, 0.2, 0.9);
	//polyhedronActor->GetProperty()->SetOpacity(0.8);
	vtkSmartPointer<vtkUnstructuredGrid> aHexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	aHexahedronGrid->Allocate(1, 1);
	//-------------------------------------------------------------------------------------
	double lonmin = -20, lonmax = 80, latmin = -70, latmax = -10;
	double dlon = 5, dlat = 5;
	int n_lon = (int)((lonmax - lonmin) / dlon);
	int n_lat = (int)((latmax - latmin) / dlat);
	//
	point3d_sphere[0].h = 0;	point3d_sphere[0].h2 = -1500;
	point3d_sphere[1].h = 0;	point3d_sphere[1].h2 = -1500;
	point3d_sphere[2].h = 0;	point3d_sphere[2].h2 = -1500;
	point3d_sphere[3].h = 0;	point3d_sphere[3].h2 = -1500;
	for (int i = 0; i < n_lat; i++)
	{
		point3d_sphere[0].lat = latmin + i*dlat;
		point3d_sphere[1].lat = point3d_sphere[0].lat;
		point3d_sphere[2].lat = point3d_sphere[0].lat + dlat;
		point3d_sphere[3].lat = point3d_sphere[2].lat;
		for (int j = 0; j < n_lon; j++)
		{
			point3d_sphere[0].lon = lonmin + j*dlon;
			point3d_sphere[1].lon = point3d_sphere[0].lon + dlon;
			point3d_sphere[2].lon = point3d_sphere[1].lon;
			point3d_sphere[3].lon = point3d_sphere[0].lon;
			ELL2XYZ(number_point, point3d_sphere, ellpar, point3d_descart);
			DisplayPolyhedronModel(hexahedronPoints.GetPointer(), aHexahedronGrid.GetPointer(), polyhedronActor.GetPointer(), point3d_descart, number_point);
		}
		
	}
	aHexahedronGrid->SetPoints(hexahedronPoints);
	//point3d_sphere[0].lon = -20;	point3d_sphere[0].lat = -70;	point3d_sphere[0].h = 0;	point3d_sphere[0].h2 = -1500000;
	//point3d_sphere[1].lon = 0;	point3d_sphere[1].lat = -70;	point3d_sphere[1].h = 0;	point3d_sphere[1].h2 = -1500000;
	//point3d_sphere[2].lon = 0;	point3d_sphere[2].lat = -50;	point3d_sphere[2].h = 0;	point3d_sphere[2].h2 = -1500000;
	//point3d_sphere[3].lon = -20;	point3d_sphere[3].lat = -50;	point3d_sphere[3].h = 0;	point3d_sphere[3].h2 = -1500000;
	//ELL2XYZ(number_point, point3d_sphere, ellpar, point3d_descart);
	//DisplayPolyhedronModel(hexahedronPoints.GetPointer(), aHexahedronGrid.GetPointer(), polyhedronActor.GetPointer(), point3d_descart, number_point);

	//point3d_sphere[0].lon = 10;	point3d_sphere[0].lat = -40;	point3d_sphere[0].h = 0;	point3d_sphere[0].h2 = -1500000;
	//point3d_sphere[1].lon = 23;	point3d_sphere[1].lat = -40;	point3d_sphere[1].h = 0;	point3d_sphere[1].h2 = -1500000;
	//point3d_sphere[2].lon = 23;	point3d_sphere[2].lat = -30;	point3d_sphere[2].h = 0;	point3d_sphere[2].h2 = -1500000;
	//point3d_sphere[3].lon = 10;	point3d_sphere[3].lat = -30;	point3d_sphere[3].h = 0;	point3d_sphere[3].h2 = -1500000;
	//ELL2XYZ(number_point, point3d_sphere, ellpar, point3d_descart);
	//DisplayPolyhedronModel(hexahedronPoints.GetPointer(), aHexahedronGrid.GetPointer(), polyhedronActor.GetPointer(), point3d_descart, number_point);
	//--------------------------------------------------------------------------------------------
	vtkSmartPointer<vtkDataSetMapper> aHexahedronMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	aHexahedronMapper->SetInputData(aHexahedronGrid);
	polyhedronActor->SetMapper(aHexahedronMapper);
	polyhedronActor->GetProperty()->SetRepresentationToWireframe();
	polyhedronActor->GetProperty()->SetLineWidth(5);
	//刷新试图显示
	vtkBoundingBox boundbox(earthActor->GetBounds());
	SetCamera(renderer.GetPointer(), boundbox, ID_CAMERA_FRONT);
	renWin->Render();

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	iren->SetInteractorStyle(style);
	iren->Start();
}

void DisplayPolyhedronModel(vtkSmartPointer<vtkPoints> hexahedronPoints, vtkSmartPointer<vtkUnstructuredGrid> aHexahedronGrid,vtkSmartPointer<vtkActor> polyhedronActor, POINT3D_Descart* point3d_descart, const int number_point)
{
	//vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
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
	for (int i = 0; i < number_point*2; i++)
	{
		aHexahedron->GetPointIds()->SetId(i, i + number_points);
	}
	aHexahedronGrid->InsertNextCell(aHexahedron->GetCellType(), aHexahedron->GetPointIds());
	
	
}

void DisplayEarthModel(vtkSmartPointer<vtkActor> earthactor, double* bounds)
{
	double R_earth = 6371*1000;
	//vtkSmartPointer<vtkSphereSource>earthsource = vtkSmartPointer<vtkSphereSource>::New();
	//earthsource->SetCenter(0, 0, 0);
	//earthsource->SetRadius(R_earth);
	//earthsource->SetPhiResolution(100);
	//earthsource->SetThetaResolution(100);
	//earthsource->SetStartPhi(bounds[2]);
	//earthsource->SetEndPhi(bounds[3]);
	//earthsource->SetStartTheta(bounds[0]);
	//earthsource->SetEndTheta(bounds[1]);
	//vtkSmartPointer<vtkPolyDataMapper> earthMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//earthMapper->SetInputConnection(earthsource->GetOutputPort());
	//earthactor->SetMapper(earthMapper);
	//earthactor->GetProperty()->SetColor(0.1, 0, 0.8);
	//earthactor->GetProperty()->SetOpacity(0.9);
	vtkEarthSource* es = vtkEarthSource::New();
	es->SetRadius(R_earth);
	es->SetOnRatio(1);
	es->SetOutline(3);
	es->OutlineOn();
	vtkPolyDataMapper* earthmapper = vtkPolyDataMapper::New();
	earthmapper->SetInputConnection(es->GetOutputPort());
	earthactor->SetMapper(earthmapper);
	earthactor->GetProperty()->SetColor(1, 0, 0);
	earthactor->GetProperty()->SetLineWidth(4);
}

int SetCamera(vtkSmartPointer<vtkRenderer> renderer, vtkBoundingBox boundingbox, int type)
{
	double center[3];
	boundingbox.GetCenter(center);
	double bounds[6];
	boundingbox.GetBounds(bounds);
	double xlength = bounds[1] - bounds[0];
	double ylength = bounds[3] - bounds[2];
	double zlength = bounds[5] - bounds[4];
	double xyBiZhi = ylength / xlength;
	double viewup_x = 0.3;
	switch (type)
	{
	case ID_CAMERA_FRONT:
		center[1] = bounds[2];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0], bounds[2] - ylength, center[2]);//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
		break;
	case ID_CAMERA_BACK:
		center[1] = bounds[3];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0], bounds[3] + ylength, center[2]);//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
		break;
	case ID_CAMERA_LEFT:
		center[0] = bounds[0];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0] - xlength, center[1], center[2]);//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
		break;
	case ID_CAMERA_RIGHT:
		center[0] = bounds[1];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0] + xlength, center[1], center[2]);//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
		break;
	case ID_CAMERA_UP:
		center[2] = bounds[5];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0], center[1], center[2] + 2 * (xlength>ylength ? xlength : ylength));//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 1, 0);//相机“上”方向
		break;
	case ID_CAMERA_DOWN:
		center[2] = bounds[4];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0], center[1], center[2] - 2 * (xlength>ylength ? xlength : ylength));//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 1, 0);//相机“上”方向
		break;
	default:
		center[2] = bounds[5];
		renderer->GetActiveCamera()->SetFocalPoint(center);//焦点
		renderer->GetActiveCamera()->SetPosition(center[0], bounds[2] - 1.5*ylength, bounds[5] + ylength / 2.0);//相机位置
		renderer->GetActiveCamera()->SetViewUp(0.0, 0, 1);//相机“上”方向
		break;
	}
	return 0;
}

