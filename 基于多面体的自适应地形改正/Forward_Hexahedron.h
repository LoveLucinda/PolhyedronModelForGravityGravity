
#ifndef FORWARD_HEXAHEDRON
#define FORWARD_HEXAHEDRON

#include "vtkSmartPointer.h"
//#include "vtkTextProperty.h"
//#include "vtkRenderingCore_AUTOINIT_vtkInteractionStyle_vtkRenderingFreeType_vtkRenderingFreeTypeOpenGL_vtkRenderingOpenGL.h"
//#include "vtkNew.h"
//#include "vtkCamera.h"
//#include "vtkLight.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkMath.h"
//#include "vtkCubeAxesActor.h"
//#include "vtkProperty.h"
//#include "vtkInteractorStyleTrackballCamera.h"

//#include <vtkUnstructuredGrid.h>
//#include <vtkDataSetMapper.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkBoundingBox.h>
//#include <vtkEarthSource.h>
//#include <vtkPointData.h>
//#include <vtkCellData.h>
#include "vtkPoints.h"
#include "vtkHexahedron.h"
#include "vtkPoints.h"
#include <vtkUnstructuredGrid.h>
#include <vtkPolygon.h>

#include "CoordinateTrans.h"

#define G 6.67E-11

void Forward_V(vtkSmartPointer<vtkHexahedron> hexahearon);

#endif