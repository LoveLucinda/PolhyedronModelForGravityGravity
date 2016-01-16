

#include "Forward_Hexahedron.h"



void Forward_V(vtkSmartPointer<vtkHexahedron> hexahearon)
{
	cout << hexahearon->GetNumberOfEdges() << endl;
	//--------------------------------------------------------
	double density = 6.27;//g/cm3
	//--------------------------------------------------------
	double V = 0; 
	vtkPolygon *polygon = static_cast<vtkPolygon*>(hexahearon->GetFace(1));
	//����㵽�������ľ���
	double P[3] = {100,100,0};
	double closestP[3];
	for (int i = 0; i < polygon->GetNumberOfPoints(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			cout << polygon->GetPoints()->GetPoint(i)[j]<<"\t";
		}
		cout << endl;
	}
	double distances=polygon->DistanceToPolygon(P, polygon->GetNumberOfPoints(),polygon->GetPoints()->GetPoint(0), polygon->GetBounds(), closestP);
	cout << distances << endl;
}