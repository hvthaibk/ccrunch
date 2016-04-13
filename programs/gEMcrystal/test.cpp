/***********************************************************************
 *  File:       test.cpp
 *
 *  Author:     Lucía N. Martín Reixach
 *
 *  Contact:    luciamartinreixach@gmail.com
 *
 *  Copyright (C) 2013 Lucía N. Martín Reixach, Inria.
 **********************************************************************/

#include "cSegment3.hpp"
#include "cPolygon.hpp"
#include "inout.hpp"

using namespace std;
using namespace gem;

int main()
{

    cVector3<float>    a = cVector3<float>(-3, 0, 0);
    cVector3<float>    b = cVector3<float>(5, 0, 0);
    cVector3<float>    c = cVector3<float>(8, 0, 0);
    cVector3<float>    d = cVector3<float>(9.2, 1, 0);
        
    cSegment3<float>    ab = cSegment3<float>(a, b);
    cSegment3<float>    cd = cSegment3<float>(c, d);
    cVector3<float>     v = ab.getIntersectPoint(cd); 

    cout << "Coincide: " << ab.isCoincide(cd) << "\n";    
    //cout << "Los segmentos se intersectan: " << ab.isIntersect(cd) << "\n";
    //cout << "Los segmentos son paralelos: " << ab.isParallel(cd) << "\n";
    //cout << "Se intersectan en el punto: " << v[0] << " " << v[1] << " " << v[2] << "\n";

    cPolygon<float> poly = cPolygon<float>();
    poly.setNpoint(5);
    poly[0][0] = 1;
    poly[2][0] = -3;
    poly[2][1] = -2;
    poly[3][1] = 5;
    poly[0][1] = -4.89;
    poly[4][1] = 2; 
	
	cSegment3<float> *seg = poly.getSegments();

	//for (int i = 0; i < poly.getNpoint(); i++)
	//	cout << "Segmento " << i << ": (" << seg[i][0][0] << ", " << seg[i][0][1] << "), (" << seg[i][1][0] << ", " << seg[i][1][1] << ")\n";

	
    //cout << "Cant de vertices: " << poly.getNpoint() << "\n";
    //cout << "Promedio: " << poly.getMean()[0] << " " << poly.getMean()[1] << " " << poly.getMean()[2] << "\n";
    //cout << "MinX: " << poly.getCoordMinX() << "\n"; 
    //cout << "MaxX: " << poly.getCoordMaxX() << "\n";
    //cout << "MinY: " << poly.getCoordMinY() << "\n";
    //cout << "MaxY: " << poly.getCoordMaxY() << "\n";

	poly.opFill();

    return 0;
}
