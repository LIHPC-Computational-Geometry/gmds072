#include <iostream>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/CAD/GeomVolume.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
int main() {
        MeshModel mod = DIM3|N|F|F2N;
        IGMesh plc(mod);

        Node n0 = plc.newNode(0,0,0);
        Node n1 = plc.newNode(0,10,0);
        Node n2 = plc.newNode(10,10,0);
        Node n3 = plc.newNode(10,0,0);

        Node n4 = plc.newNode(0,0,10);
        Node n5 = plc.newNode(0,10,10);
        Node n6 = plc.newNode(10,10,10);
        Node n7 = plc.newNode(10,0,10);

        IGMesh::surface& surf1 = plc.newSurface("s1");
        surf1.add(plc.newTriangle(n0,n1,n2));
        surf1.add(plc.newTriangle(n0,n2,n3));

        IGMesh::surface& surf2 = plc.newSurface("s2");
        surf2.add(plc.newTriangle(n4,n6,n5));
        surf2.add(plc.newTriangle(n4,n7,n6));

        IGMesh::surface& surf3 = plc.newSurface("s3");
        surf3.add(plc.newTriangle(n0,n1,n5));
        surf3.add(plc.newTriangle(n0,n5,n4));

        IGMesh::surface& surf4 = plc.newSurface("s4");
        surf4.add(plc.newTriangle(n1,n2,n6));
        surf4.add(plc.newTriangle(n1,n6,n5));

        IGMesh::surface& surf5 = plc.newSurface("s5");
        surf5.add(plc.newTriangle(n2,n3,n7));
        surf5.add(plc.newTriangle(n2,n7,n6));

        IGMesh::surface& surf6 = plc.newSurface("s6");
        surf6.add(plc.newTriangle(n3,n0,n4));
        surf6.add(plc.newTriangle(n3,n4,n7));


        IGMesh vol_mesh(DIM3|R|F|F2N|N|R2N);

        return 0;
}
