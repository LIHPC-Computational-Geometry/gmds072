// Gmsh project created on Wed Sep 30 10:05:07 2015
Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {-10, 0, 0, 1.0};
Point(4) = {-5, 0, 0, 1.0};
Point(5) = {5, 0, 0, 1.0};
Line(1) = {4, 3};
Line(2) = {5, 2};
Circle(3) = {5, 1, 4};
Circle(4) = {2, 1, 3};
Line Loop(5) = {3, 1, -4, -2};
Plane Surface(6) = {5};
