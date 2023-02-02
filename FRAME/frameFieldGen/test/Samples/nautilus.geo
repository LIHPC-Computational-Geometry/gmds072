// Gmsh project created on Fri Apr 10 09:38:29 2015
Point(1) = {0, 0, 0, 0.25};
Point(2) = {1, 0, 0, 0.25};
Point(3) = {-1, 0, 0, 0.25};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Point(4) = {-2, 0, 0, 0.25};
Point(5) = {2, 0, 0, 0.25};
Circle(3) = {4, 1, 5};
Point(6) = {-0.5, 0, 0, 0.25};
Point(7) = {-3, 0, 0, 0.25};
Circle(4) = {5, 6, 7};
Line(5) = {7, 4};
Line Loop(6) = {4, 5, 3};
Line Loop(7) = {1, 2};
Plane Surface(8) = {6, 7};
