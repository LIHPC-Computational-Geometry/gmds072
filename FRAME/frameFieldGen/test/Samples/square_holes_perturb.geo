// Gmsh project created on Fri Apr  3 15:03:32 2015
Point(1) = {10, 0, 0, 0.5};
Point(2) = {10, 10, 0, 0.5};
Point(3) = {0, 10, 0, 0.5};
Point(4) = {0, 0, 0, 0.5};
Point(5) = {5, 2.6, 0, 0.5};
Point(6) = {5, 7.4, 0, 0.5};
Point(7) = {2.5, 4.9, 0, 0.5};
Point(8) = {7.5, 5.1, 0, 0.5};
Line(1) = {3, 2};
Line(2) = {2, 1};
Line(3) = {1, 4};
Line(4) = {4, 3};
Line(5) = {7, 6};
Line(6) = {6, 8};
Line(7) = {8, 5};
Line(8) = {5, 7};
Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(11) = {9, 10};