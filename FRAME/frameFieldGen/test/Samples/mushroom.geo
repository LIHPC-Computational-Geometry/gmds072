// Gmsh project created on Wed May 27 12:37:27 2015
Point(1) = {0, 0, 0, 0.1};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {-1, 0, 0, 0.1};
Point(4) = {-2, 0, 0, 0.1};
Point(5) = {2, 0, 0, 0.1};
Point(6) = {2, 2, 0, 0.1};
Point(7) = {-2, 2, 0, 0.1};
Point(8) = {2, 2., 0, 0.1};
Line(1) = {2, 1};
Line(2) = {3, 4};
Line(3) = {4, 7};
Line(4) = {2, 5};
Line(5) = {5, 6};
Circle(6) = {2, 1, 3};
Point(9) = {3, 2, 0, 0.1};
Point(10) = {-3, 2, 0, 0.1};
Point(11) = {0, 2, 0, 0.1};
Circle(7) = {9, 11, 10};
Line(8) = {10, 7};
Line(9) = {6, 9};
Delete {
  Line{1};
}
Line Loop(10) = {7, 8, -3, -2, -6, 4, 5, 9};
Plane Surface(11) = {10};
