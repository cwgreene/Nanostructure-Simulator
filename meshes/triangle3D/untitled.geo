Point(1) = {0, 0, 0};
Point(2) = {1.0, 0, 0};
Point(3) = {.5, .86602540378, 0};
Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line Loop(4) = {1, 2, 3};
Plane Surface(5) = {4};
Point(4) = {0, 0, 1};
Point(5) = {0.5, 0.8660254, 1};
Point(6) = {1, 0, 1};
Line(6) = {2, 6};
Line(7) = {4, 1};
Line(8) = {3, 5};
Line(9) = {5, 6};
Line(10) = {6, 4};
Line(11) = {4, 5};
Line Loop(12) = {6, 10, 7, 2};
Plane Surface(13) = {12};
Line Loop(14) = {10, 11, 9};
Plane Surface(15) = {14};
Line Loop(16) = {8, -11, 7, -1};
Plane Surface(17) = {16};
Line Loop(18) = {3, 8, 9, -6};
Plane Surface(19) = {18};
Surface Loop(20) = {19, 5, 17, 15, 13};
Volume(21) = {20};
