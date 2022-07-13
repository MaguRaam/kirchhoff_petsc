//+Points:
Point(1) = {0.0, 0.0, -0.342, 1.0};
Point(2) = {0.22952, 0.0, -0.342, 1.0};
Point(3) = {-0.22952, 0.0, -0.342, 1.0};
Point(4) = {0.0, 0.22952, -0.342, 1.0};
Point(5) = {0.0, -0.22952, -0.342, 1.0};

//+Circle
Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 2};
Circle(3) = {2, 1, 5};
Circle(4) = {5, 1, 3};

//+extrude cirlce along z = 4
Extrude {0, 0, 0.68704} {
  Curve{1}; Curve{4}; Curve{3}; Curve{2}; 
}

//+create plane surface top, bottom and curved

Curve Loop(1) = {5, 17, 13, 9};
Plane Surface(21) = {1};

//+
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(22) = {2};


//+
Plane Surface(23) = {1};
