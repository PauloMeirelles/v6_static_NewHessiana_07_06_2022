//  Geometry basics, elementary entities

// Fator de escala para a malha em função do ponto construido
lc = 0.04;

// USE ANOTHER WAY TO CONTROLL MESH (SEE AGAIN THE MANUAL AND EXAMPLES, THEY HAVE THIS ANSWER)

// POINTS FOR CONSTRUCT COUNTOUR DOMAIN

Point(1) = {0, 0, 0, lc};
Point(2) = {0.0663, 0,  0, lc};
Point(3) = {0.0663, 2.0, 0, lc};
Point(4) = {0, 2.0, 0, lc};

// LINE CONSTRUCT FOR COUNTOUR DOMAIN

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// CONECTIVITY OF LINES

Curve Loop(1) = {1, 2, 3, 4};

// PLANE COUNTOUR OF DOMAIN

Plane Surface(1) = {1};

// FORMAT TO OUTPUT (SU.2)
//Mesh.Format=42;
//Save "modelo_teste.su2";




//+
//Transfinite Curve {4, 2} = 5 Using Progression 1;
//Transfinite Curve {3, 1} = 45 Using Progression 1;
//+

//+
//Physical Curve("CoundCountourNotApplied") = {1};
//Physical Curve("CoundCountourNotApplied") = {3};
//Physical Curve("CoundCountourApplied") = {2};
//Physical Curve("CoundCountourApplied") = {4};
//+
//Physical Surface("Suface") = {1};
//+

//Transfinite Curve {1,3} = 1 Using Progression 1;

//Transfinite Surface{1};
