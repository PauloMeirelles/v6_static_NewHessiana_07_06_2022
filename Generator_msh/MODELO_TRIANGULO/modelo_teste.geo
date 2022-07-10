//  Geometry basics, elementary entities

// Fator de escala para a malha em função do ponto construido
lc = 12;

// USE ANOTHER WAY TO CONTROLL MESH (SEE AGAIN THE MANUAL AND EXAMPLES, THEY HAVE THIS ANSWER)

// POINTS FOR CONSTRUCT COUNTOUR DOMAIN

Point(1) = {0, 0, 0, lc};
Point(2) = {6, 0,  0, lc};
Point(3) = {0, 6, 0, lc};


// LINE CONSTRUCT FOR COUNTOUR DOMAIN

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};


// CONECTIVITY OF LINES

Curve Loop(1) = {1, 2, 3};

// PLANE COUNTOUR OF DOMAIN

Plane Surface(1) = {1};

// FORMAT TO OUTPUT (SU.2)
//Mesh.Format=42;
//Save "modelo_teste.su2";




//+
Transfinite Surface {1};
