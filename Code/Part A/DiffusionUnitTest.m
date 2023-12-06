clear;close all;
 
%% Input parameter for test

% Input parameter

xmax = 1;  % domain

xmin = 0;

tstart =0;

tend   =1;

T = tstart - tend;

D=1;
 
%% Test 1: Check the value of local K

Ne   = 3;  % element steps

% Initialise Mesh

mesh = OneDimLinearMeshGen(xmin,xmax,Ne);

eID=1;

Lamda = 1;

J = mesh.elem(eID).J;

f =1;

gq = CreateGQScheme(2);

[matrix, diffTerm] = LocalStiffnessMatrix(eID,mesh,D,Lamda, gq);
diffTerm2 = LaplaceElemMatrix(D, eID, mesh);


result = [3*D -3*D; -3*D 3*D] - [2/3*Lamda*J 1/3*Lamda*J;1/3*Lamda*J 2/3*Lamda*J];

result2 = [3*D -3*D; -3*D 3*D];

diff = (diffTerm2-result2);

assert(sum(sum(diff*diff)) <= 1e-14)