#include<iostream>
#include<math.h>
#include<fstream>	// for file read
#include<stdlib.h>	// for casting str as int
#include<string>	// for getline / data feed
#include<vector>	// for data feed

#include <armadillo>

#include "structure.h"

using namespace std;
using namespace arma;

void parseInput(string filename, Structure *input);
Mat<double>createLocalStiffness(Element *e);
size_t split(const string &txt, vector<string> &strs, char ch);
int dirConvert (string dir);

Mat<double> constrainK(Mat<double> K, Structure *system) {
	// Modify Force Vector and Stiffness Matrix for given boundary condition
	Mat<double> Kcopy = K;
	int nBounds = (int)system->bounds.size();
	int nDof = (int)system->points.size()*3;
	for (int iBound = 0; iBound < nBounds; iBound++) {
		int dof = system->bounds[iBound].dof;
		for (int i = 0; i < nDof; i++) {
			system->globalLoad(i) -= Kcopy(dof,i)*system->bounds[iBound].disp;
			Kcopy(i,dof) = 0;
			Kcopy(dof,i) = 0;
		}
		system->globalLoad(dof) = system->bounds[iBound].disp;
		Kcopy(dof,dof) = 1;
	}
	return Kcopy;
}

int main(int argc, char *argv[]) {
	if (argc != 2)
		throw runtime_error("Please specify input file");
	Structure system;
	parseInput(argv[1], &system);
	int nEle = (int)system.elements.size();
	int nDof = (int)system.points.size()*3;
	Mat<double> K = zeros(nDof,nDof);
	
	for (int e = 0; e < nEle; e++) {
		Mat<double> KBar;
		KBar = createLocalStiffness(&system.elements[e]);
		for (int localI = 0; localI < 6; localI++) {
			for (int localJ = 0; localJ < 6; localJ++) {
				int globalI = system.elements[e].IDArray[localI];
				int globalJ = system.elements[e].IDArray[localJ];
				K(globalI, globalJ) += KBar(localI, localJ);
			}
		}
	}
	Mat<double> Kbounded = constrainK(K, &system);
	Mat<double> u = solve(Kbounded, system.globalLoad);
	cout << "Solved u\n" << u << endl;
	/* TODO
	 * Calculate Element forces from nodal disp & Solved Force Vector
	 */
	return 0;
}

void parseInput(string filename, Structure *input) {
	int nodeNum = 0;
	
	vector<string> data;
	int inputBlock = -1;
	string line;
	
	ifstream myfile (filename);
	if (myfile.is_open()) {
		while ( getline (myfile,line) ) {
			if (line == "POINTS") {
				inputBlock = 0;
				continue;
			} else if (line == "ELEMENTS") {
				inputBlock = 1;
				input->globalLoad.zeros((int)input->points.size()*3);
				continue;
			} else if (line == "BOUNDS") {
				inputBlock = 2;
				continue;
			} else if (line == "LOADING") {
				inputBlock = 3;
				continue;
			} else if (line == "") {
				inputBlock = -1;
				continue;
			};
			split(line, data, ' ');
			switch(inputBlock) {
				case 0: {
					// Points
					double x = atof(data[0].c_str());
					double y = atof(data[1].c_str());
					Point tmpPoint {nodeNum, x, y};
					input->points.push_back(tmpPoint);
					nodeNum++;
					break;
				}
				case 1: {
					// Elements
					int p1 = stoi(data[0]);
					int p2 = stoi(data[1]);
					double E = atof(data[2].c_str());
					double I = atof(data[3].c_str());
					double A = atof(data[4].c_str());
					Element tmpElement {&input->points[p1], &input->points[p2], E, I, A};
					input->elements.push_back(tmpElement);
					break;
				}
				case 2: {
					// Bounds
					int p = stoi(data[0]);
					int dof = input->points[p].n*3 + dirConvert(data[1]);
					double disp = atof(data[2].c_str());
					Bound tmpBound {dof, disp};
					input->bounds.push_back(tmpBound);
					break;
				}
				case 3: {
					// Loading
					vector<double> tmpLoad;
					if (data[0] == "ELEMENT") {
						Element *tmpElement = &input->elements[stoi(data[1])];
						double L = tmpElement->L;
						// tmpLoad = {F1, M1, F2, M2}
						if (data[2] == "DISTRIBUTED") {
							double w1 = atof(data[3].c_str());
							double w2 = atof(data[4].c_str());
							tmpLoad = { (7*w1+3*w2)*L/20, (w1/2+w2/3)*pow(L,2)/10, \
										(3*w1+7*w2)*L/20,-(w1/3+w2/2)*pow(L,2)/10};
						} else if (data[2] == "POINT") {
							double P = atof(data[3].c_str());
							double a = atof(data[4].c_str());
							double b = L-a;
							tmpLoad = { P*pow(b,2)*(3*a+b)/pow(L,3), P*a*pow(b,2)/pow(L,2), \
										P*pow(a,2)*(3*b+a)/pow(L,3),-P*b*pow(a,2)/pow(L,2)};
						} else {
							throw runtime_error("ELEMENT type loading must be either 'DISTRIBUTED' or 'POINT'");
						}
						// Add Local Element Vector to Global Force Vector
						input->addElementLoad(tmpElement, tmpLoad);
					} else {
						throw runtime_error("Can only handle Loading of type 'ELEMENT'");
					}
				}
			}
		}
		myfile.close();
	}
}

int dirConvert (string dir) {
	// Convert 'X' 'Y' or 'R' into 0 1 2
	int dof;
	if (dir == "X") {
		dof = 0;
	} else if (dir == "Y") {
		dof = 1;
	} else if (dir == "R") {
		dof = 2;
	} else {
		throw runtime_error("Direction must be 'X', 'Y', or 'Z'");
	}
	return dof;
}

Mat<double> createLocalStiffness(Element *e) {
	// Dereference into local variables for clarity
	double E = e->E;
	double I = e->I;
	double A = e->A;
	double L = e->L;
	// 1D Local Stiffness Matrix
	Mat<double>localK = {	{E*A/L,       0,       0},
							{    0, 4*E*I/L, 2*E*I/L},
							{    0, 2*E*I/L, 4*E*I/L}};
	// Rigid Body Matrix
	Mat<double>RB = {	{-1,   0, 0, 1,    0, 0},
						{ 0, 1/L, 1, 0, -1/L, 0},
						{ 0, 1/L, 0, 0, -1/L, 1}};
	// Rotational Matrix
    double c = cos(e->theta);
    double s = sin(e->theta);
	Mat<double>ROT = {	{ c, s, 0, 0, 0, 0},
						{-s, c, 0, 0, 0, 0},
						{ 0, 0, 1, 0, 0, 0},
						{ 0, 0, 0, c, s, 0},
						{ 0, 0, 0,-s, c, 0},
						{ 0, 0, 0, 0, 0, 1}};
	// Combine into Global Domain
	Mat<double> kBar = ROT.t() * RB.t() * localK * RB * ROT;
	return kBar;
}

size_t split(const string &txt, vector<string> &strs, char ch) {
	// from https://stackoverflow.com/questions/5888022/split-string-by-single-spaces
	// used to parse input file lines
    size_t pos = txt.find( ch );
    size_t initialPos = 0;
    strs.clear();

    // Decompose statement
    while( pos != std::string::npos ) {
        strs.push_back( txt.substr( initialPos, pos - initialPos ) );
        initialPos = pos + 1;
        pos = txt.find( ch, initialPos );
    }

    // Add the last one
    strs.push_back( txt.substr( initialPos, min( pos, txt.size() ) - initialPos + 1 ) );
    return strs.size();
}