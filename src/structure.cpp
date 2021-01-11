#include <vector>
#include <fstream>	// for file read
#include <stdlib.h>	// for casting str as int
#include <string>	// for getline / data feed

#include <armadillo>
#include "structure.h"

using namespace std;

Element::Element(Point* _p1, Point* _p2, double _E, double _I, double _A) {
	p1 = _p1;
	p2 = _p2;
	E = _E;
	I = _I;
	A = _A;
	L = sqrt(pow(p1->x - p2->x,2) + pow(p1->y - p2->y,2));
	theta = atan2(p2->y - p1->y, p2->x - p1->x);
	IDArray = vector<int> {p1->n*3,p1->n*3+1,p1->n*3+2,p2->n*3,p2->n*3+1,p2->n*3+2};
	Ue = arma::zeros(6);
	Fext = arma::zeros(6);
	Fe = arma::zeros(6);
}

void Element::addExtF(std::vector<double> &loadTrans) {
	// we are only dealing with perpendicular loading to the element for now
	// this will need to change to incorporate ext axial load
	Fext(1) += loadTrans[0];
	Fext(2) += loadTrans[1];
	Fext(4) += loadTrans[2];
	Fext(5) += loadTrans[3];
}

void Element::genKlocal () {
	double EAL = E * A / L;
	double EIL = E * I / L;
	double invL = 1 / L;
	arma::Mat<double>Kbar = {{EAL, 0, 0},
		{ 0, 4*EIL, 2*EIL},
		{ 0, 2*EIL, 4*EIL}};
	arma::Mat<double>RB = {	{-1, 0, 0, 1, 0, 0},
		{ 0, invL, 1, 0, -invL, 0},
		{ 0, invL, 0, 0, -invL, 1}};
	Klocal = RB.t() * Kbar * RB;
}

void Element::genROT () {
	double c = cos(theta);
	double s = sin(theta);
	ROT = {{ c, s, 0, 0, 0, 0},
		{-s, c, 0, 0, 0, 0},
		{ 0, 0, 1, 0, 0, 0},
		{ 0, 0, 0, c, s, 0},
		{ 0, 0, 0,-s, c, 0},
		{ 0, 0, 0, 0, 0, 1}};
}

void Element::genK () {
	// Split into separate functions if i want nonlinear capability in the future
	genKlocal();
	genROT();
	Ke = ROT.t() * Klocal * ROT;
}

void Element::calcFlocal(arma::Mat<double> &globalU) {
	for (int i = 0; i < 6; i++) {
		Ue[i] = globalU(IDArray[i]);
	}
	Fe = Klocal * ROT * Ue;
	Fe -= Fext;
}



void Structure::setSize(int nDof) {
	U.zeros(nDof);
	F.zeros(nDof);
	K.zeros(nDof, nDof);
}

void Structure::addPointLoad(int dof, double load) {
	F(dof) += load;
}

void Structure::addElementLoad(Element *e, std::vector<double> &loadVector) {
	// loadVector = {F1, M1, F2, M2}, use right hand rule
	arma::Mat<double> loadTrans = {
		loadVector[0]*-sin(e->theta), loadVector[0]*cos(e->theta), loadVector[1],
		loadVector[2]*-sin(e->theta), loadVector[2]*cos(e->theta), loadVector[3]};
	for (int i = 0; i < 6; i++) {
		F(e->IDArray[i]) += loadTrans(i);
	}
	// Don't need to transform because element load is already in local coords
	e->addExtF(loadVector);
}

void Structure::assembleK() {
	int nEle = (int)elements.size();
	for (int e = 0; e < nEle; e++) {
		elements[e].genK();
		for (int localI = 0; localI < 6; localI++) {
			for (int localJ = 0; localJ < 6; localJ++) {
				int globalI = elements[e].IDArray[localI];
				int globalJ = elements[e].IDArray[localJ];
				K(globalI, globalJ) += elements[e].Ke(localI, localJ);
			}
		}
	}
}

void Structure::constrainK() {
	// Modify Force Vector and Stiffness Matrix for given boundary condition
	Kbounded = K;
	int nDof = (int)points.size()*3;
	int nBounds = (int)bounds.size();
	for (int iBound = 0; iBound < nBounds; iBound++) {
		int dof = bounds[iBound].dof;
		for (int i = 0; i < nDof; i++) {
			F(i) -= Kbounded(dof,i)*bounds[iBound].disp;
			Kbounded(i,dof) = 0;
			Kbounded(dof,i) = 0;
		}
		F(dof) = bounds[iBound].disp;
		Kbounded(dof,dof) = 1;
	}
}

void Structure::solve() {
	constrainK();
	U = arma::solve(Kbounded, F);
	int nEle = (int)elements.size();
	for (int e = 0; e < nEle; e++) {
		elements[e].calcFlocal(U);
	}
	// Fsolved could be used to show reacitons in the future
	Fsolved = K * U;
}

void Structure::printNodeDisp() {
	int nPoints = (int)points.size();
	printf("Point Displacement (X,Y,R)\n" );
	for (int p = 0; p < nPoints; p++) {
		printf("Point %i\n",p);
		printf("(%f, %f, %f)\n",U(p*3),U(p*3+1),U(p*3+2));
	}
	printf("\n");
}

void Structure::printEleForce() {
	int nEle = (int)elements.size();
	printf("Local Element End Forces (A,V,M)\n");
	for (int e = 0; e < nEle; e++) {
		printf("Element %i\n",e);
		printf("(%f, %f, %f)\n",elements[e].Fe(0),elements[e].Fe(1),elements[e].Fe(2));
		printf("(%f, %f, %f)\n",elements[e].Fe(3),elements[e].Fe(4),elements[e].Fe(5));
	}
	printf("\n");
}

void parseRawInput(RawInput *input, Structure *s);
int dirConvert (string dir);
size_t split(const string &txt, vector<string> &strs, char ch);

void Structure::parseFile(string filename) {
	ifstream myfile (filename);
	if (!myfile.is_open())
		throw runtime_error("Could not open file.");
	
	string line;
	vector<string> data;
	RawInput input;
	int inputBlock = -1;
	
	while (getline (myfile,line)) {
		if (line == "POINTS") {
			inputBlock = 0;
			continue;
		} else if (line == "ELEMENTS") {
			inputBlock = 1;
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
			case 0:
				input.points.push_back(data);
				break;
			case 1:
				input.elements.push_back(data);
				break;
			case 2:
				input.bounds.push_back(data);
				break;
			case 3:
				input.loads.push_back(data);
				break;
		}
	}
	myfile.close();
	
	parseRawInput(&input, this);
}
	
void parseRawInput(RawInput *input, Structure *s) {
	int nPoint = (int)input->points.size();
	int nDof = nPoint*3;
	int nEle = (int)input->elements.size();
	int nBound = (int)input->bounds.size();
	int nLoad = (int)input->loads.size();
	
	s->setSize(nDof);
	vector<string> data;
	
	for (int i = 0; i < nPoint; i++) {
		data = input->points[i];
		double x = atof(data[0].c_str());
		double y = atof(data[1].c_str());
		Point tmpPoint {i, x, y};
		s->points.push_back(tmpPoint);
	}
	for (int i = 0; i < nEle; i++) {
		data = input->elements[i];
		int p1 = stoi(data[0]);
		int p2 = stoi(data[1]);
		double E = atof(data[2].c_str());
		double I = atof(data[3].c_str());
		double A = atof(data[4].c_str());
		Element tmpElement = Element(&s->points[p1], &s->points[p2], E, I, A);
		s->elements.push_back(tmpElement);
	}
	for (int i = 0; i < nBound; i++) {
		data = input->bounds[i];
		int p = stoi(data[0]);
		int dof = s->points[p].n*3 + dirConvert(data[1]);
		double disp = atof(data[2].c_str());
		Bound tmpBound {dof, disp};
		s->bounds.push_back(tmpBound);
	}
	for (int i = 0; i < nLoad; i++) {
		data = input->loads[i];
		vector<double> tmpLoad;
		if (data[0] == "POINT") {
			int p = stoi(data[1]);
			int dof = s->points[p].n*3 + dirConvert(data[2]);
			double load = atof(data[3].c_str());
			s->addPointLoad(dof, load);
		} else if (data[0] == "ELEMENT") {
			Element *tmpElement = &s->elements[stoi(data[1])];
			double L = tmpElement->L;
			// tmpLoad = {F1, M1, F2, M2}
			if (data[2] == "DISTRIBUTED") {
				double w1 = atof(data[3].c_str());
				double w2 = atof(data[4].c_str());
				tmpLoad = { (7*w1+3*w2)*L/20, (w1/2+w2/3)*pow(L,2)/10,
					(3*w1+7*w2)*L/20,-(w1/3+w2/2)*pow(L,2)/10};
			} else if (data[2] == "POINT") {
				double P = atof(data[3].c_str());
				double a = atof(data[4].c_str());
				double b = L-a;
				tmpLoad = { P*pow(b,2)*(3*a+b)/pow(L,3), P*a*pow(b,2)/pow(L,2),
					P*pow(a,2)*(3*b+a)/pow(L,3),-P*b*pow(a,2)/pow(L,2)};
			} else {
				throw runtime_error("ELEMENT type loading must be either 'DISTRIBUTED' or 'POINT'");
			}
			// Add Local Element Vector to structure
			s->addElementLoad(tmpElement, tmpLoad);
		} else {
			throw runtime_error("Can only handle Loading of type 'ELEMENT' or 'POINT'");
		}
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
