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

void parseFile(string filename, Structure *s);
void parseRawInput(RawInput *input, Structure *s);
int dirConvert (string dir);
size_t split(const string &txt, vector<string> &strs, char ch);

int main(int argc, char *argv[]) {
	if (argc != 2)
		throw runtime_error("Please specify input file");
	
	Structure system;
	parseFile(argv[1], &system);
	
	system.assembleK();
	system.solve();
	
	cout << "Solved displacement vector\n" << system.U << endl;
	
	int nEle = (int)system.elements.size();
	for (int e = 0; e < nEle; e++) {
		printf("Local Forces for Element %i\n",e);
		cout << system.elements[e].Fe << endl;
	}

	return 0;
}

void parseFile(string filename, Structure *s) {
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
	
	parseRawInput(&input, s);
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
		Element tmpElement {&s->points[p1], &s->points[p2], E, I, A};
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
		if (data[0] == "ELEMENT") {
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
			throw runtime_error("Can only handle Loading of type 'ELEMENT'");
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
