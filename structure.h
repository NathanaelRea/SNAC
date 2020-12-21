#include<vector>

using namespace std;
using namespace arma;

struct Point {
	int n;
	double x, y;
	double load = 0.0;
};

struct Element {
	Point  *p1, *p2;
	double E, I, A;
	double L = sqrt(pow(p1->x - p2->x,2) + pow(p1->y - p2->y,2));
	double theta = atan2(p2->y - p1->y, p2->x - p1->x);
	vector<int> IDArray {p1->n*3,p1->n*3+1,p1->n*3+2,p2->n*3,p2->n*3+1,p2->n*3+2};
};

struct Bound {
	int dof;
	double disp;
};

struct Structure {
	vector<Point> points;
	vector<Element> elements;
	vector<Bound> bounds;
	Mat<double> globalLoad;
	// loadVector = {F1, M1, F2, M2}
	// Follow right hand rule for 
	void addElementLoad(Element *e, vector<double> loadVector) {
		this->globalLoad[e->IDArray[0]] += loadVector[0]*-sin(e->theta);
		this->globalLoad[e->IDArray[1]] += loadVector[0]*cos(e->theta);
		this->globalLoad[e->IDArray[2]] += loadVector[1];
		this->globalLoad[e->IDArray[3]] += loadVector[2]*-sin(e->theta);
		this->globalLoad[e->IDArray[4]] += loadVector[2]*cos(e->theta);
		this->globalLoad[e->IDArray[5]] += loadVector[3];
	}
};
