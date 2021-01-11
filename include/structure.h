#include <vector>
#include <string>
#include <armadillo>

struct RawInput {
	std::vector<std::vector <std::string>> points;
	std::vector<std::vector <std::string>> elements;
	std::vector<std::vector <std::string>> bounds;
	std::vector<std::vector <std::string>> loads;
};

struct Point {
	int n;
	double x, y;
};

struct Bound {
	int dof;
	double disp;
};

struct Element {
	Element(Point*, Point*, double, double, double);
	Point  *p1, *p2;
	double E, I, A, L, theta;
	std::vector<int> IDArray;
	arma::Mat<double> ROT;
	arma::Mat<double> Klocal;
	arma::Mat<double> Ke;
	arma::Mat<double> Ue;
	arma::Mat<double> Fext;
	arma::Mat<double> Fe;
	void addExtF(std::vector<double> &loadTrans);
	void genKlocal();
	void genROT();
	void genK();
	void calcFlocal(arma::Mat<double> &globalU);
};

struct Structure {
	void parseFile(std::string filename);
	std::vector<Point> points;
	std::vector<Element> elements;
	std::vector<Bound> bounds;
	arma::Mat<double> U;
	arma::Mat<double> F;
	arma::Mat<double> Fsolved;
	arma::Mat<double> K;
	arma::Mat<double> Kbounded;
	void setSize(int nDof);
	void addPointLoad(int dof, double load);
	void addElementLoad(Element *e, std::vector<double> &loadVector);
	void assembleK();
	void constrainK();
	void solve();
	void printNodeDisp();
	void printEleForce();
};
