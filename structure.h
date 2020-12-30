#include<vector>

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
	Point  *p1, *p2;
	double E, I, A;
	double L = sqrt(pow(p1->x - p2->x,2) + pow(p1->y - p2->y,2));
	double theta = atan2(p2->y - p1->y, p2->x - p1->x);
	std::vector<int> IDArray {p1->n*3,p1->n*3+1,p1->n*3+2,p2->n*3,p2->n*3+1,p2->n*3+2};
	arma::Mat<double> ROT = {};
	arma::Mat<double> Klocal = {};
	arma::Mat<double> Ke = {};
	arma::Mat<double> Ue = arma::zeros(6);
	arma::Mat<double> Fext = arma::zeros(6);
	arma::Mat<double> Fe = arma::zeros(6);
	void addExtF(std::vector<double> &loadTrans) {
		// we are only dealing with perpendicular loading to the element for now
		// this will need to change to incorporate ext axial load
		Fext(1) += loadTrans[0];
		Fext(2) += loadTrans[1];
		Fext(4) += loadTrans[2];
		Fext(5) += loadTrans[3];
	}
	void genKlocal () {
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
	void genROT () {
		double c = cos(theta);
		double s = sin(theta);
		ROT = {{ c, s, 0, 0, 0, 0},
			{-s, c, 0, 0, 0, 0},
			{ 0, 0, 1, 0, 0, 0},
			{ 0, 0, 0, c, s, 0},
			{ 0, 0, 0,-s, c, 0},
			{ 0, 0, 0, 0, 0, 1}};
	}
	void genK () {
		// Split into separate functions if i want nonlinear capability in the future
		genKlocal();
		genROT();
		Ke = ROT.t() * Klocal * ROT;
	}
	void calcFlocal(arma::Mat<double> &globalU) {
		for (int i = 0; i < 6; i++) {
			Ue[i] = globalU(IDArray[i]);
		}
		Fe = Klocal * ROT * Ue;
		Fe -= Fext;
	}
};

struct Structure {
	std::vector<Point> points;
	std::vector<Element> elements;
	std::vector<Bound> bounds;
	arma::Mat<double> U;
	arma::Mat<double> F;
	arma::Mat<double> Fsolved;
	arma::Mat<double> K;
	arma::Mat<double> Kbounded;
	void setSize(int nDof) {
		U.zeros(nDof);
		F.zeros(nDof);
		K.zeros(nDof, nDof);
	}
	void addPointLoad(int dof, double load) {
		F(dof) += load;
	}
	void addElementLoad(Element *e, std::vector<double> &loadVector) {
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
	void assembleK() {
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
	void constrainK() {
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
	void solve() {
		constrainK();
		U = arma::solve(Kbounded, F);
		int nEle = (int)elements.size();
		for (int e = 0; e < nEle; e++) {
			elements[e].calcFlocal(U);
		}
		// Fsolved could be used to show reacitons in the future
		Fsolved = K * U;
	}
	void printNodeDisp() {
		int nPoints = (int)points.size();
		printf("Point Displacement (X,Y,R)\n" );
		for (int p = 0; p < nPoints; p++) {
			printf("Point %i\n",p);
			printf("(%f, %f, %f)\n",U(p*3),U(p*3+1),U(p*3+2));
		}
		printf("\n");
	}
	void printEleForce() {
		int nEle = (int)elements.size();
		printf("Local Element End Forces (A,V,M)\n");
		for (int e = 0; e < nEle; e++) {
			printf("Element %i\n",e);
			printf("(%f, %f, %f)\n",elements[e].Fe(0),elements[e].Fe(1),elements[e].Fe(2));
			printf("(%f, %f, %f)\n",elements[e].Fe(3),elements[e].Fe(4),elements[e].Fe(5));
		}
		printf("\n");
	}
};
