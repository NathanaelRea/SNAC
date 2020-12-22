#include<vector>

struct Point {
	int n;
	double x, y;
};

struct Bound {
	int dof;
	double disp;
};

arma::Mat<double> makeLocalK(double E, double I, double A, double L, double theta) {
	// 1D Local Stiffness Matrix
	arma::Mat<double>localK = {	{E*A/L,       0,       0},
								{    0, 4*E*I/L, 2*E*I/L},
								{    0, 2*E*I/L, 4*E*I/L}};
	// Rigid Body Matrix
	arma::Mat<double>RB = {	{-1,   0, 0, 1,    0, 0},
							{ 0, 1/L, 1, 0, -1/L, 0},
							{ 0, 1/L, 0, 0, -1/L, 1}};
	// Rotational Matrix
    double c = cos(theta);
    double s = sin(theta);
	arma::Mat<double>ROT = {{ c, s, 0, 0, 0, 0},
							{-s, c, 0, 0, 0, 0},
							{ 0, 0, 1, 0, 0, 0},
							{ 0, 0, 0, c, s, 0},
							{ 0, 0, 0,-s, c, 0},
							{ 0, 0, 0, 0, 0, 1}};
	// Combine into Global Domain
	arma::Mat<double> kBar = ROT.t() * RB.t() * localK * RB * ROT;
	return kBar;
}

struct Element {
	Point  *p1, *p2;
	double E, I, A;
	double L = sqrt(pow(p1->x - p2->x,2) + pow(p1->y - p2->y,2));
	double theta = atan2(p2->y - p1->y, p2->x - p1->x);
	std::vector<int> IDArray {p1->n*3,p1->n*3+1,p1->n*3+2,p2->n*3,p2->n*3+1,p2->n*3+2};
	arma::Mat<double> K = makeLocalK(E, I, A, L, theta);
};

struct Structure {
	int nDof = 0;
	int nEle = 0;
	int nBounds = 0;
	std::vector<Point> points;
	std::vector<Element> elements;
	std::vector<Bound> bounds;
	arma::Mat<double> u;
	arma::Mat<double> F;
	arma::Mat<double> Fsolved;
	arma::Mat<double> K;
	arma::Mat<double> Kbounded;
	void setSize() {
		this->nDof = (int)this->points.size()*3;
		this->nEle  = (int)this->elements.size();
		this->nBounds = (int)this->bounds.size();
		u.zeros(this->nDof);
		F.zeros(this->nDof);
		K.zeros(this->nDof, this->nDof);
	}
	void addElementLoad(Element *e, std::vector<double> loadVector) {
		// loadVector = {F1, M1, F2, M2}, use right hand rule
		this->F[e->IDArray[0]] += loadVector[0]*-sin(e->theta);
		this->F[e->IDArray[1]] += loadVector[0]*cos(e->theta);
		this->F[e->IDArray[2]] += loadVector[1];
		this->F[e->IDArray[3]] += loadVector[2]*-sin(e->theta);
		this->F[e->IDArray[4]] += loadVector[2]*cos(e->theta);
		this->F[e->IDArray[5]] += loadVector[3];
	}
	void genGlobalK() {
		for (int e = 0; e < this->nEle; e++) {
			for (int localI = 0; localI < 6; localI++) {
				for (int localJ = 0; localJ < 6; localJ++) {
					int globalI = this->elements[e].IDArray[localI];
					int globalJ = this->elements[e].IDArray[localJ];
					this->K(globalI, globalJ) += this->elements[e].K(localI, localJ);
				}
			}
		}
	}
	void constrainK() {
		// Modify Force Vector and Stiffness Matrix for given boundary condition
		this->Kbounded = this->K;
		for (int iBound = 0; iBound < this->nBounds; iBound++) {
			int dof = this->bounds[iBound].dof;
			for (int i = 0; i < nDof; i++) {
				this->F(i) -= this->Kbounded(dof,i)*this->bounds[iBound].disp;
				this->Kbounded(i,dof) = 0;
				this->Kbounded(dof,i) = 0;
			}
			this->F(dof) = this->bounds[iBound].disp;
			this->Kbounded(dof,dof) = 1;
		}
	}
	void solve() {
		this->constrainK();
		this->u = arma::solve(this->Kbounded, this->F);
	}
};
