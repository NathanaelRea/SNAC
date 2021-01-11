// GLEW & GLFW
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>

#include <armadillo>
#include "structure.h"

using namespace std;
using namespace arma;

// Window
const GLint WIDTH = 1000, HEIGHT = 1000;

void drawStructure(Structure *sys, double s);
int windowSetup(GLFWwindow *window, Structure *s);

int main(int argc, char *argv[]) {
	if (argc != 2)
		throw runtime_error("Please specify input file");
	
	// Parse Sysem Geometry from input file
	Structure system;
	system.parseFile(argv[1]);
	
	// Init window and viewport
	glfwInit();
	GLFWwindow *window = glfwCreateWindow(WIDTH, HEIGHT, "SNAC", nullptr, nullptr);
	windowSetup(window, &system);
	
	// Solve System
	system.assembleK();
	system.solve();
	system.printNodeDisp();
	system.printEleForce();
	
	// Deformed display scale
	double scale = 10;
	
	while (!glfwWindowShouldClose(window)) {
		// Input
		glfwPollEvents();
		
		// Render
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		
		// Draw
		drawStructure(&system, scale);
		
		// Swap screen buffers
		glfwSwapBuffers(window);
	}
	
	// De-allocate
	glfwTerminate();
	
	return EXIT_SUCCESS;
}

void drawStructure(Structure *sys, double s) {
	// Probably bad to call draw function e times each frame;
	// Undeformed
	glEnable(GL_LINE_STIPPLE);
	for (auto e = sys->elements.begin(); e != sys->elements.end(); ++e) {
	    glBegin(GL_LINES);
		glVertex3f(e->p1->x, e->p1->y, 0.0f);
		glVertex3f(e->p2->x, e->p2->y, 0.0f);
	    glEnd();
	}
	// Deformed
	glDisable(GL_LINE_STIPPLE);
	for (auto e = sys->elements.begin(); e != sys->elements.end(); ++e) {
	    glBegin(GL_LINES);
		glVertex3f(e->p1->x + e->Ue[0]*s, e->p1->y + e->Ue[1]*s, 0.0f);
		glVertex3f(e->p2->x + e->Ue[3]*s, e->p2->y + e->Ue[4]*s, 0.0f);
	    glEnd();
	}
}

int windowSetup(GLFWwindow *window, Structure *s) {
	// Window
	int screenWidth, screenHeight;
	glfwGetFramebufferSize(window, &screenWidth, &screenHeight);
	if (window == nullptr) {
		cout << "Failed to create GLFW window" << endl;
		glfwTerminate();
		return EXIT_FAILURE;
	}
	glfwMakeContextCurrent(window);
	
	// Viewport
	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK) {
		cout << "Failed to init GLEW" << endl;
		return EXIT_FAILURE;
	}
	glViewport(0, 0, WIDTH, HEIGHT);
	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glDepthMask(false);
	glLineWidth(3);
	glLineStipple(3, 0xAAAA);
	
	// Setup camera position
	double minx = 0;
	double maxx = 0;
	double miny = 0;
	double maxy = 0;
	for (auto &it : s->points) {
		minx = min(minx, it.x);
		maxx = max(maxx, it.x);
		miny = min(miny, it.y);
		maxy = max(maxy, it.y);
	}
	double cx = (minx + maxx)/2;
	double cy = (miny + maxy)/2;
	// 10% padding
	double box = 1.1*max(maxx-minx, maxy-miny)/2;
	
	GLdouble left = cx - box;
	GLdouble right = cx + box;
	GLdouble bottom = cy - box;
	GLdouble top = cy + box;
	
	glOrtho(left, right, bottom, top, 0.0, 1.0);
	
	return 0;
}
