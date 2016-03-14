#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"
#include "FreakFace.h"

using namespace std;

FreakFace::FreakFace() :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(true)
{

}

FreakFace::FreakFace(const shared_ptr<Shape> s) :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(true),
	shape(s)
{

}

FreakFace::~FreakFace()
{
}

void FreakFace::tare()
{
	x0 = x;
	v0 = v;
}

void FreakFace::reset()
{
	x = x0;
	v = v0;
}

void FreakFace::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if(shape) {
		MV->pushMatrix();
		MV->translate(Eigen::Vector3f(x(0), x(1), x(2)));
		MV->scale(r * objScale);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, MV->topMatrix().data());
		shape->draw(prog);
		MV->popMatrix();
	}
}
