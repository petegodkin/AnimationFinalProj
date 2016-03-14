#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Gelatin.h"
#include "Shape.h"
#include "Program.h"
#include "FreakFace.h"
#include "GLSL.h"
#include "MatrixStack.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0),
	power(25.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 2e-3;
	
	grav << 0.0, -9.8, 0.0;
	
	int rows = 6;
	int cols = 6;
	int layers = 6;
	double mass = 0.1;
	double stiffness = 1e1;
	Vector2d damping(1.0, 1.0);
	Vector3d x00(-0.25, 0.5, 0.0);
	Vector3d x01(0.25, 0.5, 0.0);
	Vector3d x10(-0.25, 0.5, -0.5);
	Vector3d x11(0.25, 0.5, -0.5);
	gelatin = make_shared<Gelatin>(rows, cols, layers, x00, x01, x10, x11, mass, stiffness, damping, Vector3d(0, 2, 8));
	
	faceShape = make_shared<Shape>();
	faceShape->loadMesh(RESOURCE_DIR + "freak_face.obj");

	sphereShape = make_shared<Shape>();
    sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");

	
	auto freak = make_shared<FreakFace>(faceShape, sphereShape);
	faces.push_back(freak);
	freak->r = 1.0;
	freak->x = Vector3d(0.0, 2.0, 0.0);
}

void Scene::init()
{
	faceShape->init();
	sphereShape->init();
	gelatin->init();
}

void Scene::tare()
{
	for(int i = 0; i < (int)faces.size(); ++i) {
		faces[i]->tare();
	}
	gelatin->tare();
}

void Scene::reset()
{
	t = 0.0;
	for(int i = 0; i < (int)faces.size(); ++i) {
		faces[i]->reset();
	}
	gelatin->reset();
}

Vector3d Scene::calcVel(Vector3d pos) {
    return pos + Vector3d(cannonDir(0), cannonDir(1), 0 - pos(2)).normalized() * power;
}

void Scene::sendAction(bool w, bool s, bool a, bool d, bool k) {
    const double speed = 0.3;
    Vector3d vel(0, 0, 0);
    if (w) {
        cannonDir(1) += speed;
    }
    if (s) {
        cannonDir(1) -= speed;
    }
    if (a) {
        cannonDir(0) -= speed;
    }
    if (d) {
        cannonDir(0) += speed;
    }
    //gelatin->move(vel);

    if (k) {
        gelatin->shoot(calcVel(gelatin->getCenter()));
    }
}

void Scene::step()
{
	t += h;
	
	// Move the sphere
	if(!faces.empty()) {
		auto s = faces.front();
		Vector3d x0 = s->x;
		double radius = 0.5;
		double a = 2.0*t;
		s->x(0) = radius * cos(a);
		s->x(2) = radius * sin(a);
		Vector3d dx = s->x - x0;
		s->v = dx/h;
	}
	
	// Simulate the gelatin
	gelatin->step(h, grav, faces);
}

void Scene::drawFaces(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	//glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	for(int i = 0; i < (int)faces.size(); ++i) {
		faces[i]->draw(MV, prog);
	}
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	gelatin->draw(MV, prog);
}

void Scene::drawNormals(shared_ptr<MatrixStack> MV, shared_ptr<MatrixStack> P) const
{
	gelatin->drawNormals(MV, P);
}

void Scene::drawTrajectory(shared_ptr<MatrixStack> MV, shared_ptr<MatrixStack> P)
{
    GLSL::checkError(GET_FILE_LINE);
	glMatrixMode(GL_PROJECTION);
	GLSL::checkError(GET_FILE_LINE);
	glPushMatrix();
	glLoadMatrixf(P->topMatrix().data());

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadMatrixf(MV->topMatrix().data());
    glBegin(GL_LINE_STRIP);

	Vector3d cur = gelatin->getCenter();
	double mass = gelatin->getMass();
	Vector3d vel = calcVel(cur);
    const int numSteps = 500;
    for (int i = 0; i < numSteps; i++) {
        glVertex3f(cur(0), cur(1), cur(2));
        cur += vel * h;
        vel += h * mass * grav;
    }

    glEnd();

    glPopMatrix();

    // Pop projection matrix
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}