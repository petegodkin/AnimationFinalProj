#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "Gelatin.h"
#include "Shape.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
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
	gelatin = make_shared<Gelatin>(rows, cols, layers, x00, x01, x10, x11, mass, stiffness, damping);
	
	sphereShape = make_shared<Shape>();
	sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	auto sphere = make_shared<Particle>(sphereShape);
	//spheres.push_back(sphere);
	sphere->r = 0.1;
	sphere->x = Vector3d(0.0, 0.2, 0.0);
}

void Scene::init()
{
	sphereShape->init();
	gelatin->init();
}

void Scene::tare()
{
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->tare();
	}
	gelatin->tare();
}

void Scene::reset()
{
	t = 0.0;
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->reset();
	}
	gelatin->reset();
}

void Scene::sendAction(bool jump) {
    if (jump) {
        gelatin->jump();
    }
}

void Scene::step()
{
	t += h;
	
	// Move the sphere
	if(!spheres.empty()) {
		auto s = spheres.front();
		Vector3d x0 = s->x;
		double radius = 0.5;
		double a = 2.0*t;
		s->x(2) = radius * sin(a);
		Vector3d dx = s->x - x0;
		s->v = dx/h;
	}
	
	// Simulate the gelatin
	gelatin->step(h, grav, spheres);
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	//glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	for(int i = 0; i < (int)spheres.size(); ++i) {
		spheres[i]->draw(MV, prog);
	}
	gelatin->draw(MV, prog);
}

void Scene::drawNormals(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const
{
	gelatin->drawNormals(MV, prog, P);
}
