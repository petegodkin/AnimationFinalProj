#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Gelatin;
class Particle;
class MatrixStack;
class Program;
class Shape;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	void drawNormals(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P) const;
	
	double getTime() const { return t; }

	void sendAction(bool jump);
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;
	
	std::shared_ptr<Shape> sphereShape;
	std::shared_ptr<Gelatin> gelatin;
	std::vector< std::shared_ptr<Particle> > spheres;
};

#endif
