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
class FreakFace;

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
	void drawFaces(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	void drawNormals(std::shared_ptr<MatrixStack> MV, std::shared_ptr<MatrixStack> P) const;
	void drawTrajectory(std::shared_ptr<MatrixStack> MV, std::shared_ptr<MatrixStack> P);
	
	double getTime() const { return t; }

	void sendAction(bool w, bool s, bool a, bool d, bool k);
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;
	Eigen::Vector2d cannonDir;
	double power;

	Eigen::Vector3d calcVel(Eigen::Vector3d pos);
	
	std::shared_ptr<Shape> faceShape;
	std::shared_ptr<Shape> sphereShape;
	std::shared_ptr<Gelatin> gelatin;
	std::vector< std::shared_ptr<FreakFace> > faces;
};

#endif
