#pragma once
#ifndef __FreakFace__
#define __FreakFace__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Shape;
class Program;
class MatrixStack;

class FreakFace
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	FreakFace();
	FreakFace(const std::shared_ptr<Shape> shape);
	virtual ~FreakFace();
	void tare();
	void reset();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p) const;

	double r; // radius
	double m; // mass
	int i;  // starting index
	Eigen::Vector3d x0; // initial position
	Eigen::Vector3d v0; // initial velocity
	Eigen::Vector3d x;  // position
	Eigen::Vector3d v;  // velocity
	bool fixed;

private:
    const double objScale = 0.1;

	const std::shared_ptr<Shape> shape;
};

#endif
