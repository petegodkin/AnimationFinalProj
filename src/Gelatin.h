#pragma once
#ifndef __Gelatin__
#define __Gelatin__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class Spring;
class MatrixStack;
class Program;

class Gelatin
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Gelatin(int rows, int cols, int layers,
		  const Eigen::Vector3d &x00,
		  const Eigen::Vector3d &x01,
		  const Eigen::Vector3d &x10,
		  const Eigen::Vector3d &x11,
		  double mass,
		  double stiffness,
		  const Eigen::Vector2d &damping);
	virtual ~Gelatin();

	void tare();
	void reset();
	void updatePosNor();
	void step(double h, const Eigen::Vector3d &grav, const std::vector< std::shared_ptr<Particle> > spheres);

	void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prgm) const;
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prgm, std::shared_ptr<MatrixStack> P) const;

private:
	int rows;
	int cols;
	int layers;
	int n;
	Eigen::Vector2d damping;
	std::vector< std::shared_ptr<Particle> > particles;
	std::vector< std::shared_ptr<Spring> > springs;

	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::MatrixXd M;
	Eigen::MatrixXd K;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

	int setNormals(int curNorIndex, int index, int adj[4]);
};

#endif