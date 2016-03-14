#include <iostream>

#include "Gelatin.h"
#include "Particle.h"
#include "Spring.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"

using namespace std;
using namespace Eigen;

shared_ptr<Spring> createSpring(const shared_ptr<Particle> p0, const shared_ptr<Particle> p1, double E)
{
	auto s = make_shared<Spring>(p0, p1);
	s->E = E;
	Vector3d x0 = p0->x;
	Vector3d x1 = p1->x;
	Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}

Gelatin::Gelatin(int rows, int cols, int layers,
			 const Vector3d &x00,
			 const Vector3d &x01,
			 const Vector3d &x10,
			 const Vector3d &x11,
			 double mass,
			 double stiffness,
			 const Vector2d &damping)
{
	assert(rows > 1);
	assert(cols > 1);
	assert(layers > 1);
	assert(mass > 0.0);
	assert(stiffness > 0.0);

	this->rows = rows;
	this->cols = cols;
	this->layers = layers;
	this->damping = damping;

	// Create particles
	n = 0;
	double r = 0.02; // Used for collisions
	int nVerts = rows*cols*layers;


	for(int i = 0; i < rows; ++i) {
		double u = i / (rows - 1.0);
		Vector3d x0 = (1 - u)*x00 + u*x10;
		Vector3d x1 = (1 - u)*x01 + u*x11;
		for(int j = 0; j < cols; ++j) {
			double v = j / (cols - 1.0);
			for (int k = 0; k < layers; ++k) {
			    Vector3d x = (1 - v)*x0 + v*x1;
			    double w = k / (layers - 1.0) * 0.5;
			    x(1) += w;
			    //std::cout << "X[" << i * rows * cols + j * cols + k << "] = " << x(0) << ", " << x(1) << ", " << x(2) << std::endl;

                auto p = make_shared<Particle>();
                particles.push_back(p);
                p->r = r;
                p->x = x;
                p->v << 0.0, 0.0, 0.0;
                p->m = mass/(nVerts);
                // Pin two particles
                if(i == 0 && (j == 0) && k == 0) {
                    p->fixed = true;
                    p->i = -1;
                } else {
                    p->fixed = false;
                    p->i = n;
                    n += 3;
                }
			}
		}
	}

	for (int k = 0; k < layers; ++k) {

        // Create x springs
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols-1; ++j) {
                int k0 = i*cols + j + k*cols*rows;
                int k1 = k0 + 1;
                springs.push_back(createSpring(particles[k0], particles[k1], stiffness));
            }
        }

        // Create y springs
        for(int j = 0; j < cols; ++j) {
            for(int i = 0; i < rows-1; ++i) {
                int k0 = i*cols + j + k*cols*rows;
                int k1 = k0 + cols;
                springs.push_back(createSpring(particles[k0], particles[k1], stiffness));
            }
        }

        // Create shear springs
        for(int i = 0; i < rows-1; ++i) {
            for(int j = 0; j < cols-1; ++j) {
                int k00 = i*cols + j + k*cols*rows;
                int k10 = k00 + 1;
                int k01 = k00 + cols;
                int k11 = k01 + 1;
                springs.push_back(createSpring(particles[k00], particles[k11], stiffness));
                springs.push_back(createSpring(particles[k10], particles[k01], stiffness));
            }
        }
	}

	for (int k = 0; k < layers-1; ++k) {
	    for(int i = 0; i < rows; ++i) {
	        for(int j = 0; j < cols; ++j) {
	            int k0 = i*cols + j + k*cols*rows;
                int k1 = k0 + cols * rows;
                springs.push_back(createSpring(particles[k0], particles[k1], stiffness));
	        }
	    }

        // Create shear springs
        for(int i = 0; i < rows-1; ++i) {
            for(int j = 0; j < cols-1; ++j) {
                int k00 = i*cols + j + k*cols*rows;
                int k10 = k00 + cols * rows;
                int k01 = k00 + 1;
                int k11 = k01 + cols * rows;
                springs.push_back(createSpring(particles[k00], particles[k11], stiffness));
                springs.push_back(createSpring(particles[k10], particles[k01], stiffness));
            }
        }

        // Create shear springs
        for(int i = 0; i < rows-1; ++i) {
            for(int j = 0; j < cols-1; ++j) {
                int k00 = i*cols + j + k*cols*rows;
                int k10 = k00 + cols * rows;
                int k01 = k00 + cols;
                int k11 = k01 + cols * rows;
                springs.push_back(createSpring(particles[k00], particles[k11], stiffness));
                springs.push_back(createSpring(particles[k10], particles[k01], stiffness));
            }
        }
	}

	// Build system matrices and vectors
	M.resize(n,n);
	K.resize(n,n);
	v.resize(n);
	f.resize(n);

	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();
	posBuf.resize(cols*layers*3*6);//nVerts*3);
	norBuf.resize(cols*layers*3*6);//nVerts*3);
	updatePosNor();

	for (int k = 0; k < 6; ++k) {
	    int initial = k * cols * rows;
	    for(int i = 0; i < cols - 1; ++i) {
            for(int j = 0; j < layers - 1; ++j) {
                int k0 = initial + i * layers + j;
                int k1 = initial + i * layers + j + 1;
                int k2 = initial + (i+1) * layers + j + 1;
                int k3 = initial + (i+1) * layers + j;
                eleBuf.push_back(k0);
                eleBuf.push_back(k1);
                eleBuf.push_back(k2);
                eleBuf.push_back(k3);
            }
        }
    }
}



Gelatin::~Gelatin()
{
}

void Gelatin::tare()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->tare();
	}
}

void Gelatin::reset()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->reset();
	}
	updatePosNor();
}

int Gelatin::setNormals(int curNorIndex, int index, int adj[4]) {
    Vector3d x = particles[index]->x;
    Vector3d nor(0.0, 0.0, 0.0);

    for (int k = 0; k < 3; k++) {
        if (adj[k] != -1 && adj[k+1] != -1) {
            Vector3d first = (particles[adj[k]]->x - x).normalized();
            Vector3d second = (particles[adj[k+1]]->x - x).normalized();

            Vector3d res = first.cross(second).normalized();
            nor += res;
        }
    }

    if (adj[3] != -1 && adj[0] != -1) {
        Vector3d first = (particles[adj[3]]->x - x).normalized();
        Vector3d second = (particles[adj[0]]->x - x).normalized();

        Vector3d res = first.cross(second).normalized();
        nor += res;
    }

    nor = nor.normalized();

    norBuf[curNorIndex++] = nor(0);
    norBuf[curNorIndex++] = nor(1);
    norBuf[curNorIndex++] = nor(2);

    return curNorIndex;
}

int frontIndex(int cols, int rows, int layers, int i, int j) {
    return i*layers + j;
}

int backIndex(int cols, int rows, int layers, int i, int j) {
    int initial = cols * layers * (rows - 1);
    return initial + i*layers + j;
}

int bottomIndex(int cols, int rows, int layers, int i, int j) {
    return i*cols + j * cols * layers;
}

int topIndex(int cols, int rows, int layers, int i, int j) {
    int initial = cols - 1;
    return initial + i*cols + j * cols * layers;
}

int leftIndex(int cols, int rows, int layers, int i, int j) {
    return i + j * rows * layers;
}

int rightIndex(int cols, int rows, int layers, int i, int j) {
    int initial = rows * (layers - 1);
    return initial + i + j * cols * layers;
}

void Gelatin::setAdj(int* adj, int i, int j, int imax, int jmax, int (*calcIndex)(int, int, int, int, int)) {
    if (i > 0) {
        adj[0] = calcIndex(cols, rows, layers, i-1, j);
    }
    if (j > 0) {
        adj[1] = calcIndex(cols, rows, layers, i, j-1);
    }
    if (i < imax) {
        adj[2] = calcIndex(cols, rows, layers, i+1, j);
    }
    if (j < jmax) {
        adj[3] = calcIndex(cols, rows, layers, i, j+1);
    }
}

void Gelatin::setAdjRev(int* adj, int i, int j, int imax, int jmax, int (*calcIndex)(int, int, int, int, int)) {
    if (i > 0) {
        adj[3] = calcIndex(cols, rows, layers, i-1, j);
    }
    if (j > 0) {
        adj[2] = calcIndex(cols, rows, layers, i, j-1);
    }
    if (i < imax) {
        adj[1] = calcIndex(cols, rows, layers, i+1, j);
    }
    if (j < jmax) {
        adj[0] = calcIndex(cols, rows, layers, i, j+1);
    }
}

void Gelatin::updatePosNor()
{
    int curIndex = 0;
    int curNorIndex = 0;
	// Position row == 0, front
    for(int i = 0; i < cols; ++i) {
        for (int j = 0; j < layers; ++j) {
            int index = frontIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);

            int adj[4] = {-1, -1, -1, -1};
            setAdj(adj, i, j, cols-1, layers-1, frontIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
        }
    }

    //row == rows - 1, back
    for(int i = 0; i < cols; ++i) {
        for (int j = 0; j < layers; ++j) {
            int index = backIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);

            int adj[4] = {-1, -1, -1, -1};
            setAdjRev(adj, i, j, cols-1, layers-1, backIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
        }
    }

    //layer == 0, bottom
    for(int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            int index = bottomIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);

            int adj[4] = {-1, -1, -1, -1};
            setAdjRev(adj, i, j, cols-1, rows-1, bottomIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
        }
    }

    //layer == layers - 1, top
    for(int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            int index =  topIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);

            int adj[4] = {-1, -1, -1, -1};
            setAdj(adj, i, j, cols-1, rows-1, topIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
        }
    }

    //col == 0, left
     for(int i = 0; i < layers; ++i) {
         for (int j = 0; j < rows; ++j) {
            int index =  leftIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);
            int adj[4] = {-1, -1, -1, -1};
            setAdj(adj, i, j, cols-1, rows-1, leftIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
         }
     }

     //col == cols - 1, right?
      for(int i = 0; i < cols; ++i) {
          for (int j = 0; j < rows; ++j) {
            int index =  rightIndex(cols, rows, layers, i, j);
            Vector3d x = particles[index]->x;
            posBuf[curIndex++] = x(0);
            posBuf[curIndex++] = x(1);
            posBuf[curIndex++] = x(2);
            int adj[4] = {-1, -1, -1, -1};
            setAdjRev(adj, i, j, cols-1, rows-1, rightIndex);
            curNorIndex = setNormals(curNorIndex, index, adj);
          }
      }
}

void Gelatin::step(double h, const Vector3d &grav, const vector< shared_ptr<Particle> > spheres)
{
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();

//
//	typedef Eigen::Triplet<double> T;
//    std::vector<T> aTrips;
//    aTrips.push_back(T(0, 0, 1));
//    aTrips.push_back(T(1, 1, 1));
//    Eigen::SparseMatrix<double> Asd(2,2);
//    Asd.setFromTriplets(aTrips.begin(), aTrips.end());

	for (int i = 0; i < particles.size(); i++) {
	    if (!particles[i]->fixed) {
	        int index = particles[i]->i;
            M.block<3,3>(index, index) = particles[i]->m * Matrix3d::Identity();
            v.block<3,1>(index, 0) = particles[i]->v;
            f.block<3,1>(index, 0) = particles[i]->m * grav;
        }
	}

	for (int i = 0; i < springs.size(); i++) {
	    double E = springs[i]->E;
	    double L = springs[i]->L;
	    auto p0 = springs[i]->p0;
	    auto p1 = springs[i]->p1;

	    //MatrixXd delta(3, 1);
        Vector3d delta = p1->x - p0->x;
	    double l = delta.norm();

	    Vector3d fs = E * (l - L) * (delta / l);

	    if (!springs[i]->p0->fixed) {
	        f.block<3,1>(p0->i, 0) += fs;
	    }
	    if (!springs[i]->p1->fixed) {
	        f.block<3,1>(p1->i, 0) -= fs;
	    }

        if (!springs[i]->p0->fixed && !springs[i]->p1->fixed) {
            Matrix3d ks =  (E / (l * l)) * ((1.0 - (l - L)/l) * (delta * delta.transpose()) + ((l - L)/l) * double(delta.transpose() * delta) * Matrix3d::Identity());
            int i0 = p0->i;
            int i1 = p1->i;
            K.block<3, 3>(i0, i0) += ks;
            K.block<3, 3>(i0, i1) -= ks;
            K.block<3, 3>(i1, i0) -= ks;
            K.block<3, 3>(i1, i1) += ks;
        }
	}

	MatrixXd b = M * v + h * f;

    MatrixXd D = h * damping(0) * M + h * h * damping(1) * K;
	MatrixXd A = M + D;

	v = A.ldlt().solve(b);

	for (int i = 0; i < particles.size(); i++) {
	    if (!particles[i]->fixed) {
            int index = particles[i]->i;
            particles[i]->v = v.block<3,1>(index, 0);
            particles[i]->x += h * particles[i]->v;
        }
	}

	for (int i = 0; i < particles.size(); i++) {
	    for (int j = 0; j < spheres.size(); j++) {
	        Vector3d delta = particles[i]->x - spheres[j]->x;
	        double dist = delta.norm();
	        if (dist < spheres[j]->r + particles[i]->r) {
	            particles[i]->v = (particles[i]->v - (particles[i]->v.dot(delta.normalized())) * delta.normalized());
	            particles[i]->x = (spheres[j]->x + (spheres[j]->r + particles[i]->r) * delta.normalized());
	        }
	    }
	}

	// Update position and normal buffers
	updatePosNor();
}

void Gelatin::init()
{
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size()*sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	assert(glGetError() == GL_NO_ERROR);
}


void Gelatin::drawNormals(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prgm, shared_ptr<MatrixStack> P) const
{

    GLSL::checkError(GET_FILE_LINE);
	glMatrixMode(GL_PROJECTION);
	GLSL::checkError(GET_FILE_LINE);
	glPushMatrix();
	glLoadMatrixf(P->topMatrix().data());

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadMatrixf(MV->topMatrix().data());
    glBegin(GL_LINES);
    for (int i = 0; i < posBuf.size(); i+=3) {
        Vector3f pos(posBuf[i], posBuf[i+1], posBuf[i+2]);
        Vector3f nor(norBuf[i], norBuf[i+1], norBuf[i+2]);
        Vector3f norTip = pos + nor * .05;
        glVertex3f(pos(0), pos(1), pos(2));
        glVertex3f(norTip(0), norTip(1), norTip(2));
    }

    glEnd();

    glPopMatrix();

    // Pop projection matrix
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
}


void Gelatin::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prgm) const
{
	MV->pushMatrix();

	glUniformMatrix4fv(prgm->getUniform("MV"), 1, GL_FALSE, MV->topMatrix().data());
	int h_pos = prgm->getAttribute("vertPos");
	GLSL::enableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	int h_nor = prgm->getAttribute("vertNor");
	GLSL::enableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size()*sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
    //	int h_tex = prgm->getAttribute("vertTex");
    //	GLSL::enableVertexAttribArray(h_tex);
    //	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
    //	glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glDrawElements(GL_QUADS, eleBuf.size(), GL_UNSIGNED_INT, (const void *)0);
	//GLSL::disableVertexAttribArray(h_tex);
	GLSL::disableVertexAttribArray(h_nor);
	GLSL::disableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

