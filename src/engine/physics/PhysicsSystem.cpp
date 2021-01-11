#include "PhysicsSystem.h"
#include "src/engine/graphics/Camera.h"
#include <numeric>

PhysicsSystem::PhysicsSystem() : colliding(false), collision_spot(glm::vec3(0)), collision_happened(false)
{
    //m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, false, false, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies.push_back(std::make_shared<EllipsoidRigidBodyComponent>(glm::vec3(1,3,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    //m_rigid_bodies.push_back(std::make_shared<CylinderRigidBodyComponent>(1, 1, 1, false, false, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[0]->setPosition(glm::vec3(0,2.5,0));
    m_rigid_bodies[0]->setEuelerAngles(glm::vec3(0,0,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, false, false, false, glm::vec3(0,.7,.7)));

    m_graphics = Graphics::getGlobalInstance();
}

void PhysicsSystem::tick(float seconds) {
    PhysicsSolver::eulerStep(this, seconds);
}

void PhysicsSystem::draw() {
    for (size_t i = 0; i < m_rigid_bodies.size(); i++) {
        //std::cout << glm::to_string(m_rigid_bodies[i]->getPosition()) << std::endl;
        //std::cout << glm::to_string(m_rigid_bodies[i]->getOrientationMatrix()) << std::endl;
        if (colliding) {
            Material m;
            m.color = glm::vec3(1,0,0);
            m_graphics->setMaterial(m);
            std::shared_ptr<Camera> ui_cam = std::make_shared<Camera>(glm::vec2(1000, 1000));
            ui_cam->setUI(true);
            auto curr_cam = m_graphics->getActiveCamera();
            m_graphics->setCamera(ui_cam);
            m_graphics->drawText("COLLISION!",100);
            m_graphics->setCamera(curr_cam);
        }
        if (collision_happened) {
            m_graphics->clearTransform();
            m_graphics->translate(collision_spot);
            m_graphics->scale(.5);
            m_graphics->drawShape("sphere");
        }
        m_rigid_bodies[i]->draw();
    }
}

void PhysicsSystem::resolveCollisions() {

    for (int i = 0; i < m_rigid_bodies.size(); i++) {
        glm::vec3 i_pos = m_rigid_bodies[i]->getPosition();
        float i_radius = m_rigid_bodies[i]->getMaxRadius();
        for (int j = i+1; j < m_rigid_bodies.size(); j++) {
            glm::vec3 j_pos = m_rigid_bodies[j]->getPosition();
            float j_radius = m_rigid_bodies[j]->getMaxRadius();
            if (glm::length(i_pos - j_pos) < i_radius + j_radius) {
                bool had_collision;
                glm::vec3 mtv;
                glm::vec3 rb1_p;
                glm::vec3 rb2_p;
                if (!collision_happened) {
                    std::tie(had_collision, mtv, rb1_p, rb2_p) = runGJKAndExpandingPolytope(m_rigid_bodies[i], m_rigid_bodies[j]);
                }
                if (had_collision) {
                    if (!collision_happened) {
                        collision_spot = rb1_p;
                    }
                    collision_happened = true;
                    std::cout << glm::to_string(mtv) << std::endl;
                    if (m_rigid_bodies[i]->getMovable()) {
                        m_rigid_bodies[i]->addToLinearMomentum(m_rigid_bodies[i]->getOrientationMatrix() * mtv);
                    }
                    if (m_rigid_bodies[j]->getMovable()) {
                        m_rigid_bodies[j]->addToLinearMomentum(m_rigid_bodies[i]->getOrientationMatrix() * -mtv);
                    }
                }
            } else {
                colliding = false;
            }
        }
    }

}

MinkowskiDifferenceResult PhysicsSystem::minkowskiDifferenceSupport(glm::vec3 dir, const std::shared_ptr<const RigidBodyComponent> rb1, const std::shared_ptr<const RigidBodyComponent> rb2) {
    // dir is in gobal space, so we need to compute dir in the object space of rb1 and rb2.
    // scaling is taken care of by the rigid body's support function, so we just need to
    // multiply by the inverse of the orientation matrix
    glm::vec3 rb1_dir = glm::inverse(rb1->getOrientationMatrix()) * dir;
    glm::vec3 rb2_dir = glm::inverse(rb2->getOrientationMatrix()) * dir;

    // support points in object space of rb1 and rb2
    glm::vec3 rb1_p_object = rb1->support(rb1_dir);
    glm::vec3 rb2_p_object = rb2->support(-rb2_dir);  // we have a negative here because we are using an identity to get the minkowski difference

    glm::vec3 rb1_p = rb1->getOrientationMatrix() * rb1_p_object + rb1->getPosition();
    glm::vec3 rb2_p = rb2->getOrientationMatrix() * rb2_p_object + rb2->getPosition();

    MinkowskiDifferenceResult ret;
    ret.m = rb1_p - rb2_p;
    ret.rb1 = rb1_p;
    ret.rb2 = rb2_p;
    return ret;
}

std::pair<bool, std::vector<MinkowskiDifferenceResult>> PhysicsSystem::runGJK(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2) {

    glm::vec3 initial_dir = glm::vec3(1,1,1);
    MinkowskiDifferenceResult s = minkowskiDifferenceSupport(initial_dir, rb1, rb2);
    std::vector<MinkowskiDifferenceResult> simplex;
    simplex.push_back(s);
    glm::vec3 d = -s.m;

    while (true) {
        s = minkowskiDifferenceSupport(d, rb1, rb2);
        if (std::isnan(s.m.x) || std::isnan(s.m.y) || std::isnan(s.m.z))
        if (glm::dot(s.m, d) < 0) {
            colliding = false;
            return std::pair<bool, std::vector<MinkowskiDifferenceResult>>(false, simplex);
        }
        simplex.push_back(s);
        std::pair<bool, glm::vec3> p = doSimplex(simplex);
        if (p.first) {
            colliding = true;
            return std::pair<bool, std::vector<MinkowskiDifferenceResult>>(true, simplex);
        }
        d = p.second;
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::doSimplex(std::vector<MinkowskiDifferenceResult> &simplex) {
    assert(simplex.size() >= 1 && simplex.size() <= 4);

    /*
    if (simplex.size() == 1 && glm::length(*simplex.begin()) < EPSILON) {
        std::cout << "finished 0" << std::endl;
        return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
    else if (simplex.size() == 2 && checkIfLineContainsOrigin(simplex)) {
        std::cout << "finished 1" << std::endl;
        //return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
    else if (simplex.size() == 3 && checkIfTriangleContainsOrigin(simplex)) {
        std::cout << "finished 2" << std::endl;
        //return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
    */
    if (simplex.size() == 4 && checkifTetrahedronContainsOrigin(simplex)) {
        //std::cout << "finished 3" << std::endl;
        return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }

    if (simplex.size() == 1) {
        glm::vec3 A = (*simplex.begin()).m;
        return std::pair<bool, glm::vec3>(false, -A);
    } else if (simplex.size() == 2) {
        //std::cout << "1" << std::endl;
        return handleSimplex1(simplex);
    } else if (simplex.size() == 3) {
        //std::cout << "2" << std::endl;
        return handleSimplex2(simplex);
    } else {
        //std::cout << "3" << std::endl;
        return handleSimplex3(simplex);
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::handleSimplex1(std::vector<MinkowskiDifferenceResult> &simplex) {
    glm::vec3 B = simplex[0].m;
    glm::vec3 A = simplex[1].m;
    if (glm::dot(B-A, -A) > 0) {
        glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
        return std::pair<bool, glm::vec3>(false, new_dir);
    } else {
        simplex.clear();
        simplex.push_back(simplex[1]);
        return std::pair<bool, glm::vec3>(false, -A);
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::handleSimplex2(std::vector<MinkowskiDifferenceResult> &simplex) {
    glm::vec3 C = simplex[0].m;
    glm::vec3 B = simplex[1].m;
    glm::vec3 A = simplex[2].m;

    if(glm::dot(glm::cross(calculateTriangleNormal(A,B,C), C-A), -A) > 0) {
        if (glm::dot(C-A, -A) > 0) {
            simplex.clear();
            simplex.push_back(simplex[2]);
            simplex.push_back(simplex[0]);
            glm::vec3 new_dir = glm::cross(glm::cross(C-A, -A), C-A);
            return std::pair<bool, glm::vec3>(false, new_dir);
        } else {
            if (glm::dot(B-A, -A) > 0) {
                simplex.clear();
                simplex.push_back(simplex[2]);
                simplex.push_back(simplex[1]);
                glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
                return std::pair<bool, glm::vec3>(false, new_dir);
            } else {
                simplex.clear();
                simplex.push_back(simplex[2]);
                return std::pair<bool, glm::vec3>(false, -A);
            }
        }
    } else {
        if (glm::dot(-glm::cross(calculateTriangleNormal(A,B,C), B-A), -A) > 0) {
            if (glm::dot(B-A, -A) > 0) {
                simplex.clear();
                simplex.push_back(simplex[2]);
                simplex.push_back(simplex[1]);
                glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
                return std::pair<bool, glm::vec3>(false, new_dir);
            } else {
                simplex.clear();
                simplex.push_back(simplex[2]);
                return std::pair<bool, glm::vec3>(false, -A);
            }
        } else {
            if (glm::dot(calculateTriangleNormal(A,B,C),-A) > 0) {
                return std::pair<bool, glm::vec3>(false, calculateTriangleNormal(A,B,C));
            } else {
                simplex.clear();
                simplex.push_back(simplex[1]);
                simplex.push_back(simplex[0]);
                simplex.push_back(simplex[2]);
                return std::pair<bool, glm::vec3>(false, -calculateTriangleNormal(A,B,C));
            }
        }
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::handleSimplex3(std::vector<MinkowskiDifferenceResult> &simplex) {
    glm::vec3 D = simplex[0].m;
    glm::vec3 C = simplex[1].m;
    glm::vec3 B = simplex[2].m;
    glm::vec3 A = simplex[3].m;

    glm::vec3 ABC_normal = calculateTriangleNormal(A,B,C);
    glm::vec3 CBA_normal = calculateTriangleNormal(C,B,A);
    bool use_CBA = false;
    if (glm::dot(ABC_normal, D-A) > 0) { // take the normal facing away from the other vertex
        use_CBA = true;
    }

    glm::vec3 ACD_normal = calculateTriangleNormal(A,C,D);
    glm::vec3 DCA_normal = calculateTriangleNormal(D,C,A);
    bool use_DCA = false;
    if (glm::dot(ACD_normal, B-A) > 0) {
        use_DCA = true;
    }

    glm::vec3 ABD_normal = calculateTriangleNormal(A,B,D);
    glm::vec3 DBA_normal = calculateTriangleNormal(D,B,A);
    bool use_DBA = false;
    if (glm::dot(ABD_normal, C-A) > 0) {
        use_DBA = true;
    }

    if ((use_CBA && glm::dot(CBA_normal,-A) > 0) || (!use_CBA && glm::dot(ABC_normal,-A) > 0)) {
        simplex.clear();
        simplex.push_back(simplex[1]);
        simplex.push_back(simplex[2]);
        simplex.push_back(simplex[3]);
        return handleSimplex2(simplex);
    } else if ((use_DCA && glm::dot(DCA_normal,-A) > 0) || (!use_DCA && glm::dot(ACD_normal,-A) > 0)) {
        simplex.clear();
        simplex.push_back(simplex[1]);
        simplex.push_back(simplex[0]);
        simplex.push_back(simplex[3]);
        return handleSimplex2(simplex);
    } else {
        simplex.clear();
        simplex.push_back(simplex[0]);
        simplex.push_back(simplex[2]);
        simplex.push_back(simplex[3]);
        return handleSimplex2(simplex);
    }
}

bool PhysicsSystem::checkIfLineContainsPoint(const glm::vec3 l1, const glm::vec3 l2, const glm::vec3 p) {
    return (glm::length(l1-p) + glm::length(l2-p) - glm::length(l1-l2)) < EPSILON;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

bool sameTriangleSide(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 n)
{
    glm::vec3 normal = glm::cross(n, v3-v2);
    float dotV4 = glm::dot(normal, v3-v1);
    float dotP = glm::dot(normal, -v1);
    return sgn(dotV4) == sgn(dotP);
}

bool PhysicsSystem::checkIfTriangleContainsOrigin(const std::vector<glm::vec3> &simplex) {
    glm::vec3 A = simplex[0];
    glm::vec3 B = simplex[1];
    glm::vec3 C = simplex[2];
    glm::vec3 n = calculateTriangleNormal(A,B,C);
    return sameTriangleSide(A,B,C,n) &&
           sameTriangleSide(B,C,A,n) &&
           sameTriangleSide(C,A,B,n);
}

bool SameSide(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 v4)
{
    glm::vec3 normal = glm::cross(v2 - v1, v3 - v1);
    float dotV4 = glm::dot(normal, v4 - v1);
    float dotP = glm::dot(normal, -v1);
    return sgn(dotV4) == sgn(dotP);
}

bool PhysicsSystem::checkifTetrahedronContainsOrigin(const std::vector<MinkowskiDifferenceResult> &simplex) {
    glm::vec3 A = simplex[0].m;
    glm::vec3 B = simplex[1].m;
    glm::vec3 C = simplex[2].m;
    glm::vec3 D = simplex[3].m;
    return SameSide(A, B, C, D) &&
           SameSide(B, C, D, A) &&
           SameSide(C, D, A, B) &&
           SameSide(D, A, B, C);
}

glm::vec3 PhysicsSystem::calculateTriangleNormal(glm::vec3 a, glm::vec3 b, glm::vec3 c) {
    return glm::cross(b - a, c - a);
}

std::tuple<glm::vec3, glm::vec3, glm::vec3> PhysicsSystem::runExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1,
                                                                                std::shared_ptr<const RigidBodyComponent> rb2, const std::vector<MinkowskiDifferenceResult> &simplex) {
    std::vector<int> faces = {0,1,2,
                              1,2,3,
                              0,2,3,
                              0,1,3};
    glm::vec3 last_v;
    glm::vec3 last_rb1;
    glm::vec3 last_rb2;
    float difference = std::numeric_limits<float>::max();
    bool first_pass = true;
    std::vector<MinkowskiDifferenceResult> polytope(simplex.begin(), simplex.end());
    while (first_pass || difference > .01) {
        std::cout << difference << std::endl;

        // find face plane closest to origin
        float smallest_distance = std::numeric_limits<float>::max();
        int face_index = 0;
        for (int i = 0; i < faces.size(); i+=3) {
            glm::vec3 A = polytope[faces[i]].m;
            glm::vec3 B = polytope[faces[i+1]].m;
            glm::vec3 C = polytope[faces[i+2]].m;
            glm::vec3 normal = glm::cross(B-A, C-A);

            glm::mat3 mat = glm::transpose(glm::mat3(1,1,1,
                          glm::dot(A, B-A), glm::dot(B, B-A), glm::dot(C, B-A),
                          glm::dot(A, C-A), glm::dot(B, C-A), glm::dot(C, C-A)));
            glm::vec3 lambda = glm::inverse(mat) * glm::vec3(1,0,0);
            float sum = lambda.x + lambda.y + lambda.z;
            glm::vec3 thing = mat * lambda;

            float dist = std::abs(glm::dot(normal, A));
            if (dist < smallest_distance) {
                if (lambda.x >= 0 && lambda.y >= 0 && lambda.z >= 0) { // need to check that projection of origin on plane appears inside triangle, should find faster way to do this
                    smallest_distance = dist;
                    face_index = i;
                }
            }
        }

        glm::vec3 A = polytope[faces[face_index]].m;
        glm::vec3 B = polytope[faces[face_index + 1]].m;
        glm::vec3 C = polytope[faces[face_index + 2]].m;

        int f1 = faces[face_index]; // A
        int f2 = faces[face_index+1]; // B
        int f3 = faces[face_index+2]; // C

        glm::mat3 mat = glm::transpose(glm::mat3(1,1,1,
                        glm::dot(A, B-A), glm::dot(B, B-A), glm::dot(C, B-A),
                        glm::dot(A, C-A), glm::dot(B, C-A), glm::dot(C, C-A)));
        glm::vec3 lambda = glm::inverse(mat) * glm::vec3(1,0,0);

        assert(lambda.x >= 0 && lambda.y >= 0 && lambda.z >= 0);

        glm::vec3 v = lambda.x*A + lambda.y*B + lambda.z*C; //glm::vec3 v = (A + B + C) / 3.0f; //


        // calculate AB


        glm::vec3 AB = B-A;
        glm::vec3 AC = C-A;
        glm::vec3 BC = C-B;

        int p = polytope.size();

        glm::vec3 projAB = glm::dot(A, AB) / glm::dot(AB,AB) * AB;
        glm::vec3 AB_v = A - projAB;

        glm::vec3 projAC = glm::dot(A, AC) / glm::dot(AC,AC) * AC;
        glm::vec3 AC_v = A - projAC;

        glm::vec3 projBC = glm::dot(B, BC) / glm::dot(BC,BC) * BC;
        glm::vec3 BC_v = B - projAB;

        polytope.push_back(minkowskiDifferenceSupport(AB_v, rb1, rb2)); // p
        polytope.push_back(minkowskiDifferenceSupport(AC_v, rb1, rb2)); // p+1
        polytope.push_back(minkowskiDifferenceSupport(BC_v, rb1, rb2)); // p+2
        glm::vec3 w;
        glm::vec3 rb1_p;
        glm::vec3 rb2_p;
        MinkowskiDifferenceResult result = minkowskiDifferenceSupport(v, rb1, rb2);
        polytope.push_back(result); // p+3

        faces[face_index] = f1; faces[face_index+1] = p; faces[face_index+2] = p+3; // A, AB, w
        faces.push_back(f1); faces.push_back(p+1); faces.push_back(p+3); // A, AC, w
        faces.push_back(f3); faces.push_back(p+1); faces.push_back(p+3); // C, AC, w
        faces.push_back(f3); faces.push_back(p+2); faces.push_back(p+3); // C, BC, w
        faces.push_back(f2); faces.push_back(p+2); faces.push_back(p+3); // B, BC, w
        faces.push_back(f2); faces.push_back(p); faces.push_back(p+3); // B, AB, w

        if (!first_pass) {
            difference = glm::length(last_v - v);
        }
        //std::cout << glm::to_string(v) << std::endl;

        last_v = v;
        MinkowskiDifferenceResult A_w = polytope[faces[face_index]];
        MinkowskiDifferenceResult B_w = polytope[faces[face_index + 1]];
        MinkowskiDifferenceResult C_w = polytope[faces[face_index + 2]];
        last_rb1 = lambda.x*A_w.rb1 + lambda.y*B_w.rb1 + lambda.z*C_w.rb1; //last_rb1 = rb1_p; //last_rb1 = lambda.x*rb1_A + lambda.y*rb1_B + lambda.z*rb1_C;
        last_rb2 = lambda.x*A_w.rb2 + lambda.y*B_w.rb2 + lambda.z*C_w.rb2; // //last_rb2 = rb2_p; //last_rb2 = lambda.x*rb2_A + lambda.y*rb2_B + lambda.z*rb2_C;
        std::cout << "mtv" << glm::to_string(last_v) << std::endl;
        first_pass = false;
    }

    return std::make_tuple(last_v,last_rb1, last_rb2);
}


std::tuple<bool, glm::vec3, glm::vec3, glm::vec3> PhysicsSystem::runGJKAndExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2) {
    std::pair<bool, std::vector<MinkowskiDifferenceResult>> pair = runGJK(rb1, rb2);
    glm::vec3 ret_vec = glm::vec3(0);
    glm::vec3 rb1_p;
    glm::vec3 rb2_p;
    if (pair.first) {
        std::tie(ret_vec, rb1_p, rb2_p) = runExpandingPolytope(rb1, rb2, pair.second);
    }

    return std::make_tuple(pair.first, ret_vec, rb1_p, rb2_p);
}

