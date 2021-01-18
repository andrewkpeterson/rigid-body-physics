#include "PhysicsSystem.h"
#include "src/engine/graphics/Camera.h"
#include <numeric>
#include <set>
#include <list>
#include <map>
#include <QFile>
#include <QTextStream>

PhysicsSystem::PhysicsSystem() : colliding(false), collision_spot(glm::vec3(0)), collision_happened(false)
{
    /*
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    //m_rigid_bodies.push_back(std::make_shared<EllipsoidRigidBodyComponent>(glm::vec3(1,3,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    //m_rigid_bodies.push_back(std::make_shared<CylinderRigidBodyComponent>(1, 1, 1, false, false, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[0]->setPosition(glm::vec3(0,2.0,0));
    m_rigid_bodies[0]->setEuelerAngles(glm::vec3(0,0,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[1]->setPosition(glm::vec3(0,3.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[2]->setPosition(glm::vec3(0,5.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[3]->setPosition(glm::vec3(0,7.5,0));


    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[4]->setPosition(glm::vec3(1.1,2.0,0));
    m_rigid_bodies[4]->setEuelerAngles(glm::vec3(0,0,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[5]->setPosition(glm::vec3(1.1,3.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[6]->setPosition(glm::vec3(1.1,5.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[7]->setPosition(glm::vec3(1.1,7.5,0));


    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[8]->setPosition(glm::vec3(-1.1,2.0,0));
    m_rigid_bodies[8]->setEuelerAngles(glm::vec3(0,0,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[9]->setPosition(glm::vec3(-1.1,3.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[10]->setPosition(glm::vec3(-1.1,5.5,0));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[11]->setPosition(glm::vec3(-1.1,7.5,0));

    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, false, false, glm::vec3(0,.7,.7)));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, false, false, glm::vec3(0,.7,.7)));
    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, false, false, glm::vec3(0,.7,.7)));

    m_rigid_bodies[12]->setPosition(glm::vec3(0,0,0));
    m_rigid_bodies[13]->setPosition(glm::vec3(1.1,0,0));
    m_rigid_bodies[14]->setPosition(glm::vec3(-1.1,0,0));

    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    m_rigid_bodies[15]->setPosition(glm::vec3(0,15,-100));
    m_rigid_bodies[15]->addToLinearMomentum(glm::vec3(0,0,20));
    m_rigid_bodies[15]->setEuelerAngles(glm::vec3(.5,0,.5));
    */

    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    //m_rigid_bodies.push_back(std::make_shared<EllipsoidRigidBodyComponent>(glm::vec3(1,3,1), 1, true, true, true, glm::vec3(.7,.7,0)));
    //m_rigid_bodies.push_back(std::make_shared<CylinderRigidBodyComponent>(1, 1, 1, false, false, true, glm::vec3(.7,.7,0)));

    m_rigid_bodies[0]->setPosition(glm::vec3(0,2.0,0));
    m_rigid_bodies[0]->setEuelerAngles(glm::vec3(0,0,.2));

    m_rigid_bodies.push_back(std::make_shared<CuboidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, false, false, glm::vec3(0,.7,.7)));


    //m_rigid_bodies.push_back(std::make_shared<EllipsoidRigidBodyComponent>(glm::vec3(1,1,1), 1, true, false, false, glm::vec3(.7,.7,0)));


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
                std::tie(had_collision, mtv, rb1_p, rb2_p) = runGJKAndExpandingPolytope(m_rigid_bodies[i], m_rigid_bodies[j]);
                if (had_collision) {
                    // apply linear and rotational impulse
                    applyImpulse(m_rigid_bodies[i], m_rigid_bodies[j], mtv, rb1_p);
                }
            } else {
                colliding = false;
            }
        }
    }

}

void PhysicsSystem::applyImpulse(std::shared_ptr<RigidBodyComponent> rb1, std::shared_ptr<RigidBodyComponent> rb2, const glm::vec3 mtv, const glm::vec3 contact_point) {

    rb1->calculateDerivatives();
    rb2->calculateDerivatives();

    // calculate the velocity of the contact points
    glm::vec3 rb1_contact_velocity = rb1->getLinearVelocity() + glm::cross(rb1->getAngularVelocity(), contact_point - rb1->getPosition());
    glm::vec3 rb2_contact_velocity = rb2->getLinearVelocity() + glm::cross(rb2->getAngularVelocity(), contact_point - rb2->getPosition());
    //std::cout << "rb1_contact_velocity: " << glm::to_string(rb1_contact_velocity) << std::endl;
    //std::cout << "rb2_contact_velocity: " << glm::to_string(rb2_contact_velocity) << std::endl;

    // the magnitude of the relative velocity in the direction of the normalized mtv
    glm::vec3 norm_mtv = glm::length(mtv) > 0 ? glm::normalize(mtv) : glm::vec3(0);
    if (glm::length(mtv) > 0) {
        float relative_velocity = glm::dot(norm_mtv, (rb1_contact_velocity - rb2_contact_velocity)); // might need a negative here!
        std::cout << "mtv: " << glm::to_string(mtv) << std::endl;
        //std::cout << "mtv.z: " << mtv.z << std::endl;
        //std::cout << "position: " << glm::to_string(rb1->getPosition()) << std::endl;
        //std::cout << "relative velocity: " << glm::to_string(relative_velocity) << std::endl;
        float numerator = -(1 + BOUNCINESS) * relative_velocity;

        float term1 = 1.0 / rb1->getMass();
        float term2 = 1.0 / rb2->getMass();
        glm::vec3 r_rb1 = contact_point - rb1->getPosition();
        glm::vec3 r_rb2 = contact_point - rb2->getPosition();
        float term3 = glm::dot(norm_mtv, glm::cross(rb1->getInverseInertiaTensor() * (glm::cross(r_rb1, norm_mtv)), r_rb1));
        float term4 = glm::dot(norm_mtv, glm::cross(rb2->getInverseInertiaTensor() * (glm::cross(r_rb2, norm_mtv)), r_rb2));

        float j = numerator / (term1 + term2 + term3 + term4);
        //std::cout << "j: " << glm::to_string(j) << std::endl;
        glm::vec3 impulse = j * norm_mtv;

        //std::cout << "impulse: " << glm::to_string(impulse) << std::endl;

        if (rb1->getMovable()) {
            rb1->addToPosition(-mtv);
            rb1->addToLinearMomentum(impulse);
            rb1->addToAngularMomentum(glm::cross(r_rb1, impulse));
        }

        if (rb2->getMovable()) {
            rb2->addToPosition(mtv);
            rb2->addToLinearMomentum(-impulse);
            rb2->addToAngularMomentum(glm::cross(r_rb2, -impulse));
        }
    }
    rb1->zeroOutDerivatives();
    rb2->zeroOutDerivatives();
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

bool PhysicsSystem::checkIfPointIsInPlane(const glm::vec3 point, const std::vector<glm::vec3> &simplex) {
    glm::vec3 normal = glm::cross(simplex[1] - simplex[0], simplex[2] - simplex[0]);
    return std::abs(glm::dot(point - simplex[0], normal)) < .001;
}

std::pair<bool, glm::vec3> PhysicsSystem::doSimplex(std::vector<MinkowskiDifferenceResult> &simplex) {
    assert(simplex.size() >= 1 && simplex.size() <= 4);

    if (simplex.size() == 1 && glm::length((*simplex.begin()).m) < EPSILON) {
        std::cout << "finished 0" << std::endl;
        return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
    else if (simplex.size() == 2 && checkIfLineContainsPoint(simplex[0].m, simplex[1].m, glm::vec3(0,0,0))) {
        std::cout << "finished 1" << std::endl;
        return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
    else if (simplex.size() == 3 && checkIfTriangleContainsPoint(glm::vec3(0,0,0), std::vector<glm::vec3>({simplex[0].m, simplex[1].m, simplex[2].m})) &&
             checkIfPointIsInPlane(glm::vec3(0,0,0), std::vector<glm::vec3>({simplex[0].m, simplex[1].m, simplex[2].m}))) {
        std::cout << "finished 2" << std::endl;
        return std::pair<bool, glm::vec3>(true, glm::vec3(0));
    }
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
    MinkowskiDifferenceResult B_r = simplex[0];
    MinkowskiDifferenceResult A_r = simplex[1];
    glm::vec3 B = simplex[0].m;
    glm::vec3 A = simplex[1].m;
    if (glm::dot(B-A, -A) > 0) {
        glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
        return std::pair<bool, glm::vec3>(false, new_dir);
    } else {
        simplex.clear();
        simplex.push_back(A_r);
        return std::pair<bool, glm::vec3>(false, -A);
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::handleSimplex2(std::vector<MinkowskiDifferenceResult> &simplex) {
    MinkowskiDifferenceResult C_r = simplex[0];
    MinkowskiDifferenceResult B_r = simplex[1];
    MinkowskiDifferenceResult A_r = simplex[2];

    glm::vec3 C = simplex[0].m;
    glm::vec3 B = simplex[1].m;
    glm::vec3 A = simplex[2].m;

    if(glm::dot(glm::cross(calculateTriangleNormal(A,B,C), C-A), -A) > 0) {
        if (glm::dot(C-A, -A) > 0) {
            simplex.clear();
            simplex.push_back(A_r);
            simplex.push_back(C_r);
            glm::vec3 new_dir = glm::cross(glm::cross(C-A, -A), C-A);
            return std::pair<bool, glm::vec3>(false, new_dir);
        } else {
            if (glm::dot(B-A, -A) > 0) {
                simplex.clear();
                simplex.push_back(A_r);
                simplex.push_back(B_r);
                glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
                return std::pair<bool, glm::vec3>(false, new_dir);
            } else {
                simplex.clear();
                simplex.push_back(A_r);
                return std::pair<bool, glm::vec3>(false, -A);
            }
        }
    } else {
        if (glm::dot(-glm::cross(calculateTriangleNormal(A,B,C), B-A), -A) > 0) {
            if (glm::dot(B-A, -A) > 0) {
                simplex.clear();
                simplex.push_back(A_r);
                simplex.push_back(B_r);
                glm::vec3 new_dir = glm::cross(glm::cross(B-A, -A), B-A);
                return std::pair<bool, glm::vec3>(false, new_dir);
            } else {
                simplex.clear();
                simplex.push_back(A_r);
                return std::pair<bool, glm::vec3>(false, -A);
            }
        } else {
            if (glm::dot(calculateTriangleNormal(A,B,C),-A) > 0) {
                return std::pair<bool, glm::vec3>(false, calculateTriangleNormal(A,B,C));
            } else {
                simplex.clear();
                simplex.push_back(B_r);
                simplex.push_back(C_r);
                simplex.push_back(A_r);
                return std::pair<bool, glm::vec3>(false, -calculateTriangleNormal(A,B,C));
            }
        }
    }
}

std::pair<bool, glm::vec3> PhysicsSystem::handleSimplex3(std::vector<MinkowskiDifferenceResult> &simplex) {
    MinkowskiDifferenceResult D_r = simplex[0];
    MinkowskiDifferenceResult C_r = simplex[1];
    MinkowskiDifferenceResult B_r = simplex[2];
    MinkowskiDifferenceResult A_r = simplex[3];

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
        simplex.push_back(C_r);
        simplex.push_back(B_r);
        simplex.push_back(A_r);
        return handleSimplex2(simplex);
    } else if ((use_DCA && glm::dot(DCA_normal,-A) > 0) || (!use_DCA && glm::dot(ACD_normal,-A) > 0)) {
        simplex.clear();
        simplex.push_back(C_r);
        simplex.push_back(D_r);
        simplex.push_back(A_r);
        return handleSimplex2(simplex);
    } else {
        simplex.clear();
        simplex.push_back(D_r);
        simplex.push_back(B_r);
        simplex.push_back(A_r);
        return handleSimplex2(simplex);
    }
}

bool PhysicsSystem::checkIfLineContainsPoint(const glm::vec3 l1, const glm::vec3 l2, const glm::vec3 p) {
    return (glm::length(l1-p) + glm::length(l2-p) - glm::length(l1-l2)) < EPSILON;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

bool sameTriangleSide(glm::vec3 point, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 n)
{
    glm::vec3 normal = glm::cross(n, v3-v2);
    float dotV4 = glm::dot(normal, v1-v3);
    float dotP = glm::dot(normal, point-v3);
    return sgn(dotV4) == sgn(dotP);
}

bool PhysicsSystem::checkIfTriangleContainsPoint(const glm::vec3 point, const std::vector<glm::vec3> &simplex) {
    glm::vec3 A = simplex[0];
    glm::vec3 B = simplex[1];
    glm::vec3 C = simplex[2];
    glm::vec3 n = calculateTriangleNormal(A,B,C);
    return sameTriangleSide(point,A,B,C,n) &&
           sameTriangleSide(point,B,C,A,n) &&
           sameTriangleSide(point,C,A,B,n);
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

glm::vec3 Barycentric(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
    glm::vec3 v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    return glm::vec3(1.0f - v - w, v, w);

}

std::tuple<glm::vec3, glm::vec3, glm::vec3> PhysicsSystem::runExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1,
                                                                                std::shared_ptr<const RigidBodyComponent> rb2, const std::vector<MinkowskiDifferenceResult> &simplex) {
    std::vector<int> faces = {0,1,2,
                              1,2,3,
                              0,2,3,
                              0,1,3};
    std::vector<int> missing = {3,0,1,2};
    glm::vec3 closest_v(std::numeric_limits<float>::max());
    glm::vec3 best_rb1;
    glm::vec3 best_rb2;
    float best_distance = std::numeric_limits<float>::max();
    bool first_pass = true;
    std::vector<MinkowskiDifferenceResult> polytope(simplex.begin(), simplex.end());
    int iteration = 0;


    // record the current polytope
    QString name = "../polytope";
    name.append("original");
    name.append(".txt");
    QFile polytope_file(name);
    if (polytope_file.open(QIODevice::ReadWrite | QFile::Truncate)) {
        QTextStream stream(&polytope_file);
        char buf[512];
        stream << "vertices=[" << endl;
        for (int polytope_i = 0; polytope_i < polytope.size(); polytope_i++) {
            std::sprintf(buf, "[%f, %f, %f],", polytope[polytope_i].m.x, polytope[polytope_i].m.y, polytope[polytope_i].m.z);
            stream << buf << endl;
        }
        stream << "]," << endl;

        stream << "faces=[" << endl;
        for (int face_i = 0; face_i < faces.size(); face_i+=3) {
            std::sprintf(buf, "[%d, %d, %d],", faces[face_i], faces[face_i+1], faces[face_i+2]);
            stream << buf << endl;
        }
        stream << "]" << endl;
    }



    while (true) {

        // find face plane closest to origin
        float smallest_distance = std::numeric_limits<float>::max();
        int face_index = 0;
        float collinear_check;
        for (int i = 0; i < faces.size(); i+=3) {
            glm::vec3 A = polytope[faces[i]].m;
            glm::vec3 B = polytope[faces[i+1]].m;
            glm::vec3 C = polytope[faces[i+2]].m;

            glm::vec3 normal = glm::normalize(glm::cross(B-A, C-A));

            glm::vec3 v = glm::dot(A, normal) / glm::dot(normal, normal) * normal;
            double dist = std::abs(glm::dot(normal, A));

            //assert(dist < 1e30);

            //if (checkIfTriangleContainsOrigin(std::vector<glm::vec3>({A,B,C}))) { //if (lambda.x >= 0 && lambda.y >= 0 && lambda.z >= 0) { // need to check that projection of origin on plane appears inside triangle, should find faster way to do this
            if (dist < smallest_distance) {
                smallest_distance = dist;
                face_index = i;
                collinear_check = glm::length(glm::cross(B-A,C-A));
            }
            //}
        }

        MinkowskiDifferenceResult A_w = polytope[faces[face_index]];
        MinkowskiDifferenceResult B_w = polytope[faces[face_index + 1]];
        MinkowskiDifferenceResult C_w = polytope[faces[face_index + 2]];


        glm::vec3 normal = glm::normalize(glm::cross(B_w.m-A_w.m, C_w.m-A_w.m));

        // search for vertex not in the current face we are looking at
        std::vector<int> current_faces = {faces[face_index], faces[face_index+1], faces[face_index+2]};
        glm::vec3 vertex_not_in_face = glm::vec3(0,0,0);
        for (int face_i = 0; face_i < faces.size(); face_i++) {
            bool in_face = false;
            for (int j = 0; j < current_faces.size(); j++) {
                bool lies_in_face = std::abs(glm::dot(polytope[faces[face_i]].m - A_w.m, normal)) < .001;
                if (faces[face_i] == current_faces[j] || lies_in_face) {
                    in_face = true;
                }
            }
            if (!in_face) {
                vertex_not_in_face = polytope[faces[face_i]].m;
                break;
            }
        }
        //assert(glm::length(vertex_not_in_face) > .0000001);

        if (glm::dot(A_w.m, normal) < 0) {
            normal = -normal;
        }

        glm::vec3 v = glm::dot(A_w.m, normal) / glm::dot(normal, normal) * normal;

        glm::vec3 rb1_p;
        glm::vec3 rb2_p;
        MinkowskiDifferenceResult result = minkowskiDifferenceSupport(v, rb1, rb2);

        MinkowskiDifferenceResult other_result = minkowskiDifferenceSupport(-v, rb1, rb2);

        glm::vec3 new_w = result.m;

        // if v and the projection of w onto v are close to the same length, then we know that the triangle contains the boundary point of the minkowski difference
        // corresponding to the mtv. We then know that v is close to the mtv, so we take the projection of w onto v, since that point is on the boundary
        glm::vec3 proj_w_onto_v = glm::dot(v,new_w) / glm::dot(v,v) * v;

        //std::cout << std::abs(glm::length(v - proj_w_onto_v)) << std::endl;

        float dot_v_and_proj_w_onto_v = glm::dot(proj_w_onto_v, v);
        //assert(dot_v_and_proj_w_onto_v > 0);

        float dot_plane_and_v = glm::dot(A_w.m-B_w.m, v);
        //assert(std::abs(dot_plane_and_v) < .000001);

        if (!first_pass && std::abs(glm::length(v - proj_w_onto_v)) < best_distance) {
            closest_v = v; //closest_v = proj_w_onto_v;
            glm::vec3 lambda = Barycentric(v, A_w.m,B_w.m,C_w.m);
            best_rb1 = lambda.x*A_w.rb1 + lambda.y*B_w.rb1 + lambda.z*C_w.rb1; //last_rb1 = rb1_p; //last_rb1 = lambda.x*rb1_A + lambda.y*rb1_B + lambda.z*rb1_C;
            best_rb2 = lambda.x*A_w.rb2 + lambda.y*B_w.rb2 + lambda.z*C_w.rb2;
            best_distance = glm::length(v - proj_w_onto_v);
        } else if (first_pass) {
            closest_v = v; //closest_v = proj_w_onto_v;
            glm::vec3 lambda = Barycentric(v, A_w.m,B_w.m,C_w.m);
            best_rb1 = lambda.x*A_w.rb1 + lambda.y*B_w.rb1 + lambda.z*C_w.rb1; //last_rb1 = rb1_p; //last_rb1 = lambda.x*rb1_A + lambda.y*rb1_B + lambda.z*rb1_C;
            best_rb2 = lambda.x*A_w.rb2 + lambda.y*B_w.rb2 + lambda.z*C_w.rb2;
            best_distance = glm::length(v - proj_w_onto_v);
        }

        if (std::abs(glm::length(v - proj_w_onto_v)) < .01) {
            closest_v = v;
            glm::vec3 lambda = Barycentric(v, A_w.m,B_w.m,C_w.m);
            best_rb1 = lambda.x*A_w.rb1 + lambda.y*B_w.rb1 + lambda.z*C_w.rb1; //last_rb1 = rb1_p; //last_rb1 = lambda.x*rb1_A + lambda.y*rb1_B + lambda.z*rb1_C;
            best_rb2 = lambda.x*A_w.rb2 + lambda.y*B_w.rb2 + lambda.z*C_w.rb2;
            return std::make_tuple(closest_v,best_rb1,best_rb2);
        }

        // make sure that we don't add a vertex already in the polytope
        // WE REALLY SHOULD BE ASSERTING THAT THIS IS NOT GOING TO HAPPEN
        bool w_already_in_polytope = false;
        for (int i = 0; i < polytope.size(); i++) {
            if (glm::length(polytope[i].m - new_w) < .0001) {
                for (int j = 0; j < faces.size(); j++) {
                    if (i == faces[j]) {
                        w_already_in_polytope = true; //return std::make_tuple(polytope[i].m, polytope[i].rb1, polytope[i].rb2);
                        //return std::make_tuple(closest_v,best_rb1,best_rb2);
                    }
                }
            }
        }

        // make sure that we are not adding a vertex that is already part of a face in the polytope
        bool w_already_on_face_of_polytope = false;
        if (!w_already_in_polytope) {
            for (int i = 0; i < faces.size(); i+=3) {
                glm::vec3 normal = glm::normalize(glm::cross(polytope[faces[i+1]].m - polytope[faces[i]].m, polytope[faces[i+2]].m - polytope[faces[i]].m));
                std::vector<glm::vec3> triangle = {polytope[faces[i]].m, polytope[faces[i+1]].m, polytope[faces[i+2]].m};
                if (glm::abs(glm::dot(polytope[faces[i]].m - new_w, normal)) < .0001 && checkIfTriangleContainsPoint(polytope[faces[i]].m,triangle)) {
                    w_already_on_face_of_polytope = true;
                    glm::vec3 face_A = polytope[faces[i]].m;
                    glm::vec3 face_B = polytope[faces[i + 1]].m;
                    glm::vec3 face_C = polytope[faces[i + 2]].m;
                    int A_i = faces[i];
                    int B_i = faces[i+1];
                    int C_i = faces[i+2];
                    int p = polytope.size();
                    polytope.push_back(result);
                    faces[i] = A_i; faces[i] = B_i; faces[i] = p;
                    faces.push_back(A_i); faces.push_back(B_i); faces.push_back(p);
                    faces.push_back(B_i); faces.push_back(C_i); faces.push_back(p);

                    // record the current polytope
                    QString name = "../polytope";
                    name.append(std::to_string(iteration).c_str());
                    name.append(".txt");
                    QFile polytope_file(name);
                    if (polytope_file.open(QIODevice::ReadWrite | QFile::Truncate)) {
                        QTextStream stream(&polytope_file);
                        char buf[512];
                        stream << "vertices=[" << endl;
                        for (int polytope_i = 0; polytope_i < polytope.size(); polytope_i++) {
                            std::sprintf(buf, "[%f, %f, %f],", polytope[polytope_i].m.x, polytope[polytope_i].m.y, polytope[polytope_i].m.z);
                            stream << buf << endl;
                        }
                        stream << "]," << endl;

                        stream << "faces=[" << endl;
                        for (int face_i = 0; face_i < faces.size(); face_i+=3) {
                            std::sprintf(buf, "[%d, %d, %d],", faces[face_i], faces[face_i+1], faces[face_i+2]);
                            stream << buf << endl;
                        }
                        stream << "]" << endl;
                    }
                    std::cout << "split face in middle" << std::endl;
                    break;
                    //return std::make_tuple(closest_v,best_rb1,best_rb2);
                }
            }
        }

        if (!w_already_in_polytope && !w_already_on_face_of_polytope) {
            int p = polytope.size();
            polytope.push_back(result);
            std::vector<int> new_faces;
            std::map<std::pair<int,int>,int> edge_counts;
            std::set<std::pair<int, int>> edges;
            for (int i = 0; i < faces.size(); i += 3) {
                glm::vec3 face_A = polytope[faces[i]].m;
                glm::vec3 face_B = polytope[faces[i + 1]].m;
                glm::vec3 face_C = polytope[faces[i + 2]].m;
                glm::vec3 normal = glm::normalize(glm::cross(face_B - face_A, face_C - face_A));

                // search for vertex not in the current face we are looking at
                std::vector<int> current_faces = {faces[i], faces[i+1], faces[i+2]};
                glm::vec3 vertex_not_in_face = glm::vec3(0,0,0);
                for (int face_i = 0; face_i < faces.size(); face_i++) {
                    bool in_face = false;
                    for (int j = 0; j < current_faces.size(); j++) {
                        bool lies_in_face = std::abs(glm::dot(polytope[faces[face_i]].m - face_A, normal)) < .001;
                        if (faces[face_i] == current_faces[j] || lies_in_face) {
                            in_face = true;
                        }
                    }
                    if (!in_face) {
                        vertex_not_in_face = polytope[faces[face_i]].m;
                        break;
                    }
                }
                //assert(glm::length(vertex_not_in_face) > .0000001);

                if (glm::dot(normal, vertex_not_in_face) > 0 && glm::dot(normal, face_A) > 0) {
                    //assert(false);
                }


                if (glm::dot(normal, face_A) < 0) {
                    //std::cout << "hello" << std::endl;
                    normal = -normal;
                }

                //bool w_on_edge = checkIfLineContainsPoint(face_A, face_B, result.m) || checkIfLineContainsPoint(face_A, face_C, result.m) ||
                        //checkIfLineContainsPoint(face_C, face_B, result.m);

                std::cout << glm::to_string(face_A) << " " << glm::to_string(face_B) << " " << glm::to_string(face_C) << std::endl;
                std::cout << glm::dot(glm::normalize(normal), glm::normalize(new_w)) << std::endl;
                std::cout << "normal " << glm::to_string(normal) << std::endl;
                if (glm::dot(glm::normalize(normal), glm::normalize(new_w)) > .54 && glm::dot(glm::normalize(normal), glm::normalize(new_w)) < .55) {
                    std::cout << "stop" << std::endl;
                }

                if (glm::dot(glm::normalize(normal), glm::normalize(new_w - face_A)) > 0) {
                    //std::cout << "found face" << std::endl;
                    //std::cout << glm::to_string(face_A) << " " << glm::to_string(face_B) << " " << glm::to_string(face_C) << std::endl;
                    //std::cout << glm::dot(glm::normalize(normal), glm::normalize(new_w)) << std::endl;
                    //std::cout << glm::to_string(result.m) << std::endl;
                    int A_i = faces[i];
                    int B_i = faces[i+1];
                    int C_i = faces[i+2];

                    int minAB = std::min(A_i, B_i);
                    int maxAB = std::max(A_i, B_i);
                    int minAC = std::min(A_i, C_i);
                    int maxAC = std::max(A_i, C_i);
                    int minBC = std::min(C_i, B_i);
                    int maxBC = std::max(C_i, B_i);


                    if (edges.find(std::pair<int,int>(minAB, maxAB)) != edges.end()) {
                        edges.erase(edges.find(std::pair<int,int>(minAB, maxAB)));
                    } else {
                        edges.insert(std::pair<int,int>(minAB, maxAB));
                    }

                    if (edges.find(std::pair<int,int>(minAC, maxAC)) != edges.end()) {
                        edges.erase(edges.find(std::pair<int,int>(minAC, maxAC)));
                    } else {
                        edges.insert(std::pair<int,int>(minAC, maxAC));
                    }

                    if (edges.find(std::pair<int,int>(minBC, maxBC)) != edges.end()) {
                        edges.erase(edges.find(std::pair<int,int>(minBC, maxBC)));
                    } else {
                        edges.insert(std::pair<int,int>(minBC, maxBC));
                    }

                } else {
                    new_faces.push_back(faces[i]); new_faces.push_back(faces[i+1]); new_faces.push_back(faces[i+2]);
                }
            }

            for (auto it = edges.begin(); it != edges.end(); it++) {
                new_faces.push_back(it->first); new_faces.push_back(it->second); new_faces.push_back(p);
            }

            std::cout << "created new faces and repaired" << std::endl;
            std::cout << "*******************************" << std::endl;
            faces = new_faces;

            // record the current polytope
            QString name = "../polytope";
            name.append(std::to_string(iteration).c_str());
            name.append(".txt");
            QFile polytope_file(name);
            if (polytope_file.open(QIODevice::ReadWrite | QFile::Truncate)) {
                QTextStream stream(&polytope_file);
                char buf[512];
                stream << "vertices=[" << endl;
                for (int polytope_i = 0; polytope_i < polytope.size(); polytope_i++) {
                    std::sprintf(buf, "[%f, %f, %f],", polytope[polytope_i].m.x, polytope[polytope_i].m.y, polytope[polytope_i].m.z);
                    stream << buf << endl;
                }
                stream << "]," << endl;

                stream << "faces=[" << endl;
                for (int face_i = 0; face_i < faces.size(); face_i+=3) {
                    std::sprintf(buf, "[%d, %d, %d],", faces[face_i], faces[face_i+1], faces[face_i+2]);
                    stream << buf << endl;
                }
                stream << "]" << endl;
            }
        }

        iteration++;
        first_pass = false;
    }

    return std::make_tuple(closest_v,best_rb1,best_rb2);

}


std::tuple<bool, glm::vec3, glm::vec3, glm::vec3> PhysicsSystem::runGJKAndExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2) {
    std::pair<bool, std::vector<MinkowskiDifferenceResult>> pair = runGJK(rb1, rb2);
    glm::vec3 ret_vec = glm::vec3(0);
    glm::vec3 rb1_p;
    glm::vec3 rb2_p;
    if (pair.first && pair.second.size() > 3) {
        std::cout << "running epa" << std::endl;
        std::tie(ret_vec, rb1_p, rb2_p) = runExpandingPolytope(rb1, rb2, pair.second);
    }

    return std::make_tuple(pair.first && pair.second.size() > 3, ret_vec, rb1_p, rb2_p);
}

