#ifndef PHYSICSSYSTEM_H
#define PHYSICSSYSTEM_H

#include "src/engine/util/CommonIncludes.h"
#include "src/engine/physics/RigidBodyComponent.h"
#include "src/engine/physics/CuboidRigidBodyComponent.h"
#include "src/engine/physics/CylinderRigidBodyComponent.h"
#include "src/engine/physics/EllipsoidRigidBodyComponent.h"

enum class PhysicsDebuggerObject;
enum class PhysicsDebuggerMode;

struct MinkowskiDifferenceResult {
    glm::vec3 m;
    glm::vec3 rb1;
    glm::vec3 rb2;
};

class PhysicsSystem
{
public:
    PhysicsSystem();

    void tick(float seconds);
    void draw(PhysicsDebuggerMode mode, glm::vec3 pos1, glm::vec3 pos2, glm::vec3 rot1, glm::vec3 rot2, glm::vec3 scale1,
              glm::vec3 scale2, PhysicsDebuggerObject obj1, PhysicsDebuggerObject obj2); // you will take this out when you port over your code
    std::vector<std::shared_ptr<RigidBodyComponent>> &getRigidBodies() { return m_rigid_bodies; }
    void resolveCollisions();
    MinkowskiDifferenceResult minkowskiDifferenceSupport(glm::vec3 dir, const std::shared_ptr<const RigidBodyComponent> rb1, const std::shared_ptr<const RigidBodyComponent> rb2);
    std::tuple<bool, glm::vec3, glm::vec3, glm::vec3> runGJKAndExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2);

    void checkCollisionForDebugging(glm::vec3 pos1, glm::vec3 pos2, glm::vec3 rot1, glm::vec3 rot2, glm::vec3 scale1, glm::vec3 scale2,
                                    PhysicsDebuggerObject obj1, PhysicsDebuggerObject obj2);

private:
    Graphics *m_graphics;
    bool colliding;
    glm::vec3 collision_spot;
    bool collision_happened;
    std::pair<bool, std::vector<MinkowskiDifferenceResult> > runGJK(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2);
    std::tuple<glm::vec3, glm::vec3, glm::vec3> runExpandingPolytope(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2, const std::vector<MinkowskiDifferenceResult> &simplex);
    std::tuple<glm::vec3, glm::vec3, glm::vec3> runExpandingPolytope2D(std::shared_ptr<const RigidBodyComponent> rb1, std::shared_ptr<const RigidBodyComponent> rb2, const std::vector<MinkowskiDifferenceResult> &simplex);
    std::pair<bool, glm::vec3> doSimplex(std::vector<MinkowskiDifferenceResult> &simplex);
    std::pair<bool, glm::vec3> handleSimplex1(std::vector<MinkowskiDifferenceResult> &simplex);
    std::pair<bool, glm::vec3> handleSimplex2(std::vector<MinkowskiDifferenceResult> &simplex);
    std::pair<bool, glm::vec3> handleSimplex3(std::vector<MinkowskiDifferenceResult> &simplex);
    glm::vec3 calculateTriangleNormal(glm::vec3 a, glm::vec3 b, glm::vec3 c);
    bool checkIfLineContainsPoint(const glm::vec3 l1, const glm::vec3 l2, const glm::vec3 p);
    bool checkIfTriangleContainsPoint(const glm::vec3 point, const std::vector<glm::vec3> &simplex);
    bool checkIfPointIsInPlane(const glm::vec3 point, const std::vector<glm::vec3> &simplex);
    bool checkifTetrahedronContainsOrigin(const std::vector<MinkowskiDifferenceResult> &simplex);
    void applyImpulse(std::shared_ptr<RigidBodyComponent> rb1, std::shared_ptr<RigidBodyComponent> rb2, const glm::vec3 mtv, const glm::vec3 contact_point);


    std::vector<std::shared_ptr<RigidBodyComponent>> m_rigid_bodies;

    const float EPSILON = .0001;
    const float BOUNCINESS = .05; // between 0 and 1
};

#endif // PHYSICSSYSTEM_H
