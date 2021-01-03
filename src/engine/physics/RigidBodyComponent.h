#ifndef RIGIDBODYCOMPONENT_H
#define RIGIDBODYCOMPONENT_H

#include "engine/util/CommonIncludes.h"
#include "src/engine/physics/TransformComponent.h"

class RigidBodyComponent
{
public:
    RigidBodyComponent(float mass, glm::vec3 euler_angles, bool collides, bool movable);
    virtual void draw() = 0;  // NOTE: get rid of this when you port the code over to your engine
private:
    // constant quantities
    float m_mass;
    glm::mat3 m_I_body;
    glm::mat3 m_I_body_inverse;  // used to calculate the inertia tensor

    // state variables
    // vec3 for position (part of transform component)
    // vec3 for euler anglesglm::mat3 m_orientation;
    glm::vec3 m_linear_momentum;
    glm::vec3 m_angular_momentum;

    // derived quantities
    glm::mat3 m_intertia_tensor_inverse;
    glm::vec3 m_linear_velocity;
    glm::vec3 m_angular_velocity;

    // computed quantities
    glm::vec3 m_force;
    glm::vec3 m_torque;
};

#endif // RIGIDBODYCOMPONENT_H
