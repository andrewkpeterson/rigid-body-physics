#ifndef RIGIDBODYCOMPONENT_H
#define RIGIDBODYCOMPONENT_H

#include "engine/util/CommonIncludes.h"
#include "src/engine/physics/TransformComponent.h"
#include "PhysicsSolver.h"

class RigidBodyComponent
{
public:
    RigidBodyComponent(float mass, bool collides, glm::mat3 I_body, bool movable, bool use_gravity, glm::vec3 color);

    glm::mat3 getOrientationMatrix() const;

    virtual void draw() = 0;  // NOTE: TODO: get rid of this when you port the code over to your engine
    virtual glm::vec3 support(glm::vec3 dir) const = 0;
    virtual float getMaxRadius() const = 0;

    void calculateDerivatives();
    void applyDerivatives(float seconds);
    void zeroOutDerivatives();

    static glm::mat3 gramSchmidt(const glm::mat3 m);

    void addToPosition(glm::vec3 p) { m_transform_component->setPosition(m_transform_component->getPosition() + p); }
    void setPosition(glm::vec3 p) { m_transform_component->setPosition(p); }
    glm::vec3 getPosition() const { return m_transform_component->getPosition(); }

    void setOrientation(glm::mat3 m) { m_transform_component->setEulerAngles(m_transform_component->rotationMat2EulerAngles(m)); }
    void setEuelerAngles(glm::vec3 euler) { m_transform_component->setEulerAngles(euler); }

    void addToLinearMomentum(glm::vec3 p) { m_linear_momentum += p; }
    glm::vec3 getLinearMomentum() const { return m_linear_momentum; }

    void addToAngularMomentum(glm::vec3 p) { m_angular_momentum += p; }
    glm::vec3 getAngularMomentum() const { return m_angular_momentum; }

    glm::vec3 getLinearVelocity() const { return m_linear_velocity; }
    glm::mat3 getOrientationDerivative() const { return m_orientation_derivative; }
    glm::vec3 getForce() const { return m_force; }
    glm::vec3 getTorque() const { return m_torque; }

    bool getCollides() { return m_collides; }
    bool getMovable() { return m_movable; }

    // constants
    const float GRAVITY_ACCELERATION = -.5;

protected:


    void applyGravity();

    // helper methods
    glm::mat3 star(glm::vec3 &v);

    // constant quantities
    float m_mass;
    bool m_collides;
    bool m_movable;
    bool m_use_gravity;
    glm::mat3 m_I_body;
    glm::mat3 m_I_body_inverse;  // used to calculate the inertia tensor

    // state vector
    // vec3 for position (part of transform component)
    // glm::mat3 m_orientation (euler angles part of transform component)
    glm::vec3 m_linear_momentum;
    glm::vec3 m_angular_momentum;

    glm::mat3 m_intertia_tensor_inverse; // this is used to calculate angular velocity
    glm::vec3 m_angular_velocity; // this isused to calculate orientation derivative

    // derivative of state vector
    glm::vec3 m_linear_velocity;
    glm::mat3 m_orientation_derivative;
    glm::vec3 m_force;
    glm::vec3 m_torque;

    // transform component
    std::unique_ptr<TransformComponent> m_transform_component;

    glm::vec3 m_color;

};

#endif // RIGIDBODYCOMPONENT_H
