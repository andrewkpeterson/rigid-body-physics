#include "RigidBodyComponent.h"
#include "glm/gtx/string_cast.hpp"

RigidBodyComponent::RigidBodyComponent(float mass, bool collides, glm::mat3 I_body, bool movable, bool use_gravity, glm::vec3 color) :
    m_mass(mass), m_collides(collides), m_movable(movable), m_I_body(I_body), m_use_gravity(use_gravity), m_color(color)
{
    // TODO: when you port this code over, you should be adding this component in your gameobject constructor
    m_transform_component = std::make_unique<TransformComponent>(glm::vec3(0,0,0), glm::vec3(0,0,0), glm::vec3(1,1,1));

    // initialize member variables
    m_I_body_inverse = glm::inverse(m_I_body);
    glm::mat3 orientation = getOrientationMatrix();
    m_intertia_tensor_inverse = orientation * m_I_body * glm::transpose(orientation);
    m_linear_momentum = glm::vec3(0);
    m_angular_momentum = glm::vec3(0);
    m_linear_velocity = glm::vec3(0);
    m_angular_velocity = glm::inverse(m_intertia_tensor_inverse) * m_angular_momentum;
    m_orientation_derivative = star(m_angular_velocity) * orientation;
    m_force = glm::vec3(0);
    m_torque = glm::vec3(0);
}

glm::mat3 RigidBodyComponent::getOrientationMatrix() const {
    return glm::mat3(glm::rotate(m_transform_component->getEulerAngles().x, glm::vec3(1,0,0))
           * glm::rotate(m_transform_component->getEulerAngles().z, glm::vec3(0,0,1))
           * glm::rotate(m_transform_component->getEulerAngles().y, glm::vec3(0,1,0)));
}

glm::mat3 RigidBodyComponent::star(glm::vec3 &v) {
    return glm::transpose(glm::mat3(0, -v[2], v[1],
                                    v[2], 0, -v[0],
                                    -v[1], v[0], 0));
}

void RigidBodyComponent::calculateDerivatives() {

    // calculate linear velocity
    m_linear_velocity = m_linear_momentum / m_mass;

    //std::cout << "linear velocity: " << glm::to_string(m_linear_velocity) << std::endl;

    // calculate orientation derivative
    glm::mat3 orientation = getOrientationMatrix();
    m_intertia_tensor_inverse = orientation * m_I_body_inverse * glm::transpose(orientation);
    m_angular_velocity = m_intertia_tensor_inverse * m_angular_momentum;
    m_orientation_derivative = star(m_angular_velocity) * orientation;

    //std::cout << "euler " << glm::to_string(m_transform_component->getEulerAngles()) << std::endl;
    //std::cout << "I_body" << glm::to_string(m_I_body) << std::endl;
    //std::cout << "orientation " << glm::to_string(orientation) << std::endl;
    //std::cout << "inverse inertia tensor " << glm::to_string(m_intertia_tensor_inverse) << std::endl;
    //std::cout << "m_angular_velocity " << glm::to_string(m_angular_velocity) << std::endl;
    //std::cout << "m_orientation_derivative " << glm::to_string(m_orientation_derivative) << std::endl;


    // calculate m_force
    if (m_use_gravity) {
        applyGravity();
    }

    // calculate m_torque

}

glm::mat3 RigidBodyComponent::gramSchmidt(const glm::mat3 m) {
    glm::vec3 v0 = m[0];
    glm::vec3 v1 = m[1];
    glm::vec3 v2 = m[2];

    glm::vec3 w0 = glm::normalize(v0);
    glm::vec3 w1 = glm::normalize(v1 - glm::dot(v0, v1) / glm::dot(v0, v0) * v0);
    glm::vec3 w2 = v2 - (glm::dot(v0, v2) / glm::dot(v0, v0) * v0) - (glm::dot(v1, v2) / glm::dot(v1, v1) * v1);

    glm::mat3 ret;
    ret[0] = w0;
    ret[1] = w1;
    ret[2] = w2;

    return ret;
}


void RigidBodyComponent::applyDerivatives(float seconds) {
    addToPosition(seconds * m_linear_velocity);
    setOrientation(gramSchmidt(getOrientationMatrix() + m_orientation_derivative * seconds));
    addToLinearMomentum(seconds * m_force);
    addToAngularMomentum(seconds * m_torque);
}

void RigidBodyComponent::zeroOutDerivatives() {
    m_linear_velocity = glm::vec3(0);
    m_orientation_derivative = glm::mat3(0);
    m_force = glm::vec3(0);
    m_torque = glm::vec3(0);
}

void RigidBodyComponent::applyGravity() {
    m_force += glm::vec3(0,GRAVITY_ACCELERATION,0) * m_mass;
}
