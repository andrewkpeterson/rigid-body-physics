#include "EllipsoidRigidBodyComponent.h"

EllipsoidRigidBodyComponent::EllipsoidRigidBodyComponent(glm::vec3 dimensions, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color) :
    RigidBodyComponent(mass, collides, calculateIBody(dimensions, mass), movable, use_gravity, color), m_dims(dimensions)
{
    m_graphics = Graphics::getGlobalInstance();
}

glm::mat3 EllipsoidRigidBodyComponent::calculateIBody(const glm::vec3 &d, const float mass) {
    float Ix = 1.0/5.0 * mass * (d.y * d.y + d.z * d.z);
    float Iy = 1.0/5.0 * mass * (d.x * d.x + d.z * d.z);
    float Iz = 1.0/5.0 * mass * (d.x * d.x + d.y * d.y);
    glm::mat3 mat = glm::mat3(Ix, 0.f, 0.f,
                              0.f, Iy, 0.f,
                              0.f,0.f, Iz);
    return mat;
}

glm::vec3 EllipsoidRigidBodyComponent::support(glm::vec3 dir) const {
    dir = glm::normalize(dir);
    return glm::vec3(dir.x * m_dims.x / 2.0f, dir.y * m_dims.y / 2.0f, dir.z * m_dims.z / 2.0f);
}

float EllipsoidRigidBodyComponent::getMaxRadius() const {
    return std::max(std::max(m_dims[0] / 2.0, m_dims[1] / 2.0), m_dims[2] / 2.0);
}

void EllipsoidRigidBodyComponent::draw() {
    m_graphics->setShader("default");
    m_graphics->clearTransform();
    m_graphics->setColor(m_color);
    m_transform_component->setTransform();
    m_graphics->drawShape("sphere");
}

