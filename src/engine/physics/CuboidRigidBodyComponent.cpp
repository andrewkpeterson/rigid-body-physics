#include "CuboidRigidBodyComponent.h"

CuboidRigidBodyComponent::CuboidRigidBodyComponent(glm::vec3 dimensions, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color) :
    RigidBodyComponent(mass, collides, calculateIBody(dimensions, mass), movable, use_gravity, color, dimensions), m_dims(dimensions)
{
    m_graphics = Graphics::getGlobalInstance();
}

glm::mat3 CuboidRigidBodyComponent::calculateIBody(const glm::vec3 &d, const float mass) {
    glm::mat3 mat =  (mass / 12.f) * glm::mat3(d.y * d.y + d.z * d.z, 0.f, 0.f,
                                               0.f, d.x * d.x + d.z * d.z, 0.f,
                                               0.f,0.f, d.x * d.x + d.y * d.y);
    return mat;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

glm::vec3 CuboidRigidBodyComponent::support(glm::vec3 dir) const {
    return glm::vec3(float(sgn(dir.x)) * m_dims.x / 2.0f, float(sgn(dir.y)) * m_dims.y / 2.0f, float(sgn(dir.z)) * m_dims.z / 2.0f);
}

float CuboidRigidBodyComponent::getMaxRadius() const {
    float max_dim = std::max(std::max(m_dims[0] / 2.0, m_dims[1] / 2.0), m_dims[2] / 2.0);
    return max_dim * std::sqrt(2.0);
}

void CuboidRigidBodyComponent::draw() {
    m_graphics->setShader("default");
    m_graphics->clearTransform();
    m_graphics->setColor(m_color);
    m_transform_component->setTransform();
    m_graphics->drawShape("cube");
}
