#include "CylinderRigidBodyComponent.h"

CylinderRigidBodyComponent::CylinderRigidBodyComponent(float diameter, float height, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color) :
   RigidBodyComponent(mass, collides, calculateIBody(diameter, height, mass), movable, use_gravity, color), m_height(height), m_diameter(diameter)
{
    m_graphics = Graphics::getGlobalInstance();
}

glm::mat3 CylinderRigidBodyComponent::calculateIBody(const float diameter, const float height, const float mass) {
    float r = diameter / 2.0;
    float Ix = mass * (3.0/20.0*r*r + 3.0/80.0*height*height);
    float Iz = Ix;
    float Iy = 3.0/10.0*mass*r*r;
    glm::mat3 mat = glm::mat3(Ix, 0.f, 0.f,
                              0.f, Iy, 0.f,
                              0.f,0.f, Iz);
    return mat;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

glm::vec3 CylinderRigidBodyComponent::support(glm::vec3 dir) const {
    float r = m_diameter / 2.0;
    float sigma = std::sqrt(dir.x*dir.x + dir.z*dir.z);
    if (sigma > 0) {
        return glm::vec3(r / sigma * dir.x, sgn(dir.y) * m_height / 2.0, r / sigma * dir.z);
    } else {
        return glm::vec3(0, sgn(dir.y) * m_height/2.0,0);
    }
}

/*
glm::vec3 CylinderRigidBodyComponent::support(glm::vec3 dir) const {
    float r = m_diameter / 2.0;
    float sin_alpha = r / std::sqrt(r*r + m_height*m_height);
    float sigma = std::sqrt(dir.x*dir.x + dir.z*dir.z);
    if (dir.y  > glm::length(dir) * sin_alpha) {
        return glm::vec3(0,m_height/2.0,0);
    } else if (sigma > 0) {
        return glm::vec3(r / sigma * dir.x, -m_height / 2.0, r / sigma * dir.z);
    } else {
        return glm::vec3(0,-m_height/2.0,0);
    }
}
*/

float CylinderRigidBodyComponent::getMaxRadius() const {
    return glm::length(glm::vec2(m_diameter / 2.0, .5*m_height));
}

void CylinderRigidBodyComponent::draw() {
    m_graphics->setShader("default");
    m_graphics->clearTransform();
    m_graphics->setColor(m_color);
    m_transform_component->setTransform();
    m_graphics->translate(glm::vec3(0,-.5,0));
    m_graphics->drawShape("cylinder");
}
