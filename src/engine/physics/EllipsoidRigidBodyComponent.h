#ifndef ELLIPSOIDRIGIDBODYCOMPONENT_H
#define ELLIPSOIDRIGIDBODYCOMPONENT_H

#include "RigidBodyComponent.h"

class EllipsoidRigidBodyComponent : public RigidBodyComponent
{
public:
    EllipsoidRigidBodyComponent(glm::vec3 dimensions, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color);
    glm::mat3 calculateIBody(const glm::vec3 &d, const float mass);
    void draw() override;
    glm::vec3 support(glm::vec3 dir) const override;
    float getMaxRadius() const override;

private:
    Graphics *m_graphics;
    glm::vec3 m_dims;
};

#endif // ELLIPSOIDRIGIDBODYCOMPONENT_H
