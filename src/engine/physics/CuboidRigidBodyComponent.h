#ifndef CUBOIDRIGIDBODYCOMPONENT_H
#define CUBOIDRIGIDBODYCOMPONENT_H

#include "RigidBodyComponent.h"

class CuboidRigidBodyComponent : public RigidBodyComponent
{
public:
    CuboidRigidBodyComponent(glm::vec3 dimensions, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color);
    glm::mat3 calculateIBody(const glm::vec3 &dimensions, const float mass);
    void draw() override;
    glm::vec3 support(glm::vec3 dir) const override;
    float getMaxRadius() const override;

private:
    Graphics *m_graphics;
    glm::vec3 m_dims;
};

#endif // CUBOIDRIGIDBODYCOMPONENT_H
