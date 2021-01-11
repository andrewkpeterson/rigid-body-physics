#ifndef CYLINDERRIGIDBODYCOMPONENT_H
#define CYLINDERRIGIDBODYCOMPONENT_H

#include "RigidBodyComponent.h"

class CylinderRigidBodyComponent : public RigidBodyComponent
{
public:
    CylinderRigidBodyComponent(float diameter, float height, float mass, bool collides, bool movable, bool use_gravity, glm::vec3 color);
    glm::mat3 calculateIBody(const float diameter, const float height, const float mass);
    void draw() override;
    glm::vec3 support(glm::vec3 dir) const override;
    float getMaxRadius() const override;

private:
    Graphics *m_graphics;
    float m_diameter;
    float m_height;
};

#endif // CYLINDERRIGIDBODYCOMPONENT_H
