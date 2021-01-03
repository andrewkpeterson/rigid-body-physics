#ifndef TRANSFORM_COMPONENT_H
#define TRANSFORM_COMPONENT_H

#include "engine/util/CommonIncludes.h"
#include "src/engine/graphics/Graphics.h"

class TransformComponent
{

public:
    TransformComponent(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale);
    void setTransform();

private:
    glm::vec3 m_position;
    // The gameobject is looking down the x axis in object space, and the up axis is y.
    // So, roll = m_rotation.x, pitch = m_rotation.z, and yaw = m_rotation.y.
    // We apply the intrinsic rotations in the order roll, pitch, and yaw.
    // Using matrices r_x, r_y and r_z, where r_a rotates about the global axis a,
    // we know that vertex_in_world_space = r_y * r_z * r_x * vertex_in_object_space, where rx is the
    // rotation matrix that rotates about the global x axis
    // in terms of gimbals, ry is the outer gimbal and rx is the inner gimbal
    glm::vec3 m_rotation;
    glm::vec3 m_scale;

    Graphics *g;
};

#endif // TRANSFORM_COMPONENT_H
