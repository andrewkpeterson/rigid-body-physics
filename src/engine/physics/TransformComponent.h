#ifndef TRANSFORM_COMPONENT_H
#define TRANSFORM_COMPONENT_H

#include "engine/util/CommonIncludes.h"
#include "src/engine/graphics/Graphics.h"

class TransformComponent
{

public:
    TransformComponent(glm::vec3 position, glm::vec3 euler_angles, glm::vec3 scale);

    void setTransform();

    glm::vec3 getPosition() const { return m_position; }
    glm::vec3 getEulerAngles() const { return m_euler_angles; }
    glm::vec3 getScale() const { return m_scale; }

    void setPosition(glm::vec3 p) { m_position = p; }
    void setEulerAngles(glm::vec3 angles) { m_euler_angles = angles; }
    void setScale(glm::vec3 scale) { m_scale = scale; }

    glm::vec3 rotationMat2EulerAngles(const glm::mat3 &column_major_mat);

private:
    glm::vec3 m_position;
    // The gameobject is looking down the x axis in object space, and the up axis is y.
    // So, roll = m_rotation.x, pitch = m_rotation.z, and yaw = m_rotation.y.
    // We apply the intrinsic rotations in the order roll, pitch, and yaw.
    // Using matrices r_x, r_y and r_z, where r_a rotates about the global axis a,
    // we know that vertex_in_world_space = r_x * r_z * r_y * vertex_in_object_space.
    // in terms of gimbals, rx is the outer gimbal and ry is the inner gimbal.
    glm::vec3 m_euler_angles;
    glm::vec3 m_scale;

    Graphics *g;
};

#endif // TRANSFORM_COMPONENT_H
