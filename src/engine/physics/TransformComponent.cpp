#include "TransformComponent.h"

TransformComponent::TransformComponent(glm::vec3 position, glm::vec3 euler_angles, glm::vec3 scale) :
    m_position(position), m_euler_angles(euler_angles), m_scale(scale)
{
    g = Graphics::getGlobalInstance();
}

void TransformComponent::setTransform() {
    g->translate(m_position);
    g->rotate(m_euler_angles.x, glm::vec3(1,0,0));
    g->rotate(m_euler_angles.z, glm::vec3(0,0,1));
    g->rotate(m_euler_angles.y, glm::vec3(0,1,0));
    g->scale(m_scale);
}

glm::vec3 TransformComponent::rotationMat2EulerAngles(const glm::mat3 &column_major_mat) {
    // intrinsic XZY euler angles
    const glm::mat3 &m = column_major_mat;
    float thetaX;
    float thetaY;
    float thetaZ;

    if (m[1][0] < 1) {
        if (m[1][0] > -1) {
            thetaZ = std::asin(-m[1][0]);
            thetaX = std::atan2(m[1][2], m[1][1]);
            thetaY = std::atan2(m[2][0], m[0][0]);
        } else { // m[1][0] == -1
            // not a unique solution
            thetaZ = M_PI / 2.f;
            thetaX = -std::atan2(-m[0][2], m[2][2]);
            thetaY = 0;
        }
    } else {
        // not a unique solution
        thetaZ = -M_PI / 2.f;
        thetaX = std::atan2(-m[0][2], m[2][2]);
        thetaY = 0;
    }

    return glm::vec3(thetaX, thetaY, thetaZ);
}
