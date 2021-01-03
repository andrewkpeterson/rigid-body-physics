#include "TransformComponent.h"

TransformComponent::TransformComponent(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale) :
    m_position(position), m_rotation(rotation), m_scale(scale)
{
    g = Graphics::getGlobalInstance();
}

void TransformComponent::setTransform() {
    g->translate(m_position);
    g->rotate(m_rotation.x, glm::vec3(1,0,0));
    g->rotate(m_rotation.z, glm::vec3(0,0,1));
    g->rotate(m_rotation.y, glm::vec3(0,1,0));
    g->scale(m_scale);
}
