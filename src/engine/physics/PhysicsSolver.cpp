#include "PhysicsSolver.h"

#include "PhysicsSystem.h"

void PhysicsSolver::eulerStep(PhysicsSystem *s, const float seconds) {
    // resolve collisions by applying linear or rotational impulse
    s->resolveCollisions();

    // compute the derivative of the state vector for all rigid bodies and apply euler step
    auto &rbs = s->getRigidBodies();
    for (auto it = rbs.begin(); it != rbs.end(); it++) {
        auto &rb = *it;
        rb->zeroOutDerivatives();
        rb->calculateDerivatives();
        rb->applyDerivatives(seconds);
        rb->zeroOutDerivatives();
    }
}
