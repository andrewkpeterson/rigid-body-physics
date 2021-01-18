#ifndef PHYSICSSOLVER_H
#define PHYSICSSOLVER_H

#include "src/engine/util/CommonIncludes.h"

class PhysicsSystem;

class PhysicsSolver
{
public:
    static void eulerStep(PhysicsSystem *s, const float seconds);
    static void midpointMethodStep(PhysicsSystem *s, const float seconds);
};

#endif // PHYSICSSOLVER_H
