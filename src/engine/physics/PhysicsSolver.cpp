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

/*
void PhysicsSolver::midpointMethodStep(PhysicsSystem *s, const float seconds) {

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


    // save the state of the system
    std::vector<Vector3f> saved_pos;
    std::vector<Vector3f> saved_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        saved_pos.push_back(node->pos);
        saved_vel.push_back(node->vel);
    }

    // move halfway in the direction of the derivative
    system->accumulateForces(ground_y, spheres);
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        //std::cout << node->force_acc << std::endl;
        node->pos += (node->vel * dt) / 2.0f;
        node->vel += ((node->force_acc / node->mass) * dt) / 2.0f;
    }
    system->zeroOutForces();

    // evaluate the forces (we already have the velocity) at the midpoint
    system->accumulateForces(ground_y, spheres);

    // return back to the state at the beginning of the timestep and take a full step in direction of derivative at midpoint
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        //std::cout << node->force_acc << std::endl;
        node->pos = saved_pos[i] + (node->vel * dt); // node->vel is the velocity at the midpoint
        node->vel = saved_vel[i] + ((node->force_acc / node->mass) * dt); // node->force_acc is the force at the midpoint
    }

    system->zeroOutForces();
}
*/
/*
void Solver::RK4Step(std::shared_ptr<System> system, float dt, float ground_y,
                     std::vector<std::shared_ptr<Sphere> > spheres) {

    const std::vector<std::shared_ptr<Node>> nodes = system->getNodes();

    // save the state of the system
    std::vector<Vector3f> saved_pos;
    std::vector<Vector3f> saved_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        saved_pos.push_back(node->pos);
        saved_vel.push_back(node->vel);
    }

    // move halfway in the direction of the derivative and save the increment k1
    system->accumulateForces(ground_y, spheres);
    std::vector<Vector3f> k1_pos;
    std::vector<Vector3f> k1_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        k1_pos.push_back(node->vel * dt);
        k1_vel.push_back((node->force_acc / node->mass) * dt);
        node->pos += (node->vel * dt) / 2.0f;
        node->vel += ((node->force_acc / node->mass) * dt) / 2.0f;
    }
    system->zeroOutForces();

    // evaluate the derivative at the midpoint calculated with k1 and save the increment k2
    system->accumulateForces(ground_y, spheres);
    std::vector<Vector3f> k2_pos;
    std::vector<Vector3f> k2_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        k2_pos.push_back(node->vel * dt); // node->vel is velocity at midpoint using k1
        k2_vel.push_back((node->force_acc / node->mass) * dt);
        node->pos = saved_pos[i] + (node->vel * dt) / 2.0f;
        node->vel = saved_vel[i] + ((node->force_acc / node->mass) * dt) / 2.0f;
    }
    system->zeroOutForces();

    // evaluate the derivative at the midpoint calculated with k2 and save the increment k3
    system->accumulateForces(ground_y, spheres);
    std::vector<Vector3f> k3_pos;
    std::vector<Vector3f> k3_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        k3_pos.push_back(node->vel * dt); // node->vel is velocity at midpoint using k2
        k3_vel.push_back((node->force_acc / node->mass) * dt);
        node->pos = saved_pos[i] + (node->vel * dt);
        node->vel = saved_vel[i] + ((node->force_acc / node->mass) * dt);
    }
    system->zeroOutForces();

    // evaluate the derivative at the endpoint calculated with k3 and save the increment k4
    system->accumulateForces(ground_y, spheres);
    std::vector<Vector3f> k4_pos;
    std::vector<Vector3f> k4_vel;
    for (int i = 0; i < nodes.size(); i++) {
        std::shared_ptr<Node> node = nodes[i];
        k4_pos.push_back(node->vel * dt); // node->vel is velocity at midpoint using k3
        k4_vel.push_back((node->force_acc / node->mass) * dt);
        node->pos = saved_pos[i] + 1/6.0f * (k1_pos[i] + 2.0*k2_pos[i] + 2.0*k3_pos[i] + k4_pos[i]);
        node->vel = saved_vel[i] + 1/6.0f * (k1_vel[i] + 2.0*k2_vel[i] + 2.0*k3_vel[i] + k4_vel[i]);
    }
    system->zeroOutForces();
}
*/
