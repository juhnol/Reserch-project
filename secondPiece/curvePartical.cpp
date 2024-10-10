#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>

// Particle tracking
struct Particle {
    double x, y;
};

std::vector<Particle> generateGaussianParticles(int n, double mean, double stddev) {
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::normal_distribution<> dist(mean, stddev); // Define the distribution

    std::vector<Particle> particles;
    for(int i = 0; i < n; ++i) {
        Particle p;
        p.x = dist(gen); // Generate x position
        p.y = dist(gen); // Generate y position
        particles.push_back(p);
    }
    return particles;
}

int main() {
    int N = 1000; // Number of particles
    double mean = 0.0; // Mean position
    double stddev = 1.0; // Standard deviation

    auto particles = generateGaussianParticles(N, mean, stddev);

    // Open a file for output
    std::ofstream outFile("particles.csv");
    for (const auto& p : particles) {
        outFile << p.x << "," << p.y << "\n";
    }
    outFile.close();

    return 0;
}
