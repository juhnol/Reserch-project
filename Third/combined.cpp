#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

// Particle structure
struct Particle {
    double x, y;

    // Apply a Brownian step to the particle, including constant velocity effect
    void applyBrownianStep(double stepStdDev, double timeStep, double diffusion, double constantVelocityX, std::mt19937& gen) {
        std::normal_distribution<> stepDist(0.0, stepStdDev);
        // Apply step in x direction with Brownian motion and constant velocity
        x += stepDist(gen) * std::sqrt(2 * diffusion * timeStep) + constantVelocityX * timeStep; 
        // Apply step in y direction with Brownian motion
        y += stepDist(gen) * std::sqrt(2 * diffusion * timeStep);
    }
};

//now I need to add the movements from the contours
//plot this data 

std::vector<Particle> generateGaussianParticles(int n, double mean, double stddev, std::mt19937& gen) {
    std::normal_distribution<> dist(mean, stddev);
    std::vector<Particle> particles;
    for (int i = 0; i < n; ++i) {
        Particle p;
        p.x = dist(gen);
        p.y = dist(gen);
        particles.push_back(p);
    }
    return particles;
}

void calculateAndPrintStats(const std::vector<Particle>& particles, int N, double currentTime, std::ofstream& varianceOutFile) {
    double sumX = 0, sumY = 0;
    for (const auto& p : particles) {
        sumX += p.x;
        sumY += p.y;
    }
    double avgX = sumX / N;
    double avgY = sumY / N;

    double sumSqX = 0, sumSqY = 0;
    for (const auto& p : particles) {
        sumSqX += (p.x - avgX) * (p.x - avgX);
        sumSqY += (p.y - avgY) * (p.y - avgY);
    }
    double stdDevX = std::sqrt(sumSqX / N);
    double stdDevY = std::sqrt(sumSqY / N);

    double varianceX = sumSqX / N;
    double varianceY = sumSqY / N;

    if (varianceOutFile.is_open()) {
        varianceOutFile << currentTime << "," << varianceX << "," << varianceY << std::endl;
    }

    std::cout << "Average X: " << avgX << ", Average Y: " << avgY << std::endl;
    std::cout << "StdDev X: " << stdDevX << ", StdDev Y: " << stdDevY << std::endl;
}

int main() {
    const int N = 1000; // Number of particles
    const double mean = 0.0;
    const double stddev = 0.05;
    const double diffusion = 0.1;
    const double timeStep = 0.05;
    const double endTime = 10.0;
    const double constantVelocityX = 1.0; // Constant velocity in the x direction
    const double checkInterval = 2.0; // Interval at which to output data

    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate initial set of particles
    auto particles = generateGaussianParticles(N, mean, stddev, gen);

    double nextCheckTime = checkInterval; // Next time to check and output data
    std::ofstream varianceOutFile("variance_data1.csv"); // New variance output file name

    for(double t = 0; t <= endTime; t += timeStep) {
        // Apply Brownian motion step to each particle
        for (auto& particle : particles) {
            particle.applyBrownianStep(std::sqrt(2 * diffusion * timeStep), timeStep, diffusion, constantVelocityX, gen);
        }

        // Check if it's time to output data
        if (t >= nextCheckTime || std::abs(t - endTime) < 1e-6) {
            // Output data to CSV file
            std::ostringstream filename;
            filename << "particlesWithBrownian1" << std::fixed << std::setprecision(2) << t << ".csv";
            std::ofstream outFile(filename.str());

            if (outFile.is_open()) {
                for (const auto& particle : particles) {
                    outFile << particle.x << "," << particle.y << "\n";
                }
                outFile.close();
            } else {
                std::cerr << "Unable to open file " << filename.str() << " for writing." << std::endl;
            }

            std::cout << "Time = " << t << std::endl;
            calculateAndPrintStats(particles, N, t, varianceOutFile);

            nextCheckTime += checkInterval; // Schedule next output
        }
    }
    varianceOutFile.close();

    return 0;

}