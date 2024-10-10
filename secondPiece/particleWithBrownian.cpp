#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream> // Make sure this include is present for std::ostringstream
#include <iomanip> // For std::setprecision

// Particle tracking
struct Particle {
    double x, y;

    // Apply a Brownian step to the particle
    void applyBrownianStep(double stepStdDev, double timeStep, double diffusion, std::mt19937& gen) {
        std::normal_distribution<> stepDist(0.0, stepStdDev);
        x += stepDist(gen) * std::sqrt(2 * diffusion * timeStep); // Apply step in x direction
        y += stepDist(gen) * std::sqrt(2 * diffusion * timeStep); // Apply step in y direction
    }
};

std::vector<Particle> generateGaussianParticles(int n, double mean, double stddev, std::mt19937& gen) {
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

void calculateAndPrintStats(const std::vector<Particle>& particles, int N) {
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

    std::cout << "Average X: " << avgX << ", Average Y: " << avgY << std::endl;
    std::cout << "StdDev X: " << stdDevX << ", StdDev Y: " << stdDevY << std::endl;
}

int main() {
    int N = 1000;
    double mean = 0.0;
    double stddev = 0.05;
    double diffusion = 0.1;
    double timeStep = 0.05;
    double endTime = 10.0;
    double stepStdDev = std::sqrt(2 * diffusion * timeStep);
    double checkInterval = 2.0;

    std::random_device rd;  
    std::mt19937 gen(rd());

    auto particles = generateGaussianParticles(N, mean, stddev, gen);


    double nextCheckTime = checkInterval;
    for(double t = 0; t <= endTime; t += timeStep) {
       
        for (auto& p : particles) {
            p.applyBrownianStep(stepStdDev, timeStep, diffusion, gen);
        }

        if (t >= nextCheckTime || std::abs(t - endTime) < 1e-6) {
            double stdDevTheoretical = stepStdDev; 
            std::normal_distribution<> theoreticalDist(0.0, stdDevTheoretical);

            std::vector<double> theoreticalX(N), theoreticalY(N);
            for(int i = 0; i < N; ++i) {
                theoreticalX[i] = theoreticalDist(gen);
                theoreticalY[i] = theoreticalDist(gen);
            }

            std::ostringstream theoreticalFileName;
            theoreticalFileName << "theoreticalParticlesWithBrownian" << std::fixed << std::setprecision(2) << t << ".csv";
            std::ofstream theoreticalOutFile(theoreticalFileName.str());
             for(int i = 0; i < N; ++i) {
                theoreticalOutFile << theoreticalX[i] << "," << theoreticalY[i] << "\n";
            }
            theoreticalOutFile.close();


            std::cout << "Time = " << t << std::endl;
            calculateAndPrintStats(particles, N);

            // CSV file output
            std::ostringstream fileName;
            fileName << "particlesWithBrownian" << std::fixed << std::setprecision(2) << t << ".csv";
            std::ofstream outFile(fileName.str());
            for (const auto& p : particles) {
                outFile << p.x << "," << p.y << "\n";
            }
            outFile.close();

            nextCheckTime += checkInterval;
        }
    }
  

    return 0;
}




