#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

// Particle structure

const double stddev = 0.05;

struct Vector2D {
    double x, y;
    
    Vector2D(double x = 0.0, double y = 0.0) : x(x), y(y) {}

    Vector2D operator+(const Vector2D& other) const {
        return {x + other.x, y + other.y};
    }

    Vector2D operator*(double scalar) const {
        return {x * scalar, y * scalar};
    }

    Vector2D operator-(const Vector2D& other) const {
        return {x - other.x, y - other.y};
    }
};


Vector2D calculateVelocityAtPosition(double x, double y, const std::vector<Vector2D>& waveVectors, const std::vector<double>& phases);

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

     void applyVelocityField(const std::vector<Vector2D>& waveVectors, const std::vector<double>& phases, double timeStep) {
        Vector2D velocity = calculateVelocityAtPosition(x, y, waveVectors, phases);
        x += velocity.x * timeStep;
        y += velocity.y * timeStep;
    }

};

//now I need to add the movements from the contours

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

Vector2D calculateVelocityAtPosition(double x, double y, const std::vector<Vector2D>& waveVectors, const std::vector<double>& phases) {
    Vector2D velocity(0.0, 0.0);
    for (size_t i = 0; i < waveVectors.size(); ++i) {
        Vector2D k = waveVectors[i];
        double phase = phases[i];
        double magnitude = std::sqrt(k.x * k.x + k.y * k.y);
        if (magnitude > 0.0) {
            Vector2D unitK = k * (1 / magnitude);
            Vector2D orthoK(-unitK.y, unitK.x);
            double contribution = std::cos(k.x * x + k.y * y + phase);
            velocity = velocity + orthoK * contribution;
        }
    }
    return velocity * stddev;
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
    double varianceX = sumSqX / N;
    double varianceY = sumSqY / N;

    // Output variance data to the CSV file
    if (varianceOutFile.is_open()) {
        varianceOutFile << currentTime << "," << varianceX << "," << varianceY << std::endl;
    }
}

int main() {
    const int N = 1000; // Number of particles
    const double mean = 0.0;
    const double diffusion = 0.1;
    const double timeStep = 0.05;
    const double endTime = 10.0;
    const double constantVelocityX = 1.0; // Constant velocity in the x direction
    const double checkInterval = 2.0; // Interval at which to output data
    const double pi = 3.14; 
    const int modes = 20; 
    

    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gaussian_dist(0.0, stddev);
    std::uniform_real_distribution<> uniform_dist(0.0, 2 * pi);

   



    // Generate initial set of particles
    auto particles = generateGaussianParticles(N, mean, stddev, gen);

    std::vector<Vector2D> waveVectors(N);
    std::vector<double> phases(N);
    for (int i = 0; i < modes; ++i) {
        waveVectors[i] = {gaussian_dist(gen), gaussian_dist(gen)};
        phases[i] = uniform_dist(gen);
    }

    double nextCheckTime = checkInterval; // Next time to check and output data
     std::ofstream varianceOutFile("variance_data.csv");

    for(double t = 0; t <= endTime; t += timeStep) {
        // Apply Brownian motion and velocity field effects
        for (auto& particle : particles) {
            particle.applyBrownianStep(std::sqrt(2 * diffusion * timeStep), timeStep, diffusion, constantVelocityX, gen);
            particle.applyVelocityField(waveVectors, phases, timeStep);
        }

        if (t >= nextCheckTime || std::abs(t - endTime) < 1e-6) {
            std::ostringstream filename;
            filename << "particlesWithBrownian" << std::fixed << std::setprecision(2) << t << ".csv";
            std::ofstream outFile(filename.str());

            if (outFile.is_open()) {
                for (const auto& particle : particles) {
                    outFile << particle.x << "," << particle.y << "\n";
                }
                outFile.close();
            } else {
                std::cerr << "Unable to open file " << filename.str() << " for writing." << std::endl;
            }

            // Calculate and output variance data for the current timestep
            calculateAndPrintStats(particles, N, t, varianceOutFile);

            nextCheckTime += checkInterval;
        }
    }

    varianceOutFile.close(); // Close the variance data file after the loop
    return 0;
}

