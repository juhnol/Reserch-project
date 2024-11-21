# Flow in Gaussian random fields

This project provides a toy model for simulating flow of passive solute in gaussian random fields. 

## Principle

- Fokker-Planck equation describes the evolution of particle concentration:
  
  $$\frac{\partial}{\partial t} C + \nabla \cdot (\vec{u}C - D \nabla C) = 0$$

  where $C = C(\vec{x}, t)$ is concentration field, $\vec{u} = \vec{u}(\vec{x})$ is static gaussian velocity field, and $D$ is diffusion coefficient.
  
- Langevin equation describes the motion of a single particle:

  $$\frac{d}{dt} \vec{X}(t|\vec{x}_0) = \vec{u}(\vec{X}(t|\vec{x}_0)) + \vec{\xi}(t)$$

  where $\vec{X}(t|\vec{x}_0)$ is the trajetory of a particle whose initial position is $\vec{x}_0$, and $\vec{\xi}(t)$ is gaussian white noise.

- Dispersion coefficents

## Methods

See `flow_in_gaussian_velocity_field_v2.1.ipynb`.

- Realization of gaussian random velocity fields.
  
- Runge-Kutta method for random walk particle tracking.

- Calculation of dispersion coefficients.

## Usage

See `flow_in_gaussian_velocity_field_v2.1.ipynb`. 

- The whole simulation are based on a class `Simulation`. After specifying parameters for simulation, call `Simulation.initialize_velocity_field()` and `Simulation.run_simulation()` to run the simulation.

- After that, call `simulation.plot_particle_trajectories()` for particle trajectory visualization, and `simulation.calculate_dispersion_coefficients()` and `simulation.plot_dispersion_coefficients()` to calculate and visualize disperison coefficients.
