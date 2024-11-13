# Dispersion in gaussian random fields

This project is a toy model for simulating dispersion of passive solute in gaussian random fields. 

## Principle

- Fokker-Planck equation describes the evolution of particle concentration:
  
  $$\frac{\partial}{\partial t} C + \nabla \cdot (\vec{u}C - D \nabla C) = 0$$

  where $C = C(\vec{x}, t)$ is concentration field, $\vec{u} = \vec{u}(\vec{x})$ is static gaussian velocity field, and $D$ is diffusion coefficient.
  
- Langevin equation describes the motion of a single particle:

  $$\frac{d}{dt} \vec{X}(t|\vec{x}_0) = \vec{u}(\vec{X}(t|\vec{x}_0)) + \vec{\xi}(t)$$

  where $\vec{X}(t|\vec{x}_0)$ is the trajetory of a particle whose initial position is $\vec{x}_0$, and $\vec{\xi}(t)$ is gaussian white noise.

- Dispersion coefficents

## Method

- Realization of gaussian random velocity fields.
  
- Runge-Kutta method for random walk particle tracking. 
