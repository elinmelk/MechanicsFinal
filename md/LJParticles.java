/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.sip.ch08.md; // Ensure this package matches your directory structure
import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;

import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;

/**
 * LJParticles evolves a two-dimensional system of interacting particles
 * via the Lennard-Jones potential using a Verlet ODESolver.
 * Can now use repulsive walls for confinement instead of PBC.
 * Also draws a Minimum Bounding Box (MBB) around all particles.
 *
 * @author Jan Tobochnik, Wolfgang Christian, Harvey Gould
 * @version 1.4 revised 05/15/2025 (added MBB drawing)
 */
public class LJParticles implements Drawable, ODE {
  public double state[];
  public double ax[], ay[];
  public int N, nx, ny;
  public double Lx, Ly;
  public double rho;
  public double initialKineticEnergy;
  public int steps = 0;
  public double dt = 0.01;
  public double t;
  public double totalPotentialEnergyAccumulator;
  public double totalKineticEnergyAccumulator, totalKineticEnergySquaredAccumulator;
  public double virialAccumulator;
  public String initialConfiguration;
  Verlet odeSolver = new Verlet(this);
  public double[] radii = new double[0];
  HashMap<Point, ArrayList<Integer>> grid; // Unused currently
  public double cellSize = 2.0; // Unused currently

  // Central force and damping parameters from user's version
  public double gravityForce = 20.89;
  public double dampingCoeff = 220;
  public double springConstant = 30;

  // Parameters for Wall Forces
  public boolean useRepulsiveWalls = false;
  public double wallStiffness = 200.0;

  // --- NEW: Fields for Minimum Bounding Box ---
  public double mbbMinX, mbbMaxX, mbbMinY, mbbMaxY;

  public void initialize() {
    if (nx <= 0 || ny <= 0) { N = 0; } else { N = nx*ny; }
    if (Lx == 0 || Ly == 0) {
        rho = (N > 0) ? Double.NaN : 0;
        if (N > 0) System.err.println("Lx or Ly is zero. Density calculation will fail.");
    } else {
        rho = N / (Lx * Ly);
    }

    t = 0;
    resetAverages(); // This will also reset MBB
    state = new double[1 + 4 * N];
    ax = new double[N];
    ay = new double[N];

    if (N > 0) {
        radii = new double[N];
        for (int i = 0; i < N; i++) {
          radii[i] = 0.3 + 0.4 * Math.random();
        }
    } else {
        radii = new double[0];
    }

    if (initialConfiguration == null) {
        initialConfiguration = "random";
        System.err.println("Warning: initialConfiguration is null. Defaulting to random.");
    }

    if ("triangular".equals(initialConfiguration)) setTriangularLattice();
    else if ("rectangular".equals(initialConfiguration)) setRectangularLattice();
    else setRandomPositions();
    
    setVelocities();
    computeAcceleration();
    odeSolver.setStepSize(dt);
    updateMinimumBoundingBox(); // Initial MBB calculation
  }

  // --- NEW: Method to update MBB, also called from resetAverages ---
  private void updateMinimumBoundingBox() {
      if (N > 0 && state != null && radii != null && state.length >= 4*N && radii.length == N) {
          mbbMinX = Double.MAX_VALUE;
          mbbMaxX = -Double.MAX_VALUE; // or Double.MIN_VALUE for positive numbers
          mbbMinY = Double.MAX_VALUE;
          mbbMaxY = -Double.MAX_VALUE;

          for (int i = 0; i < N; i++) {
              double x = state[4*i];
              double y = state[4*i+2];
              double r = radii[i];
              mbbMinX = Math.min(mbbMinX, x - r);
              mbbMaxX = Math.max(mbbMaxX, x + r);
              mbbMinY = Math.min(mbbMinY, y - r);
              mbbMaxY = Math.max(mbbMaxY, y + r);
          }
      } else { // No particles or not ready, reset MBB to box boundaries or zero size
          mbbMinX = 0; mbbMaxX = Lx; // Default to full simulation box
          mbbMinY = 0; mbbMaxY = Ly;
          if (N==0) { // Or make it zero size if no particles
            mbbMinX = Lx / 2.0; mbbMaxX = Lx / 2.0;
            mbbMinY = Ly / 2.0; mbbMaxY = Ly / 2.0;
          }
      }
  }
  
  // setRandomPositions, setRectangularLattice, setTriangularLattice, setVelocities as before
  // (ensure setRandomPositions uses direct distances for overlap check if useRepulsiveWalls is true, as in previous response)

  public void setRandomPositions() {
    if (N == 0) return;
    boolean overlap;
      for (int i = 0; i < N; ++i) {
          int attempts = 0;
          do {
              overlap = false;
              state[4 * i] = Lx * Math.random();
              state[4 * i + 2] = Ly * Math.random();
              for (int j = 0; j < i; ++j) {
                  double dx_val, dy_val;
                  if (useRepulsiveWalls) {
                      dx_val = state[4*i] - state[4*j];
                      dy_val = state[4*i+2] - state[4*j+2];
                  } else {
                      dx_val = pbcSeparation(state[4 * i] - state[4 * j], Lx);
                      dy_val = pbcSeparation(state[4 * i + 2] - state[4 * j + 2], Ly);
                  }
                  double rMin = radii[i] + radii[j];
                  if (dx_val * dx_val + dy_val * dy_val < rMin * rMin) {
                      overlap = true;
                      break;
                  }
              }
              attempts++;
              if (attempts > 20000 && N > 1) {
                  System.err.println("Warning: Could not place particle " + i + " without overlap after many attempts.");
                  break;
              }
          } while (overlap);
      }
  }

  public void setRectangularLattice() {
    if (N == 0 || nx == 0 || ny == 0) return;
    double dx_spacing = Lx/nx;
    double dy_spacing = Ly/ny;
    for(int ix = 0;ix<nx;++ix) {
      for(int iy = 0;iy<ny;++iy) {
        int i = ix+iy*nx;
        if (i < N) {
            state[4*i] = dx_spacing*(ix+0.5);
            state[4*i+2] = dy_spacing*(iy+0.5);
        }
      }
    }
  }

  public void setTriangularLattice() {
    if (N == 0 || nx == 0 || ny == 0) return;
    double dx_spacing = Lx/nx;
    double dy_spacing = Ly/ny;
    for(int ix = 0;ix<nx;++ix) {
      for(int iy = 0;iy<ny;++iy) {
        int i = ix+iy*nx;
        if (i < N) {
            state[4*i+2] = dy_spacing*(iy+0.5);
            if(iy%2==0) {
              state[4*i] = dx_spacing*(ix+0.25);
            } else {
              state[4*i] = dx_spacing*(ix+0.75);
            }
        }
      }
    }
  }

  public void setVelocities() {
    if (N == 0) return;
    double vxSum = 0.0;
    double vySum = 0.0;
    for(int i = 0;i<N;++i) {
      state[4*i+1] = Math.random()-0.5;
      state[4*i+3] = Math.random()-0.5;
      vxSum += state[4*i+1];
      vySum += state[4*i+3];
    }

    double vxcm = (N==0) ? 0 : vxSum/N;
    double vycm = (N==0) ? 0 : vySum/N;
    for(int i = 0;i<N;++i) {
      state[4*i+1] -= vxcm;
      state[4*i+3] -= vycm;
    }
    double v2sum = 0;
    for(int i = 0;i<N;++i) {
      v2sum += state[4*i+1]*state[4*i+1]+state[4*i+3]*state[4*i+3];
    }
     if (v2sum == 0 && initialKineticEnergy > 0 && N > 0) {
        System.err.println("Warning: Kinetic energy is zero before rescaling. Re-randomizing velocities slightly.");
        for(int i = 0; i < N; ++i) {
            state[4*i+1] += (Math.random()-0.5) * 1e-6;
            state[4*i+3] += (Math.random()-0.5) * 1e-6;
        }
        v2sum = 0;
        for(int i = 0;i<N;++i) v2sum += state[4*i+1]*state[4*i+1]+state[4*i+3]*state[4*i+3];
    }
    if (N == 0 || v2sum == 0) { return; }

    double kineticEnergyPerParticle = 0.5*v2sum/N;
    double rescale = (kineticEnergyPerParticle == 0) ? 1.0 : Math.sqrt(initialKineticEnergy/kineticEnergyPerParticle);
    if (Double.isInfinite(rescale) || Double.isNaN(rescale)) rescale = 1.0;

    for(int i = 0;i<N;++i) {
      state[4*i+1] *= rescale;
      state[4*i+3] *= rescale;
    }
  }


  public double getMeanTemperature() {
    if (N == 0 || steps == 0) return 0;
    return totalKineticEnergyAccumulator/(N*steps);
  }

  public double getMeanEnergy() {
    if (steps == 0) return 0;
    return totalKineticEnergyAccumulator/steps+totalPotentialEnergyAccumulator/steps;
  }

  public double getMeanPressure() {
    if (N == 0 || steps == 0) return 0;
    double meanTemp = getMeanTemperature();
    if (meanTemp == 0 && N > 0 && !useRepulsiveWalls) return Double.POSITIVE_INFINITY;
    if (meanTemp == 0 && N == 0) return 1.0;
    if (meanTemp == 0 && N > 0 && useRepulsiveWalls) return 0; // Or calculate from wall collisions
    double meanVirial = virialAccumulator/steps;
    // For a walled system, the virial pressure definition needs careful handling of wall forces.
    // The current virialAccumulator only includes inter-particle forces.
    // A more general form for pressure P = rho * T + Virial / (Dim * Volume)
    return rho * meanTemp + meanVirial / (2.0 * Lx * Ly); // Assuming 2 dimensions
  }

  public double getHeatCapacity() {
     if (N == 0 || steps == 0) return 0;
     double meanTemperature = getMeanTemperature();
     if (meanTemperature == 0) return 0;
     double meanKineticEnergySquared = totalKineticEnergySquaredAccumulator/steps;
     double meanKineticEnergy = totalKineticEnergyAccumulator/steps;
     double sigma2 = meanKineticEnergySquared-meanKineticEnergy*meanKineticEnergy; // Var(Total KE)
     // Cv/kB = N_dof / (2 - (4/N_dof) * Var(KE_total)/<KE_total>^2) -- for NVT
     // The formula used: N_particles / (1 - Var(Total KE) / (N_particles * T^2))
     double denom = 1.0 - sigma2/(N * meanTemperature * meanTemperature); // Assuming k_B = 1
     if (denom == 0) return (sigma2 == 0 && N>0) ? N : Double.POSITIVE_INFINITY;
     return N/denom;
  }

  public void resetAverages() {
    steps = 0;
    virialAccumulator = 0;
    totalPotentialEnergyAccumulator = 0;
    totalKineticEnergyAccumulator = 0;
    totalKineticEnergySquaredAccumulator = 0;
    updateMinimumBoundingBox(); // Reset/recalculate MBB when averages are reset
  }

  public void computeAcceleration() {
    // ... (as in previous response, with conditional LJ separation and wall/central forces) ...
    // This method remains largely the same as the one you provided in the prior prompt,
    // which correctly implements the switch for useRepulsiveWalls.
    if (N == 0) return;
    for(int i = 0;i<N;i++) {
      ax[i] = 0;
      ay[i] = 0;
    }

    for(int i = 0;i<N-1;i++) {
      for(int j = i+1;j<N;j++) {
        double dx_lj, dy_lj;
        if (useRepulsiveWalls) { 
            dx_lj = state[4*i] - state[4*j];
            dy_lj = state[4*i+2] - state[4*j+2];
        } else { 
            dx_lj = pbcSeparation(state[4*i]-state[4*j], Lx);
            dy_lj = pbcSeparation(state[4*i+2]-state[4*j+2], Ly);
        }
        double r2 = dx_lj*dx_lj + dy_lj*dy_lj;
        if (r2 == 0) continue;
        if (r2 < 0.0001 && N > 1) r2 = 0.0001; 
        double oneOverR2 = 1.0/r2;
        double oneOverR6 = oneOverR2*oneOverR2*oneOverR2;
        double fOverR = 48.0*oneOverR6*(oneOverR6-0.5)*oneOverR2;
        double fx = fOverR*dx_lj;
        double fy = fOverR*dy_lj;
        ax[i] += fx; ay[i] += fy;
        ax[j] -= fx; ay[j] -= fy;
        totalPotentialEnergyAccumulator += 4.0*(oneOverR6*oneOverR6-oneOverR6);
        virialAccumulator += dx_lj*fx + dy_lj*fy;
      }
    }

    for(int i = 0; i < N; i++) {
        if (useRepulsiveWalls) {
            double x = state[4*i]; double y = state[4*i+2];
            double R_particle = radii[i];
            if (x < R_particle) ax[i] += wallStiffness * (R_particle - x);
            if (x > Lx - R_particle) ax[i] -= wallStiffness * (x - (Lx - R_particle));
            if (y < R_particle) ay[i] += wallStiffness * (R_particle - y);
            if (y > Ly - R_particle) ay[i] -= wallStiffness * (y - (Ly - R_particle));
        } else {
            double centerX = Lx / 2.0; double centerY = Ly / 2.0;
            double dx_center = pbcSeparation(state[4*i] - centerX, Lx);
            double dy_center = pbcSeparation(state[4*i+2] - centerY, Ly);
            double effectiveK = springConstant - gravityForce; 
            ax[i] += -effectiveK * dx_center;
            ay[i] += -effectiveK * dy_center;
        }
        ax[i] += -dampingCoeff * state[4*i + 1]; 
        ay[i] += -dampingCoeff * state[4*i + 3]; 
    }
  }
  
  // buildSpatialGrid, pbcSeparation, pbcPosition as before

  private void buildSpatialGrid() { 
      grid = new HashMap<>();
      for (int i = 0; i < N; i++) {
          int cellX = (int)(state[4 * i] / cellSize);
          int cellY = (int)(state[4 * i + 2] / cellSize);
          Point cell = new Point(cellX, cellY);
          grid.computeIfAbsent(cell, k -> new ArrayList<>()).add(i);
      }
  }

  private double pbcSeparation(double ds, double L) {
    if (L == 0) return ds; 
    if (ds > L/2.0) ds -= L; 
    else if (ds < -L/2.0) ds += L;
    return ds;
  }

  private double pbcPosition(double s, double L) {
    if (L == 0) return s; 
    if (s >= L) s -= L; 
    else if (s < 0) s += L; 
    return s;
  }


  public void getRate(double[] state, double[] rate) {
    // ... (as in previous response) ...
    if (N == 0) { 
        if (rate.length > 0 && rate.length == state.length && state.length == 1+4*N && 4*N < rate.length) rate[4*N] = 1; 
        return;
    }
    if(odeSolver.getRateCounter()==1) {
      computeAcceleration();
    }
    for(int i = 0;i<N;i++) {
      rate[4*i] = state[4*i+1];
      rate[4*i+2] = state[4*i+3];
      rate[4*i+1] = ax[i];
      rate[4*i+3] = ay[i];
    }
    if (4*N < rate.length) rate[4*N] = 1; 
  }

  public double[] getState() {
    return state;
  }

  public void step(HistogramFrame xVelocityHistogram) {
    // ... (check for N==0, null state) ...
    if (state == null || state.length == 0 ) return;
    if (N == 0 && state.length == 1 ) { 
        if (state.length > 0) state[0] += dt; 
        t = (state.length > 0) ? state[0] : t + dt;
        steps++;
        return;
    }
    if (N == 0) return;


    odeSolver.step();
    double currentStepKineticEnergy = 0;
    for(int i = 0;i<N;i++) {
      currentStepKineticEnergy += (state[4*i+1]*state[4*i+1]+state[4*i+3]*state[4*i+3]);
      if (xVelocityHistogram != null) {
          xVelocityHistogram.append(state[4*i+1]);
      }
      
      if (!useRepulsiveWalls) {
          state[4*i] = pbcPosition(state[4*i], Lx);
          state[4*i+2] = pbcPosition(state[4*i+2], Ly);
      }
    }
    currentStepKineticEnergy *= 0.5;
    steps++;
    totalKineticEnergyAccumulator += currentStepKineticEnergy;
    totalKineticEnergySquaredAccumulator += currentStepKineticEnergy*currentStepKineticEnergy;
    
    if (4*N < state.length) t = state[4*N];
    else t += dt;

    // --- NEW: Update MBB after positions are updated ---
    updateMinimumBoundingBox();
  }

  public void draw(DrawingPanel panel, Graphics g) {
	  if (state == null || radii == null || N == 0) {
          return;
      }

      // Draw particles as black outlines
      g.setColor(Color.black);
      for (int i = 0; i < N; i++) {
          double particleRadius = radii[i];
          int pxRadius = Math.abs(panel.xToPix(particleRadius) - panel.xToPix(0));
          int pyRadius = Math.abs(panel.yToPix(particleRadius) - panel.yToPix(0)); // Can be different for non-square aspect
          int xCenterPix = panel.xToPix(state[4 * i]);
          int yCenterPix = panel.yToPix(state[4 * i + 2]);
          int xpix = xCenterPix - pxRadius;
          int ypix = yCenterPix - pyRadius;
          g.drawOval(xpix, ypix, 2 * pxRadius, 2 * pyRadius);
      }

      // Draw the main simulation box boundary
      g.setColor(Color.black);
      int xBoxPix = panel.xToPix(0);
      int yBoxPix = panel.yToPix(Ly);
      int wBoxPix = panel.xToPix(Lx) - panel.xToPix(0);
      int hBoxPix = panel.yToPix(0) - panel.yToPix(Ly);
      g.drawRect(xBoxPix, yBoxPix, wBoxPix, hBoxPix);

      // --- NEW: Draw the Minimum Bounding Box for particles ---
      if (N > 0 && mbbMaxX > mbbMinX && mbbMaxY > mbbMinY) { // Check if MBB is valid
          g.setColor(Color.green); // Distinct color for the MBB

          // Convert MBB world coordinates to pixel coordinates
          int p_mbbMinX = panel.xToPix(mbbMinX);
          int p_mbbMaxY = panel.yToPix(mbbMaxY); // Top-left y in pixels
          int p_mbbMaxX = panel.xToPix(mbbMaxX);
          int p_mbbMinY = panel.yToPix(mbbMinY); // Bottom-right y in pixels (usually larger value for screen y)

          int mbbPixelWidth = p_mbbMaxX - p_mbbMinX;
          int mbbPixelHeight = p_mbbMinY - p_mbbMaxY; // This should be positive for typical screen coordinates

          if (mbbPixelWidth > 0 && mbbPixelHeight > 0) {
             g.drawRect(p_mbbMinX, p_mbbMaxY, mbbPixelWidth, mbbPixelHeight);
          }
      }
  }
}
// LICENSE TEXT