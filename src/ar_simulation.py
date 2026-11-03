import numpy as np
import matplotlib.pyplot as plt

class AwRascleModel:
    def __init__(self, gamma=1.0):
        """
        Aw-Rascle model with pressure function p(rho) = rho^gamma
        """
        self.gamma = gamma
    
    def p(self, rho):
        """Pressure function"""
        return rho ** self.gamma
    
    def p_prime(self, rho):
        """Derivative of pressure function"""
        return self.gamma * rho ** (self.gamma - 1)
    
    def lambda1(self, rho, v):
        """First characteristic speed"""
        return v - rho * self.p_prime(rho)
    
    def lambda2(self, rho, v):
        """Second characteristic speed"""
        return v
    
    def intermediate_state(self, rho_l, v_l, rho_r, v_r):
        """Compute intermediate state (rho_m, v_m)"""
        v_m = v_r  # From contact discontinuity

        # v_m + p(rho_m) = v_l + p(rho_l) => p(rho_m) = v_l + p(rho_l) - v_m
        p_rho_m = v_l + self.p(rho_l) - v_m
        
        # Solve p(rho) = target_p for rho
        if p_rho_m <= 0:
            return 1e-6  # Small positive value to avoid division by zero

        rho_m = p_rho_m ** (1/self.gamma)
        return rho_m, v_m
    
    def shock_speed(self, rho_l, v_l, rho_m):
        """Shock speed for 1-wave"""
        return v_l - (rho_m * (self.p(rho_m) - self.p(rho_l))) / (rho_m - rho_l)
    
    def rarefaction_fan(self, rho_l, v_l, xi):
        """Solution in rarefaction fan for given xi = x/t"""
        # Solve: xi = v - rho * p'(rho) and v + p(rho) = v_l + p(rho_l)
        def equation(rho):
            v = v_l + self.p(rho_l) - self.p(rho)
            return v - rho * self.p_prime(rho) - xi
        
        # Use bisection to solve for rho
        rho_min, rho_max = 1e-6, 1.0
        for _ in range(50):
            rho_mid = (rho_min + rho_max) / 2
            if equation(rho_mid) > 0:
                rho_min = rho_mid
            else:
                rho_max = rho_mid
        
        rho = (rho_min + rho_max) / 2
        v = v_l + self.p(rho_l) - self.p(rho)
        return rho, v

def solve_riemann_aw_rascle(model, rho_l, v_l, rho_r, v_r, x, t):
    """
    Solve Riemann problem for Aw-Rascle model at given positions x and time t
    """
    # Compute intermediate state
    rho_m, v_m = model.intermediate_state(rho_l, v_l, rho_r, v_r)
    
    # Initialize solution arrays
    rho_sol = np.zeros_like(x)
    v_sol = np.zeros_like(x)
    
    # Determine wave type for 1-wave
    if rho_m < rho_l:
        # 1-Rarefaction wave
        lambda1_l = model.lambda1(rho_l, v_l)
        lambda1_m = model.lambda1(rho_m, v_m)
        lambda2_m = model.lambda2(rho_m, v_m)  # = v_m = v_r
        
        for i, xi in enumerate(x / t):
            if xi <= lambda1_l:
                # Left state
                rho_sol[i] = rho_l
                v_sol[i] = v_l
            elif lambda1_l < xi < lambda1_m:
                # Rarefaction fan
                rho_sol[i], v_sol[i] = model.rarefaction_fan(rho_l, v_l, xi)
            elif lambda1_m <= xi <= lambda2_m:
                # Intermediate state
                rho_sol[i] = rho_m
                v_sol[i] = v_m
            else:
                # Right state
                rho_sol[i] = rho_r
                v_sol[i] = v_r
    else:
        # 1-Shock wave
        s = model.shock_speed(rho_l, v_l, rho_m)
        lambda2_m = model.lambda2(rho_m, v_m)  # = v_m = v_r
        
        for i, xi in enumerate(x / t):
            if xi <= s:
                # Left state
                rho_sol[i] = rho_l
                v_sol[i] = v_l
            elif s < xi <= lambda2_m:
                # Intermediate state
                rho_sol[i] = rho_m
                v_sol[i] = v_m
            else:
                # Right state
                rho_sol[i] = rho_r
                v_sol[i] = v_r
    
    return rho_sol, v_sol

def plot_waves(model, rho_l, v_l, rho_r, v_r, case_name, times=[0.1, 0.5, 1.0, 2.0]):
    """Plot solutions at multiple time steps for a given case"""
    x = np.linspace(-2, 2, 1000)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
    
    for i, t in enumerate(times):
        rho, v = solve_riemann_aw_rascle(model, rho_l, v_l, rho_r, v_r, x, t)
        
        ax1.plot(x, rho, color=colors[i], linewidth=2, label=f't = {t}')
        ax2.plot(x, v, color=colors[i], linewidth=2, label=f't = {t}')
    
    ax1.set_xlabel('Position x')
    ax1.set_ylabel(r'Density $\rho$')
    ax1.set_title(f'{case_name}\nDensity Evolution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Position x')
    ax2.set_ylabel('Velocity v')
    ax2.set_title(f'{case_name}\nVelocity Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    # Initialize model
    model = AwRascleModel(gamma=3.0)
    
    # Case 1: Rarefaction + Contact (rho_m < rho_l)
    rho_l1, v_l1 = 0.8, 0.2
    rho_r1, v_r1 = 0.2, 0.6
    rho_m1, v_m1 = model.intermediate_state(rho_l1, v_l1, rho_r1, v_r1)
    
    # Case 2: Shock + Contact (rho_m > rho_l)
    rho_l2, v_l2 = 0.2, 0.6
    rho_r2, v_r2 = 0.8, 0.2
    rho_m2, v_m2 = model.intermediate_state(rho_l2, v_l2, rho_r2, v_r2)
    
    # Plot multiple time steps for each case
    times = [0.1, 0.5, 1.0, 2.0]
    
    print("\nPlotting multiple time steps for rarefaction case...")
    fig_rare = plot_waves(model, rho_l1, v_l1, rho_r1, v_r1, 
                            "Rarefaction + Contact Case", times)
    
    print("Plotting multiple time steps for shock case...")
    fig_shock = plot_waves(model, rho_l2, v_l2, rho_r2, v_r2, 
                            "Shock + Contact Case", times)
    
    plt.show()
