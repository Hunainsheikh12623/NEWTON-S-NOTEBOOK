#include<iostream>
#include<math.h>
#include<emscripten/emscripten.h>
using namespace std;

// 9th physics:

// chapter 1 Kinematics:

// horizontal motion:
extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double finalvelocity(double u, double a, double t) {
        return u + (a * t);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double displacement(double u, double a, double t) {
        return (u * t) + (0.5 * a * t * t);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double acceleration(double v, double u, double s) {
        double squa_diff = (v*v) - (u*u);
        return squa_diff / (2*s);
        }
}

// motion under gravity:

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double finalvelocity_gravity(double u, double t) {
        return u + (9.8 * t);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double displacement_gravity(double v) {
        return (v * v) / (2 * 9.8) ;
    }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double time_gravity(double s) {
        return sqrt(2 * s / 9.8);
    }
}


// chapter 2 Dynamics:

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double force(double m, double a) {
        return (m * a);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double mass(double f, double a) {
        return (f / a);
    }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double acceleration_force(double f, double m) {
        return (f / m);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double momentum(double m, double v) {
        return (m * v);
    }
}


// chapter 3 work, power and energy:

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double work(double f, double d) {
        return (f * d);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double power(double w, double t) {
        return (w / t);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double energy(double m, double v) {
        return (0.5 * m * v * v);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double p_energy(double m, double h) {
        return (9.8 * m * h);
        }
}


// chapter 4 gravitation:

extern "C"{
    EMSCRIPTEN_KEEPALIVE 
    double force_gravitation(double m1, int e1, double m2, int e2, double d, int e3) {
        const int Ge = -11;
        const double G = 6.674;
        int total = (e1 + e2 + Ge) - (e3*2);   
        return (G * (m1 * m2) / (d * d)) * pow(10, total);
        }
}

// for m1 and m2:
extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double m1_gravitation(double m2, int e1, double d, int e2, double f, int e3) {
        const int Ge = -11;
        const double G = 6.674;
        int total = (e3 + (2*e2)) - (e1 + Ge);

        return ((f * (d * d)) / ((G * m2))) * pow(10, total);
}
}

extern "C"{
        EMSCRIPTEN_KEEPALIVE
        double distance_gravitation(double m1, int e1, double m2, int e2, double f, int e3) {
            const int Ge = -11;
            const double G = 6.674;
            int total = (e1 + e2 + Ge) - (e3);

            return sqrt(((G * (m1 * m2)) / f) * pow(10, total));

        }
}


// class 10th:

// chapter 01: Electricity;

// ohm's law:

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double voltage(double I, double R) {
        return (I * R);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double current(double V, double R) {
        return (V / R);
        }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double resistance(double V, double I) {
        return (V / I);
    }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double electricalpower(double V, double I, double R) {
        if(R == 0) {
            return (V * I);
        }
        else if(V == 0) {
            return (I * I * R);
        }
        else if(I == 0) {
            return (V * V / R);
        }
        else {
            return 0;
        }
    }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double electricalenergy(double V, double I, double R, double t) {
        if(R == 0) {
            return (V * I * t);
        }
        else if(V == 0) {
            return (I * I * R * t);
        }
        else if(I == 0) {
            return (V * V * t / R);
        }
        else {
            return 0;
        }
    }
}


// Class 11

// chapter 01: Two Dimensional Motion:

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double hori_velocity_comp(double u, double theta) {
        return (u * cos(theta));
    }
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double vert_velocity_comp(double u, double theta) {
        return (u * sin(theta));
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double resultant_velocity(double Vx, double Vy) {
        return sqrt(pow(Vx, 2) + pow(Vy, 2));
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double horizontal_range(double u, double theta) {
        const float g = 9.8;
        return (u * u * sin(2 * theta)) / g;
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double maximum_height(double u, double theta) {
        const float g = 9.8;
        return (u * u * power(sin(theta), 2)) / (2 * g);
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double time_of_flight(double u, double theta) {
        const float g = 9.8;
        return (2 * u * sin(theta)) / g;
}
}

extern "C"{
    EMSCRIPTEN_KEEPALIVE
    double position_x_t(double u, double theta, double t) {
        const float g = 9.8;
        return (u * cos(theta) * t);
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double position_y_t(double u, double theta, double t) {
        const float g = 9.8;
        return ((u * sin(theta) * t) - (0.5 * g * t * t));
    }
}


// chapter 2, 3 and 4 sound, oscillation and waves

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double wave_speed(double frequency, double wavelength) {
        return frequency * wavelength; // v = f * λ
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double frequency_from_period(double time_period) {
        return 1.0 / time_period; // f = 1 / T
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double kinetic_energy_shm(double m, double w, double A, double X) {
        return 0.5 * m * pow(w, 2) * (pow(A, 2) - pow(X, 2)); // KE = 1/2 m ω^2 (A^2 - x^2)
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double potential_energy_shm(double m, double w, double X) {
        return 0.5 * m * pow(w, 2) * pow(X, 2); // PE = 1/2 m ω^2 x^2
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double wave_power(double md, double w, double A, double V) {
        return 0.5 * md * pow(w, 2) * pow(A, 2) * V; // P = 0.5 * μ * ω^2 * A^2 * v
    }
}

extern "C" {
    // towards the source
    EMSCRIPTEN_KEEPALIVE
    double de_Tsource(double fs, double V, double Vo) {
        return fs * (V + Vo) / V; // f' = f (v + vo) / v
    }
}

extern "C" {
    // away from the source
    EMSCRIPTEN_KEEPALIVE
    double de_Asource(double fs, double V, double Vo) {
        return fs * (V - Vo) / V; // f' = f (v - vo) / v
    }
}

extern "C" {
    // towards the observer
    EMSCRIPTEN_KEEPALIVE
    double de_Tobserver(double fs, double V, double Vs) {
        return fs * V / (V - Vs); // f' = f v / (v - vs)
    }
}

extern "C" {
    // away from the observer
    EMSCRIPTEN_KEEPALIVE
    double de_Aobserver(double fs, double V, double Vs) {
        return fs * V / (V + Vs); // f' = f v / (v + vs)
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double beat_frequency(double f1, double f2) {
        return fabs(f1 - f2); // f_beat = |f1 - f2|
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double standing_wave_wavelength(double l, int n) {
        return 2.0 * l / n; // λ_n = 2L / n
    }
}    

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double resonance_frequency(double K, double m) {
        return 1.0 / (2.0 * M_PI) * sqrt(K / m); // f = 1 / 2π √(k / m)
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double sound_intensity(double P, double A) {
        return P / A; // I = P / A
    }
}

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    double sound_level(double I, double Io) {
        return 10.0 * log10(I / Io); // L = 10 log10(I / I0)
    }
}



const double k = 8.99e9;       // Coulomb's constant
const double epsilon_0 = 8.854e-12; // Permittivity of free space
const double h = 6.626e-34;    // Planck's constant
const double c = 3e8;          // Speed of light
const double m_e = 9.109e-31;  // Mass of an electron

extern "C" {
    // Electric Field
    EMSCRIPTEN_KEEPALIVE
    double electric_field(double F, double q, double r) {
        if (q != 0 && F != 0) {
            return F / q; // E = F / q
        } else if (q != 0) {
            return k * q / (r * r); // E = k * q / r²
        }
        return 0.0;
    }

    // Electric Potential
    EMSCRIPTEN_KEEPALIVE
    double electric_potential(double q, double r) {
        return k * q / r; // V = k * q / r
    }

    // Capacitance
    EMSCRIPTEN_KEEPALIVE
    double capacitance(double Q, double V, double A, double d) {
        if (Q != 0 && V != 0) {
            return Q / V; // C = Q / V
        } else if (A != 0 && d != 0) {
            return epsilon_0 * A / d; // C = ε₀ * A / d
        }
        return 0.0;
    }

    // Energy in Capacitor
    EMSCRIPTEN_KEEPALIVE
    double energy_in_capacitor(double C, double V) {
        return 0.5 * C * V * V; // U = 1/2 * C * V²
    }

    // Energy of Photon
    EMSCRIPTEN_KEEPALIVE
    double energy_photon(double f, double lambda) {
        if (f > 0) {
            return h * f; // E = h * f
        } else if (lambda > 0) {
            return (h * c) / lambda; // E = h * c / λ
        }
        return 0.0;
    }

    // de Broglie Wavelength
    EMSCRIPTEN_KEEPALIVE
    double de_broglie_wavelength(double m, double v) {
        return h / (m * v); // λ = h / mv
    }

    // Time Dilation
    EMSCRIPTEN_KEEPALIVE
    double time_dilation(double t0, double v) {
        return t0 / sqrt(1 - (v * v) / (c * c)); // Δt = Δt₀ / √(1 - v²/c²)
    }

    // Length Contraction
    EMSCRIPTEN_KEEPALIVE
    double length_contraction(double l0, double v) {
        return l0 * sqrt(1 - (v * v) / (c * c)); // L = L₀ √(1 - v²/c²)
    }

    // Compton Wavelength Shift
    EMSCRIPTEN_KEEPALIVE
    double compton_wavelength_shift(double theta) {
        return (h / (m_e * c)) * (1 - cos(theta)); // Δλ = h/(mₑc) * (1 - cos θ)
    }

    // Nuclear Decay
    EMSCRIPTEN_KEEPALIVE
    double nuclear_decay(double N0, double lambda, double t) {
        return N0 * exp(-lambda * t); // N = N₀ e^(-λt)
    }
}

EMSCRIPTEN_KEEPALIVE
int main() {
    return 0;
}





