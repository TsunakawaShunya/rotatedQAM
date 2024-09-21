// 回転角に対するBER
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// 回転角
double angles_deg[3];
double angles_rad[3];
double lambda(double theta_rad) {
    return 1 / (pow(cos(theta_rad), 2) * pow(sin(theta_rad), 2)) + 1 / pow(cos(2 * theta_rad), 2);
}

int main() {
    Simulator sim;

    switch(sim.NUMBER_OF_BIT) {
        case 2:

        // optimal_angle
        angles_deg[0] = 27.4;
        angles_rad[0] = 2 * atan(sqrt(5 + 2 * sqrt(3) - 2 * sqrt(9+ 5 * sqrt(3))));

        // prev reserach's result
        angles_deg[1] = 31.7;
        angles_rad[1] = angles_deg[1] * M_PI / 180.0;

        angles_deg[2] = 29.6;
        angles_rad[2] = angles_deg[2] * M_PI / 180.0;

        for(double angle : angles_rad) {
            std::cout << angle * 180.0 / M_PI << ":" << lambda(angle) << std::endl;
        }
        break;
        case 4:
        break;
    }

    return 0;
}