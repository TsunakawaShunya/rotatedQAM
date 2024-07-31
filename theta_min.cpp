// 最小解を算出
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// 最適化の目的関数
std::function<double(double)>func;

// 初期値
double start_deg;
double end_deg;
double step_deg;
double optimal_deg;


int main() {
    Simulator sim;

    // SNR設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "EbN0dB [dB]?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> EbN0dB;

    sim.setSymbol();        // 従来QAMシンボル設計

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            sim.set_QPSKNoiseSD(EbN0dB);

            // 最適化の目的関数を定義
            func = [&sim](double theta) { return sim.getBerUpperBoundVec(theta).sum(); };

            // 初期値
            start_deg = 27.3;
            end_deg = 27.4;
            step_deg = 0.0001;
            optimal_deg = sim.get_optimizeTheta_deg(func, start_deg, end_deg, step_deg);

            // 最適な回転角を標準出力
            std::cout << "Optimal theta (deg): " << optimal_deg << std::endl;
        break;
        case 4:
            sim.set_16QAMNoiseSD(EbN0dB);

            // 最適化の目的関数を定義
            func = [&sim](double theta) { return sim.getBerUpperBoundVec(theta).sum(); };

            // 初期値
            start_deg = 20.0;
            end_deg = 22.0;
            step_deg = 0.0001;
            optimal_deg = sim.get_optimizeTheta_deg(func, start_deg, end_deg, step_deg);
            
            // 最適な回転角を標準出力
            std::cout << "Optimal theta (deg): " << optimal_deg << std::endl;
        break;
    }
    return 0;
}