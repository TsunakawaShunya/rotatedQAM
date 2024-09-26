// 回転角に対するBER
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// BER
double ber;                                              // BERシミュレーション値

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

        break;
        case 4:
            sim.set_16QAMNoiseSD(EbN0dB);
        break;
    }

    // ------------------------------------最適な回転角をコンピュータサーチによって導出------------------------------------
    // 最適化の目的関数を定義
    std::function<double(double)>func = 
    [&sim](double theta) { 
        return sim.getBerUpperBoundVec(theta).sum(); 
    };

    // 初期推定値と最適化の実行（ニュートン法）
    double initial_deg = 20.0;
    double optimal_deg_byNewton = sim.get_optimizeTheta_deg_byNewton(func, initial_deg);
    double optimal_deg = sim.get_optimizeTheta_deg(func, optimal_deg_byNewton - 5.0, optimal_deg_byNewton + 5.0);

    // 標準出力
    std::cout << "Optimal theta (deg): " << optimal_deg << std::endl;
    return 0;
}