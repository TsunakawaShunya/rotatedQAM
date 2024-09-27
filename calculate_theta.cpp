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

    sim.setSymbol();            // 従来QAMシンボル設計
    sim.setNoiseSD(EbN0dB);     // Eb/N0

    // ------------ 最適な回転角を探す ------------
    // 最適化の目的関数を定義
    std::function<double(double)>func = 
    [&sim](double theta_rad) { 
        return sim.getBerUpperBoundVec(theta_rad).sum(); 
    };

    // 初期推定値と最適化の実行（ニュートン法）
    double initial_deg = 20.0;
    //double optimal_deg_byNewton = sim.get_optimizeTheta_deg_byNewton(func, initial_deg);
    //double optimal_deg = sim.get_optimizeTheta_deg(func, optimal_deg_byNewton - 1.0, optimal_deg_byNewton + 1.0, 1e-03);
    double optimal_deg = sim.get_optimizeTheta_deg(func, 20.5, 21.5, 1e-04);

    // 標準出力
    std::cout << "Optimal theta (deg): " << optimal_deg << std::endl;
    std::cout << "f(x) = : " << func(optimal_deg * 180.0 / M_PI) << std::endl;
    return 0;
}