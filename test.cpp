// 回転角に対するBER
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// 回転角
static const double theta_min = 0.0;                     // 回転角の最小値 [rad]
static const double theta_max = 45 * M_PI / 180 ;        // 回転角の最大値 [rad]
static const double theta_stp = 5.0 * M_PI / 180;        // 回転角の間隔 [rad]
double theta_deg = 0.0;                                  // csvファイルの第1カラム(角度 [°])

// ファイル
std::string filename;
std::ofstream ofs;

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
    sim.setNoiseSD(EbN0dB);

    // 最適な回転角の結果を追加（64QAM）
    // 20.9291
    theta_deg = 20.9291;
    sim.setRotationSymbol(theta_deg * 180 / M_PI);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;

    return 0;
}