// SNに対するBER
/*  
    最適な回転角
    QPSK: 27.3678°
    16QAM: 21.0148°
    64QAM: 20.7135°
    256QAM: 20.6784°
*/

#include "simulator.h"

// SNR
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 40.1;       // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 5.0;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// 回転角
double theta_deg;
double theta;

// ファイル
std::string filename;
std::ofstream ofs;

// BER
double ber;                                 // BERシミュレーション値

int main() {   
    Simulator sim;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "SN - BER" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "theta? [deg]" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> theta_deg;
    theta = theta_deg * M_PI / 180;

    // シンボル設計
    sim.setSymbol();                        // 従来QAMでのシンボル設計
    sim.setRotationSymbol(theta);           // 回転

    // ファイルの初期化
    int M = std::pow(2, sim.NUMBER_OF_BIT);     // 多値数
    filename = "Rotated" + std::to_string(M) + "QAM_sim_" + std::to_string(theta_deg) + "deg.csv";
    ofs.open(filename);

    // Eb/N0[dB]でループ
    for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        sim.setNoiseSD(EbN0dB);

        // シミュレーション
        ber = sim.getBerSimulation();

        // 標準出力
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

        // ファイル出力
        ofs << EbN0dB << "," << ber << std::endl;
    }

    ofs.close();

    return 0;
}
