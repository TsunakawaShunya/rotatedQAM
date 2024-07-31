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
std::string filenameUpperBoundSum;      // 厳密な BER 上界の合計
std::vector<std::string> filenameUpperBound;        // 誤りビットごとの厳密な BER 上界の合計
std::string filenameNearlyUpperBoundSum;        // BER 上界近似の合計
std::vector<std::string> filenameNearlyUpperBound;      // 誤りビット数ごとの BER 上界近似

std::ofstream ofs;
std::ofstream ofsUpperBoundSum;     // 厳密な BER 上界の合計
std::vector<std::ofstream> ofsUpperBound;      // 誤りビット数ごとの 厳密な BER 上界
std::ofstream ofsNearlyUpperBoundSum;        // BER 上界近似の合計
std::vector<std::ofstream> ofsNearlyUpperBound;      // 誤りビット数ごとのBER 上界近似

// BER
double ber;     // BERシミュレーション値
Eigen::VectorXd berUpperBoundVec;       // BER上界
Eigen::VectorXd berNearlyUpperBoundVec;     // BER上界の近似値

int main() {
    Simulator sim;

    // SNR設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "EbN0dB [dB]?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> EbN0dB;

    sim.setSymbol();        // 従来QAMシンボル設計

    // 配列の初期化
    filenameUpperBound.resize(sim.NUMBER_OF_BIT);
    ofsUpperBound.resize(sim.NUMBER_OF_BIT);
    berUpperBoundVec.resize(sim.NUMBER_OF_BIT);
    filenameNearlyUpperBound.resize(sim.NUMBER_OF_BIT);
    ofsNearlyUpperBound.resize(sim.NUMBER_OF_BIT);
    berNearlyUpperBoundVec.resize(sim.NUMBER_OF_BIT);

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            /*        
            // 最適な回転角の結果を追加
            // 27.3678
            theta_deg = 27.3678;
            //std::cout << theta_opt * 180 / M_PI << std::endl;
            sim.setRotationSymbol(theta_opt);
            ber = sim.getBerSimulation();
            std::cout << theta_deg << "," << ber << std::endl;
            ofs << theta_deg << "," << ber << std::endl;

            theta_deg = 29.6;
            sim.setRotationSymbol(theta_deg * M_PI / 180);
            ber = sim.getBerSimulation();
            std::cout << theta_deg << "," << ber << std::endl;
            ofs << theta_deg << "," << ber << std::endl;

            theta_deg = 31.7;
            sim.setRotationSymbol(theta_deg * M_PI / 180);
            ber = sim.getBerSimulation();
            std::cout << theta_deg << "," << ber << std::endl;
            ofs << theta_deg << "," << ber << std::endl;
            ofs.close();
            */
        break;
        case 4:
            sim.set_16QAMNoiseSD(EbN0dB);

            double optimal_deg;
            optimal_deg = 21.0148;

            sim.setRotationSymbol(optimal_deg * M_PI / 180);       // 回転 

            ber = sim.getBerSimulation();

            // 標準出力
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << "simulation : " << optimal_deg << "," << ber << std::endl;
        break;
    }

    return 0;
}