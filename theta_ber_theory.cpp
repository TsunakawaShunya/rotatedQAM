// 回転角に対するBER
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// 回転角
static const double theta_min = 0.0;                    // 回転角の最小値 [rad]
static const double theta_max = 45 * M_PI / 180 ;       // 回転角の最大値 [rad]
static const double theta_stp = 5.0 * M_PI / 180;       // 回転角の間隔 [rad]
double theta_deg = 0.0;                                 // csvファイルの第1カラム(角度 [°])

// ファイル
std::string filenameUpperBoundSum;                      // 厳密な BER 上界の合計
std::vector<std::string> filenameUpperBound;            // 誤りビットごとの厳密な BER 上界の合計
std::string filenameNearlyUpperBoundSum;                // BER 上界近似の合計
std::vector<std::string> filenameNearlyUpperBound;      // 誤りビット数ごとの BER 上界近似

std::ofstream ofsUpperBoundSum;                         // 厳密な BER 上界の合計
std::vector<std::ofstream> ofsUpperBound;               // 誤りビット数ごとの 厳密な BER 上界
std::ofstream ofsNearlyUpperBoundSum;                   // BER 上界近似の合計
std::vector<std::ofstream> ofsNearlyUpperBound;         // 誤りビット数ごとのBER 上界近似

// BER
Eigen::VectorXd berUpperBoundVec;                       // BER上界
Eigen::VectorXd berNearlyUpperBoundVec;                 // BER上界の近似値

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
            filenameUpperBoundSum = "RotatedQPSK_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsUpperBoundSum.open(filenameUpperBoundSum);
            filenameNearlyUpperBoundSum = "RotatedQPSK_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

            for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                filenameUpperBound[i] = "RotatedQPSK_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsUpperBound[i].open(filenameUpperBound[i]);
                filenameNearlyUpperBound[i] = "RotatedQPSK_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsNearlyUpperBound[i].open(filenameNearlyUpperBound[i]);
            }

            sim.set_QPSKNoiseSD(EbN0dB);

            for(double theta = theta_min; theta <= theta_max; theta += theta_stp) {
                sim.setRotationSymbol(theta);       // 回転

                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                std::cout << "--------------------------------------------" << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                // ファイル出力
                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsUpperBound[i] << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                ofsUpperBoundSum << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsNearlyUpperBound[i] << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                ofsNearlyUpperBoundSum << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                theta_deg += 0.5;     // csvファイルの第1カラムをインクリメント
            }
        break;
        case 4:
            // ファイル作成
            filenameUpperBoundSum = "Rotated16QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsUpperBoundSum.open(filenameUpperBoundSum);
            filenameNearlyUpperBoundSum = "Rotated16QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

            for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                filenameUpperBound[i] = "Rotated16QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsUpperBound[i].open(filenameUpperBound[i]);
                filenameNearlyUpperBound[i] = "Rotated16QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsNearlyUpperBound[i].open(filenameNearlyUpperBound[i]);
            }

            sim.set_16QAMNoiseSD(EbN0dB);

            for(double theta = theta_min; theta <= theta_max; theta += theta_stp) { 
                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                std::cout << "--------------------------------------------" << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                // ファイル出力
                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsUpperBound[i] << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                ofsUpperBoundSum << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsNearlyUpperBound[i] << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                ofsNearlyUpperBoundSum << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                theta_deg += 5.0;     // csvファイルの第1カラムをインクリメント
            }
        break;
    }
    
    ofsUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsUpperBound[i].close();
    }
    ofsNearlyUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsNearlyUpperBound[i].close();
    }

    return 0;
}