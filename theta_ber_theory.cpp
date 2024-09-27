#include "simulator.h"
#include <vector>
#include <cmath>    // std::pow

// SNR
double EbN0dB;

// 回転角
static const double theta_min = 0.0 * M_PI / 180;       // 回転角の最小値 [rad]
static const double theta_max = 45.0 * M_PI / 180;      // 回転角の最大値 [rad]
static const double theta_stp = 0.5 * M_PI / 180;       // 回転角の間隔 [rad]
double theta_deg = 0.0;                                 // csvファイルの第1カラム(角度 [°])

// BER
Eigen::VectorXd berUpperBoundVec;                       // BER上界
Eigen::VectorXd berNearlyUpperBoundVec;                 // BER上界の近似値

// ファイルを初期化して開く関数
void initializeFiles(Simulator& sim, std::ofstream& ofsUpperBoundSum,
                     std::ofstream& ofsNearlyUpperBoundSum, std::vector<std::ofstream>& ofsUpperBound, std::vector<std::ofstream>& ofsNearlyUpperBound) {

    // 2^NUMBER_OF_BITの計算 (QAMの階数)
    int modulationOrder = static_cast<int>(std::pow(2, sim.NUMBER_OF_BIT));

    // ファイル名の作成
    std::string filenameUpperBoundSum = "Rotated" + std::to_string(modulationOrder) + "QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
    ofsUpperBoundSum.open(filenameUpperBoundSum);

    std::string filenameNearlyUpperBoundSum = "Rotated" + std::to_string(modulationOrder) + "QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
    ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

    // 各ビットごとのファイル名を作成
    for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        std::string filenameUpperBound = "Rotated" + std::to_string(modulationOrder) + "QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
        ofsUpperBound[i].open(filenameUpperBound);

        std::string filenameNearlyUpperBound = "Rotated" + std::to_string(modulationOrder) + "QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
        ofsNearlyUpperBound[i].open(filenameNearlyUpperBound);
    }
}

// ファイル出力と標準出力を行う共通関数
void processTheta(Simulator& sim, std::ofstream& ofsUpperBoundSum, std::ofstream& ofsNearlyUpperBoundSum,
                  std::vector<std::ofstream>& ofsUpperBound, std::vector<std::ofstream>& ofsNearlyUpperBound) {

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
}

int main() {
    Simulator sim;

    // SNR設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "EbN0dB [dB]?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> EbN0dB;

    sim.setSymbol();        // 従来QAMシンボル設計

    // 配列の初期化
    std::vector<std::ofstream> ofsUpperBound(sim.NUMBER_OF_BIT);
    std::vector<std::ofstream> ofsNearlyUpperBound(sim.NUMBER_OF_BIT);
    berUpperBoundVec.resize(sim.NUMBER_OF_BIT);
    berNearlyUpperBoundVec.resize(sim.NUMBER_OF_BIT);

    std::ofstream ofsUpperBoundSum;
    std::ofstream ofsNearlyUpperBoundSum;

    // ファイル生成と処理の呼び出し
    initializeFiles(sim, ofsUpperBoundSum, ofsNearlyUpperBoundSum, ofsUpperBound, ofsNearlyUpperBound);
    sim.set_NoiseSD(EbN0dB);
    processTheta(sim, ofsUpperBoundSum, ofsNearlyUpperBoundSum, ofsUpperBound, ofsNearlyUpperBound);

    // ファイルを閉じる
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