// SNに対するBER
/*  
    最適な回転角
    QPSK: 27.3678°
    16QAM: 21.0148°
    64QAM: 20.9291°
    256QAM: 21.0098°
*/
#include "simulator.h"
#include <vector>

// SNR
static const double EbN0dBmin = 0.0;   // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 40.0;  // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 0.1;   // Eb/N0 の間隔 [dB]

// 回転角
double theta_deg;
double theta;

// ファイル名およびストリーム
std::ofstream ofsUpperBoundSum;
std::vector<std::ofstream> ofsUpperBound;
std::ofstream ofsNearlyUpperBoundSum;
std::vector<std::ofstream> ofsNearlyUpperBound;

// BER
Eigen::VectorXd berUpperBoundVec;
Eigen::VectorXd berNearlyUpperBoundVec;

// ファイル作成の関数
void initializeFiles(Simulator& sim, std::ofstream& ofsUpperBoundSum,
                std::ofstream& ofsNearlyUpperBoundSum, 
                std::vector<std::ofstream>& ofsUpperBound, 
                std::vector<std::ofstream>& ofsNearlyUpperBound) {
    // 多値数
    int M = static_cast<int>(std::pow(2, sim.NUMBER_OF_BIT));

    // ファイル名の作成
    std::string filenameUpperBoundSum = "Rotated" + std::to_string(M) + "QAM_upperbound_" + std::to_string(theta_deg) + "deg.csv";
    ofsUpperBoundSum.open(filenameUpperBoundSum);

    std::string filenameNearlyUpperBoundSum = "Rotated" + std::to_string(M) + "QAM_nearly_upperbound_" + std::to_string(theta_deg) + "deg.csv";
    ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

    // 各ビットごとのファイル名を作成
    for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        std::string filenameUpperBound = "Rotated" + std::to_string(M) + "QAM_upperbound_" + std::to_string(theta_deg) + "deg.csv" + std::to_string(i + 1) + ".csv";
        ofsUpperBound[i].open(filenameUpperBound);

        std::string filenameNearlyUpperBound = "Rotated" + std::to_string(M) + "QAM_nearly_upperbound_" + std::to_string(theta_deg) + "deg.csv" + std::to_string(i + 1) + ".csv";
        ofsNearlyUpperBound[i].open(filenameNearlyUpperBound);
    }
}

// Eb/N0 でループし、結果を出力およびファイル書き込み
void process(Simulator& sim, std::ofstream& ofsUpperBoundSum, std::ofstream& ofsNearlyUpperBoundSum,
                  std::vector<std::ofstream>& ofsUpperBound, std::vector<std::ofstream>& ofsNearlyUpperBound) {
    for (double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        sim.setNoiseSD(EbN0dB);

        // BER の取得
        berUpperBoundVec = sim.getBerUpperBoundVec(theta);
        berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

        // 標準出力
        std::cout << "--------------------------------------------" << std::endl;
        for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
            std::cout << "upperbound(" << i << ") : " << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
            std::cout << "nearly_upperbound(" << i << ") : " << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
        }
        std::cout << "upperbound : " << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;
        std::cout << "nearly_upperbound : " << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        // ファイル書き込み
        for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
            ofsUpperBound[i] << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
            ofsNearlyUpperBound[i] << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
        }
        ofsUpperBoundSum << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;
        ofsNearlyUpperBoundSum << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
    }
}

// main部
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
    sim.setSymbol();
    sim.setRotationSymbol(theta);

    // 配列の初期化
    std::vector<std::ofstream> ofsUpperBound(sim.NUMBER_OF_BIT);
    std::vector<std::ofstream> ofsNearlyUpperBound(sim.NUMBER_OF_BIT);
    berUpperBoundVec.resize(sim.NUMBER_OF_BIT);
    berNearlyUpperBoundVec.resize(sim.NUMBER_OF_BIT);

    std::ofstream ofsUpperBoundSum;
    std::ofstream ofsNearlyUpperBoundSum;

    // ファイル生成と処理の呼び出し
    initializeFiles(sim, ofsUpperBoundSum, ofsNearlyUpperBoundSum, ofsUpperBound, ofsNearlyUpperBound);
    process(sim, ofsUpperBoundSum, ofsNearlyUpperBoundSum, ofsUpperBound, ofsNearlyUpperBound);

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