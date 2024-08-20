// SNに対するBER
#include "simulator.h"
#include <vector>

// SNR
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 40.0;        // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 5.0;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// 回転角
double theta_deg;
double theta;

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
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "SN - BER" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "theta? [deg]" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> theta_deg;
    theta = theta_deg * M_PI / 180;

    // シンボル設計
    sim.setSymbol();        // 従来QAMでのシンボル設計
    sim.setRotationSymbol(theta);       // 回転


    // 配列の初期化
    filenameUpperBound.resize(sim.NUMBER_OF_BIT);
    ofsUpperBound.resize(sim.NUMBER_OF_BIT);
    berUpperBoundVec.resize(sim.NUMBER_OF_BIT);
    filenameNearlyUpperBound.resize(sim.NUMBER_OF_BIT);
    ofsNearlyUpperBound.resize(sim.NUMBER_OF_BIT);
    berNearlyUpperBoundVec.resize(sim.NUMBER_OF_BIT);

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            // ファイル作成
            //filename = "RotatedQPSK_sim_" + std::to_string(theta_deg) + "deg.csv";
            //ofs.open(filename);
            filenameUpperBoundSum = "RotatedQPSK_upperbound_" + std::to_string(theta_deg) + "deg.csv";
            ofsUpperBoundSum.open(filenameUpperBoundSum);
            filenameNearlyUpperBoundSum = "RotatedQPSK_nearly_upperbound_" + std::to_string(theta_deg) + "deg.csv";
            ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

            for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                filenameUpperBound[i] = "RotatedQPSK_upperbound_" + std::to_string(theta_deg) + "deg_" + std::to_string(i + 1) + ".csv";
                ofsUpperBound[i].open(filenameUpperBound[i]);
                filenameNearlyUpperBound[i] = "RotatedQPSK_nearly_upperbound_" + std::to_string(theta_deg) + "deg_" + std::to_string(i + 1) + ".csv";
                ofsNearlyUpperBound[i].open(filenameNearlyUpperBound[i]);
            }

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                sim.set_QPSKNoiseSD(EbN0dB);

                //ber = sim.getBerSimulation();
                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                std::cout << "--------------------------------------------" << std::endl;
                //std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
                std::cout << "--------------------------------------------" << std::endl;

                // ファイル出力
                //ofs << EbN0dB << "," << ber << std::endl;
                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsUpperBound[i] << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
                }
                ofsUpperBoundSum << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;

                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsNearlyUpperBound[i] << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                ofsNearlyUpperBoundSum << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
            }
        break;
        case 4:
            // ファイル作成
            /*
            filename = "Rotated16QAM_sim_" + std::to_string(theta_deg) + "deg.csv";
            ofs.open(filename);
            */
            filenameUpperBoundSum = "Rotated16QAM_upperbound_" + std::to_string(theta_deg) + "deg.csv";
            ofsUpperBoundSum.open(filenameUpperBoundSum);
            filenameNearlyUpperBoundSum = "Rotated16QAM_nearly_upperbound_" + std::to_string(theta_deg) + "deg.csv";
            ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

            for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                filenameUpperBound[i] = "Rotated16QAM_upperbound_" + std::to_string(theta_deg) + "deg_" + std::to_string(i + 1) + ".csv";
                ofsUpperBound[i].open(filenameUpperBound[i]);
                filenameNearlyUpperBound[i] = "Rotated16QAM_nearly_upperbound_" + std::to_string(theta_deg) + "deg_" + std::to_string(i + 1) + ".csv";
                ofsNearlyUpperBound[i].open(filenameNearlyUpperBound[i]);
            }
            


            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                sim.set_16QAMNoiseSD(EbN0dB);

                //ber = sim.getBerSimulation();
                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                /*
                std::cout << "--------------------------------------------" << std::endl;
                std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;
                */
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
               
                // ファイル出力
                //ofs << EbN0dB << "," << ber << std::endl;

                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsUpperBound[i] << EbN0dB << "," << berUpperBoundVec(i) << std::endl;
                }
                ofsUpperBoundSum << EbN0dB << "," << berUpperBoundVec.sum() << std::endl;

                for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    ofsNearlyUpperBound[i] << EbN0dB << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                ofsNearlyUpperBoundSum << EbN0dB << "," << berNearlyUpperBoundVec.sum() << std::endl;
            }
        break;
    }

    // 終了
    ofs.close();

    /*
    ofsUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsUpperBound[i].close();
    }
    ofsNearlyUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsNearlyUpperBound[i].close();
    }
    */

    return 0;
}
