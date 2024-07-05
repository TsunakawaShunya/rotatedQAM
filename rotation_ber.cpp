// 回転角に対するBER
#include "simulator.h"
#include <vector>

// SNR
double EbN0dB;

// 回転角
static const double theta_min = 0.5 * M_PI / 180;        // 回転角の最小値 [rad]
static const double theta_max = 45 * M_PI / 180 ;        // 回転角の最大値 [rad]
static const double theta_stp = 0.5 * M_PI / 180;        // 回転角の間隔 [rad]
double theta_deg = 0.5;       // csvファイルの第1カラム(角度 [°])

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
            // ファイル作成
            //filename = "RotatedQPSK_sim_" + std::to_string((int)EbN0dB) + "dB.csv";
            //ofs.open(filename);
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

                //ber = sim.getBerSimulation();
                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                std::cout << "--------------------------------------------" << std::endl;
                //std::cout << "simulation : " << theta_deg << "," << ber << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                // ファイル出力
                //ofs << theta_deg << "," << ber << std::endl;
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
            // ファイル作成
            //filename = "RotatedQPSK_sim_" + std::to_string((int)EbN0dB) + "dB.csv";
            //ofs.open(filename);
            filenameUpperBoundSum = "Rotated16QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsUpperBoundSum.open(filenameUpperBoundSum);
            filenameNearlyUpperBoundSum = "Rotated16QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofsNearlyUpperBoundSum.open(filenameNearlyUpperBoundSum);

            sim.set_16QAMNoiseSD(EbN0dB);

            for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                filenameUpperBound[i] = "Rotated16QAM_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsUpperBound[i].open(filenameUpperBound[i]);
                filenameNearlyUpperBound[i] = "Rotated16QAM_nearly_upperbound_" + std::to_string((int)EbN0dB) + "dB_" + std::to_string(i + 1) + ".csv";
                ofsNearlyUpperBound[i].open(filenameNearlyUpperBound[i]);
            }

            for(double theta = theta_min; theta <= theta_max; theta += theta_stp) { 
                sim.setRotationSymbol(theta);       // 回転 

                //ber = sim.getBerSimulation();
                berUpperBoundVec = sim.getBerUpperBoundVec(theta);
                berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(theta);

                // 標準出力
                std::cout << "--------------------------------------------" << std::endl;
                //std::cout << "simulation : " << theta_deg << "," << ber << std::endl;
                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "upperbound(" << i << ") : " << theta_deg << "," << berUpperBoundVec(i) << std::endl;
                }
                std::cout << "upperbound : " << theta_deg << "," << berUpperBoundVec.sum() << std::endl;

                for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
                    std::cout << "nearly_upperbound(" << i << ") : " << theta_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
                }
                std::cout << "nearly_upperbound : " << theta_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

                // ファイル出力
                //ofs << theta_deg << "," << ber << std::endl;
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
    }

    // ------------------------------------最適な回転角をコンピュータサーチによって導出------------------------------------
    // 最適化の目的関数を定義
    std::function<double(double)>func = [&sim](double theta) { return sim.getBerUpperBoundVec(theta).sum(); };

    // 初期推定値と最適化の実行
    double initial_rad = 20 * M_PI / 180;
    double optimal_rad = sim.get_optimizeTheta_rad_byNewton(func, initial_rad);
    double optimal_deg = optimal_rad * 180 / M_PI;

    /*
    // 最適な回転角での結果
    //ber = sim.getBerSimulation();
    berUpperBoundVec = sim.getBerUpperBoundVec(optimal_rad);
    berNearlyUpperBoundVec = sim.get_nearlyBerUpperBoundVec(optimal_rad);

    // 標準出力
    //std::cout << "simulation : " << optimal_deg << "," << ber << std::endl;
    for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        std::cout << "upperbound(" << i << ") : " << optimal_deg << "," << berUpperBoundVec(i) << std::endl;
    }
    std::cout << "upperbound : " << optimal_deg << "," << berUpperBoundVec.sum() << std::endl;
    for(int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        std::cout << "nearly_upperbound(" << i << ") : " << optimal_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
    }
    std::cout << "nearly_upperbound : " << optimal_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;

    // ファイル出力
    //ofs << optimal_deg << "," << ber << std::endl;
    for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        ofsUpperBound[i] << optimal_deg << "," << berUpperBoundVec(i) << std::endl;
    }
    ofsUpperBoundSum << optimal_deg << "," << berUpperBoundVec.sum() << std::endl;
    for (int i = 0; i < sim.NUMBER_OF_BIT; i++) {
        ofsNearlyUpperBound[i] << optimal_deg << "," << berNearlyUpperBoundVec(i) << std::endl;
    }
    ofsNearlyUpperBoundSum << optimal_deg << "," << berNearlyUpperBoundVec.sum() << std::endl;
    */

    // 終了
    ofsUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsUpperBound[i].close();
    }
    ofsNearlyUpperBoundSum.close();
    for (int i = 0; i < sim.NUMBER_OF_BIT; ++i) {
        ofsNearlyUpperBound[i].close();
    }

    // 最適な回転角を標準出力
    std::cout << "Optimal theta (deg): " << optimal_deg << std::endl;

    return 0;
}