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

    // ファイルの初期化
    int M = std::pow(sim.NUMBER_OF_BIT, 2);     // 多値数
    filename = "Rotated" + std::to_string(M) + "QAM_sim_" + std::to_string((int)EbN0dB) + "dB.csv";
    ofs.open(filename);

    for(double theta = theta_min; theta <= theta_max; theta += theta_stp) {
        // 回転
        sim.setRotationSymbol(theta);

        // シミュレーション
        ber = sim.getBerSimulation();

        // 標準出力
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "simulation : " << theta_deg << "," << ber << std::endl;

        // ファイル出力
        ofs << theta_deg << "," << ber << std::endl;

        theta_deg += 5.0;     // csvファイルの第1カラムをインクリメント
    }

    /*        
    // 最適な回転角の結果を追加（4QAM）
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


    /*        
    // 最適な回転角の結果を追加（16QAM）
    // 21.0148
    theta_deg = 21.0148;
    //std::cout << theta_opt * 180 / M_PI << std::endl;
    sim.setRotationSymbol(theta_opt);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;

    theta_deg = 31.7;
    sim.setRotationSymbol(theta_deg * M_PI / 180);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;

    theta_deg = 35.8;
    sim.setRotationSymbol(theta_deg * M_PI / 180);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;
    ofs.close();
    */

       /*        
    // 最適な回転角の結果を追加（16QAM）
    // 21.0148
    theta_deg = 21.0148;
    //std::cout << theta_opt * 180 / M_PI << std::endl;
    sim.setRotationSymbol(theta_opt);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;

    theta_deg = 31.7;
    sim.setRotationSymbol(theta_deg * M_PI / 180);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;

    theta_deg = 35.8;
    sim.setRotationSymbol(theta_deg * M_PI / 180);
    ber = sim.getBerSimulation();
    std::cout << theta_deg << "," << ber << std::endl;
    ofs << theta_deg << "," << ber << std::endl;
    ofs.close();
    */

    return 0;
}