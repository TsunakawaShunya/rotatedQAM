// L次ダイバーシチ（シミュレーション）

#include "simulator.h"

// SNR
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 40.1;       // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 5.0;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// ダイバーシチ
int L;

// ファイル
std::string filename;
std::ofstream ofs;

// BER
double ber;                                 // BERシミュレーション値

int main() {   
    Simulator sim;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "L?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> L;

    // シンボル設計
    sim.setSymbol();                        // 従来QAMでのシンボル設計

    // ファイルの初期化
    int M = std::pow(sim.NUMBER_OF_BIT, 2);     // 多値数
    filename = std::to_string(M) + "QAM_" + std::to_string(L) + "diversity_sim.csv";
    ofs.open(filename);

    // Eb/N0[dB]でループ
    for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        sim.setNoiseSD(EbN0dB);

        // シミュレーション
        ber = sim.getBerSimulation_Ldiversity(L);

        // 標準出力
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << "simulation : " << EbN0dB << "," << ber << std::endl;

        // ファイル出力
        ofs << EbN0dB << "," << ber << std::endl;
    }
    ofs.close();

    //sim.checkSymbolPower();

    return 0;
}
