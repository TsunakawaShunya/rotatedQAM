// L次ダイバーシチ（ガウスの超幾何関数を使用）
#include "simulator.h"

// SNR
double EbN0dB;

// ファイル
std::string filename;
std::ofstream ofs;

double ber;

// ダイバーシチ次数
int Lmax = 30;

int main() {
    Simulator sim;
    // SNR設定
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "EbN0dB [dB]?" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> EbN0dB;

    sim.setSymbol();        // 従来QAMでのシンボル設計

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            filename = "4QAM_LDiversity_theory_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofs.open(filename);

            for(int L = 1; L <= Lmax; L++) {

                // 標準出力
                ber = sim.get_4QAMTheory_Ldiversity(EbN0dB, L);
                std::cout << "Theory" << L << "," << ber << std::endl;
                // ファイル出力
                ofs << L << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 4:
            filename = "16QAM_LDiversity_hyp_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofs.open(filename);

            for(int L = 1; L <= Lmax; L++) {

                // 標準出力
                ber = sim.get_16QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp" << L << "," << ber << std::endl;

                // ファイル出力
                ofs << L << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 6:
            filename = "64QAM_LDiversity_hyp_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofs.open(filename);

            for(int L = 1; L <= Lmax; L++) {

                // 標準出力
                ber = sim.get_64QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp" << L << "," << ber << std::endl;

                // ファイル出力
                ofs << L << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 8:
            filename = "256QAM_LDiversity_hyp_" + std::to_string((int)EbN0dB) + "dB.csv";
            ofs.open(filename);

            for(int L = 1; L <= Lmax; L++) {

                // 標準出力
                ber = sim.get_256QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp" << L << "," << ber << std::endl;

                // ファイル出力
                ofs << L << "," << ber << std::endl;
            }
            ofs.close();
        break;
    }
    return 0;
}