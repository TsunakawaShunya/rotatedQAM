// L次ダイバーシチ（ガウスの超幾何関数を使用）
#include "simulator.h"

// SNR
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 30.1;        // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 0.1;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// ファイル
std::string filename;
std::ofstream ofs;

double ber;

// ダイバーシチ次数
int L;

int main() {
    Simulator sim;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Number of Diversity? (L)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> L;

    sim.setSymbol();        // 従来QAMでのシンボル設計

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            filename = "4QAM_" + std::to_string(L) + "Diversity_hyp.csv";
            ofs.open(filename);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {

                // 標準出力
                ber = sim.get_4QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp(" << L << "-Diversity) : " << EbN0dB << "," << ber << std::endl;

                // ファイル出力
                ofs << EbN0dB << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 4:
            filename = "16QAM_" + std::to_string(L) + "Diversity_hyp.csv";
            ofs.open(filename);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {

                // 標準出力
                ber = sim.get_16QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp(" << L << "-Diversity) : " << EbN0dB << "," << ber << std::endl;

                // ファイル出力
                ofs << EbN0dB << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 6:
            filename = "64QAM_" + std::to_string(L) + "Diversity_hyp.csv";
            ofs.open(filename);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {

                // 標準出力
                ber = sim.get_64QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp(" << L << "-Diversity) : " << EbN0dB << "," << ber << std::endl;

                // ファイル出力
                ofs << EbN0dB << "," << ber << std::endl;
            }
            ofs.close();
        break;
        case 8:
            filename = "256QAM_" + std::to_string(L) + "Diversity_hyp.csv";
            ofs.open(filename);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {

                // 標準出力
                ber = sim.get_256QAMTheory_Ldiversity_hyp(EbN0dB, L);
                std::cout << "Theory hyp(" << L << "-Diversity) : " << EbN0dB << "," << ber << std::endl;

                // ファイル出力
                ofs << EbN0dB << "," << ber << std::endl;
            }
            ofs.close();
        break;
    }
    return 0;
}