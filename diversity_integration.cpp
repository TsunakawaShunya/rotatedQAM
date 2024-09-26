// L次のダイバーシチ（数値積分）
#include "simulator.h"

// SNR
static const double EbN0dBmin = 0.0;        // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 40.0;        // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 5.0;        // Eb/N0 の間隔 [dB]
double EbN0dB;

// ファイル
std::string filenameTheory_1;
std::string filenameTheory_2;
std::ofstream ofsTheory_1;
std::ofstream ofsTheory_2;

double berTheory_1;     // 1次ダイバーシチ
double berTheory_2;        // 2次ダイバーシチ

int main() {
    Simulator sim;

    sim.setSymbol();        // 従来QAMでのシンボル設計

    switch(sim.NUMBER_OF_BIT) {
        case 2:
            filenameTheory_1 = "QPSK_theory_1Diversity_int.csv";
            filenameTheory_2 = "QPSK_theory_2Diversity_int.csv";
            ofsTheory_1.open(filenameTheory_1);
            ofsTheory_2.open(filenameTheory_2);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                sim.set_QPSKNoiseSD(EbN0dB);

                // 標準出力
                berTheory_1 = sim.get_4QAMTheory_Ldiversity_int(EbN0dB, 1);
                std::cout << "Theory int(1-Diversity) : " << EbN0dB << "," << berTheory_1 << std::endl;
                berTheory_2 = sim.get_4QAMTheory_Ldiversity_int(EbN0dB, 2);
                std::cout << "Theory int(2-Diversity) : " << EbN0dB << "," << berTheory_2 << std::endl;

                // ファイル出力
                ofsTheory_1 << EbN0dB << "," << berTheory_1 << std::endl;
                ofsTheory_2 << EbN0dB << "," << berTheory_2 << std::endl;
            }
            ofsTheory_1.close();
            ofsTheory_2.close();
        break;
        case 4:
            filenameTheory_1 = "16QAM_theory_1Diversity_int.csv";
            filenameTheory_2 = "16QAM_theory_2Diversity_int.csv";
            ofsTheory_1.open(filenameTheory_1);
            ofsTheory_2.open(filenameTheory_2);

            for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
                sim.set_QPSKNoiseSD(EbN0dB);
                // 標準出力
                berTheory_1 = sim.get_16QAMTheory_Ldiversity_int(EbN0dB, 1);
                std::cout << "Theory int(1-Diversity) : " << EbN0dB << "," << berTheory_1 << std::endl;
                berTheory_2 = sim.get_16QAMTheory_Ldiversity_int(EbN0dB, 2);
                std::cout << "Theory int(2-Diversity) : " << EbN0dB << "," << berTheory_2 << std::endl;

                // ファイル出力
                ofsTheory_1 << EbN0dB << "," << berTheory_1 << std::endl;
                ofsTheory_2 << EbN0dB << "," << berTheory_2 << std::endl;
            }
            ofsTheory_1.close();
            ofsTheory_2.close();
        break;
    }
    return 0;
}