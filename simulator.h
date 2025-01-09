/*
 * File:   Simulator.h
 * Author: Tsunakawa
 *
 * Created on 2023/04/29, 15:49
*/

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES
#include <C:eigen-3.4.0/Eigen/Dense>
#include <C:eigen-3.4.0/Eigen/Eigen>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <bitset>
#include <cmath>
#include <functional>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <initializer_list>
#include <stdexcept>
#include "random_collection.h"

class Simulator {
    public:
    int NUMBER_OF_BIT;      // ビット数(QAM方式)

    // コンストラクタ
    Simulator() {
        // 次元を設定
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Number of Bit? (QPSK:2, 16QAM:4, 64QAM:6, 256QAM:8, 1024QAM:10)" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cin >> NUMBER_OF_BIT;

        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "TRIAL?" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cin >> N_Tri;

        numberOfSymbols_ = pow(2, NUMBER_OF_BIT);

        // リサイズ
        num_.resize(numberOfSymbols_);
        grayNum_.resize(numberOfSymbols_);
        symbol_.resize(numberOfSymbols_);
        rotatedSymbol_.resize(numberOfSymbols_);
        phi_.resize(numberOfSymbols_);
        distance_.resize(numberOfSymbols_);
        // rotationMatrix.resize(2, 2);

        // データ（0 ~ 2^M - 1）をセット
        setNum();
        
        // 乱数の設定 
        unitIntUniformRand_.init(0, numberOfSymbols_ - 1, seed);
        unitCNormalRand_.init(0, M_SQRT1_2, seed);
        unitNormalRand_.init(0, 1, seed);
    }

    // デストラクタ
    virtual ~Simulator() {
    }

    // シンボル設計
    void setSymbol() {
        int M = numberOfSymbols_;               // シンボル数 (M = 2^NUMBER_OF_BIT)
        int sqrtM = sqrt(M);                    // 実部/虚部のレベル数 (例: 16QAMならsqrtM=4)
        double P = 1.0 / (2.0 * (M - 1) / 3.0);

        // シンボル設計
        int i = 0;
        for (int v1 = 0; v1 < sqrtM; v1++) {
            for (int v2 = 0; v2 < sqrtM; v2++) {
                symbol_(i).real((2 * v1 - (sqrtM - 1)) * sqrt(P));  // 実部
                if (v1 % 2 == 0) {  // v1が偶数のときは通常の配置
                    symbol_(i).imag((2 * v2 - (sqrtM - 1)) * sqrt(P));  // 虚部
                } else {  // v1が奇数のとき、虚部の値を逆順にする
                    symbol_(i).imag(((sqrtM - 1) - 2 * v2) * sqrt(P));
                }
                i++;
            }
        }
    }

    // 一般化する前のシンボル設計メソッド
    void setSymbol_notGene() {
        switch(NUMBER_OF_BIT) {
            // 4QAM(QPSK)
            case 2:
                P = 1.0 / 2.0;
                for(int i = 0; i < numberOfSymbols_; i++) {
                    // 1ビット目
                    if((grayNum_[i] & 1) == 0) {
                        symbol_(i).real(sqrt(P));
                    } else {
                        symbol_(i).real(-sqrt(P));
                    }
                    // 2ビット目
                    if(((grayNum_[i] >> 1) & 1) == 0) {
                        symbol_(i).imag(sqrt(P));
                    } else {
                        symbol_(i).imag(-sqrt(P));
                    }
                }
                break;

            // 16QAM
            case 4:
                P = 1.0 / 10.0;
                for(int i = 0; i < numberOfSymbols_; i++) {
                    // 実部の設計(3, 4ビット目で決まる)
                    switch((grayNum_[i] >> 2) & 0b11) {
                        case 0b00: symbol_(i).real(-3 * sqrt(P)); break;
                        case 0b01: symbol_(i).real(-sqrt(P)); break;
                        case 0b11: symbol_(i).real(sqrt(P)); break;
                        case 0b10: symbol_(i).real(3 * sqrt(P)); break;
                    }
                    // 虚部の設計(1, 2ビット目で決まる)
                    switch(grayNum_[i] & 0b11) {
                        case 0b00: symbol_(i).imag(-3 * sqrt(P)); break;
                        case 0b01: symbol_(i).imag(-sqrt(P)); break;
                        case 0b11: symbol_(i).imag(sqrt(P)); break;
                        case 0b10: symbol_(i).imag(3 * sqrt(P)); break;
                    }
                }
                break;

            // 64QAM
            case 6:
                P = 1.0 / 42.0;
                for(int i = 0; i < numberOfSymbols_; i++) {
                    // 実部の設計(5, 6ビット目で決まる)
                    switch((grayNum_[i] >> 4) & 0b11) {
                        case 0b00: symbol_(i).real(-7 * sqrt(P)); break;
                        case 0b01: symbol_(i).real(-5 * sqrt(P)); break;
                        case 0b11: symbol_(i).real(-3 * sqrt(P)); break;
                        case 0b10: symbol_(i).real(-sqrt(P)); break;
                        case 0b100: symbol_(i).real(sqrt(P)); break;
                        case 0b101: symbol_(i).real(3 * sqrt(P)); break;
                        case 0b111: symbol_(i).real(5 * sqrt(P)); break;
                        case 0b110: symbol_(i).real(7 * sqrt(P)); break;
                    }
                    // 虚部の設計(1, 2, 3ビット目で決まる)
                    switch(grayNum_[i] & 0b11) {
                        case 0b00: symbol_(i).imag(-7 * sqrt(P)); break;
                        case 0b01: symbol_(i).imag(-5 * sqrt(P)); break;
                        case 0b11: symbol_(i).imag(-3 * sqrt(P)); break;
                        case 0b10: symbol_(i).imag(-sqrt(P)); break;
                        case 0b100: symbol_(i).imag(sqrt(P)); break;
                        case 0b101: symbol_(i).imag(3 * sqrt(P)); break;
                        case 0b111: symbol_(i).imag(5 * sqrt(P)); break;
                        case 0b110: symbol_(i).imag(7 * sqrt(P)); break;
                    }
                }
                break;

            // 256QAM
            case 8:
                P = 1.0 / 170.0;
                for(int i = 0; i < numberOfSymbols_; i++) {
                    // 実部の設計(7, 8ビット目で決まる)
                    switch((grayNum_[i] >> 6) & 0b11) {
                        case 0b00: symbol_(i).real(-15 * sqrt(P)); break;
                        case 0b01: symbol_(i).real(-13 * sqrt(P)); break;
                        case 0b11: symbol_(i).real(-11 * sqrt(P)); break;
                        case 0b10: symbol_(i).real(-9 * sqrt(P)); break;
                        case 0b100: symbol_(i).real(-7 * sqrt(P)); break;
                        case 0b101: symbol_(i).real(-5 * sqrt(P)); break;
                        case 0b111: symbol_(i).real(-3 * sqrt(P)); break;
                        case 0b110: symbol_(i).real(-sqrt(P)); break;
                        case 0b1000: symbol_(i).real(sqrt(P)); break;
                        case 0b1001: symbol_(i).real(3 * sqrt(P)); break;
                        case 0b1011: symbol_(i).real(5 * sqrt(P)); break;
                        case 0b1010: symbol_(i).real(7 * sqrt(P)); break;
                        case 0b1100: symbol_(i).real(9 * sqrt(P)); break;
                        case 0b1101: symbol_(i).real(11 * sqrt(P)); break;
                        case 0b1111: symbol_(i).real(13 * sqrt(P)); break;
                        case 0b1110: symbol_(i).real(15 * sqrt(P)); break;
                    }
                    // 虚部の設計(1, 2, 3ビット目で決まる)
                    switch(grayNum_[i] & 0b11) {
                        case 0b00: symbol_(i).imag(-15 * sqrt(P)); break;
                        case 0b01: symbol_(i).imag(-13 * sqrt(P)); break;
                        case 0b11: symbol_(i).imag(-11 * sqrt(P)); break;
                        case 0b10: symbol_(i).imag(-9 * sqrt(P)); break;
                        case 0b100: symbol_(i).imag(-7 * sqrt(P)); break;
                        case 0b101: symbol_(i).imag(-5 * sqrt(P)); break;
                        case 0b111: symbol_(i).imag(-3 * sqrt(P)); break;
                        case 0b110: symbol_(i).imag(-sqrt(P)); break;
                        case 0b1000: symbol_(i).imag(sqrt(P)); break;
                        case 0b1001: symbol_(i).imag(3 * sqrt(P)); break;
                        case 0b1011: symbol_(i).imag(5 * sqrt(P)); break;
                        case 0b1010: symbol_(i).imag(7 * sqrt(P)); break;
                        case 0b1100: symbol_(i).imag(9 * sqrt(P)); break;
                        case 0b1101: symbol_(i).imag(11 * sqrt(P)); break;
                        case 0b1111: symbol_(i).imag(13 * sqrt(P)); break;
                        case 0b1110: symbol_(i).imag(15 * sqrt(P)); break;
                    }
                }
                break;
        }
    }

    // シンボルの平均電力チェック
    void checkSymbolPower() {
        double avgPower = 0.0;
        for (int i = 0; i < numberOfSymbols_; i++) {
            avgPower += std::norm(symbol_(i));  // シンボルの電力 = 実部^2 + 虚部^2
        }
        avgPower /= numberOfSymbols_;
        std::cout << "Average Symbol Power: " << avgPower << std::endl;
    }

    // シンボル回転
    void setRotationSymbol(double theta) {
        // 回転行列の設定
        rotationMatrix_(0, 0) = cos(theta);
        rotationMatrix_(0, 1) = -sin(theta);
        rotationMatrix_(1, 0) = sin(theta);
        rotationMatrix_(1, 1) = cos(theta);

        for(int i = 0; i < numberOfSymbols_; i++) {
            rotatedSymbol_(i).real(cos(theta) * symbol_(i).real() - sin(theta) * symbol_(i).imag());
            rotatedSymbol_(i).imag(sin(theta) * symbol_(i).real() + cos(theta) * symbol_(i).imag());
        }
    }

    // 加法性雑音の標準偏差を設定
    void setNoiseSD(double EbN0dB) {
        noiseSD_ = sqrt(pow(10.0, -0.1 * EbN0dB) / (double)NUMBER_OF_BIT);
    }

    // 結果
    // シミュレーション
    double getBerSimulation() {
        int berCnt = 0;     // ビット誤りの総数

        for(int tri = 0; tri <= N_Tri; tri++) {
            berCnt += getBitErrorCount();
        }

        return (double)berCnt / (double)N_Tri / (double)NUMBER_OF_BIT / 2.0;
    }

    // L次ダイバーシチシミュレーション
    double getBerSimulation_Ldiversity(int L) {
        int txData;                 // 送信データ
        int rxData;                 // 受信データ
        Eigen::VectorXcd x(L);      // 送信シンボルベクトル
        Eigen::VectorXcd y(L);      // 受信シンボルベクトル
        Eigen::VectorXcd h(L);      // 伝送路応答行列
        Eigen::VectorXcd n(L);      // 雑音
        Eigen::VectorXd obj(numberOfSymbols_);  // 目的関数ベクトル
        int berCnt = 0;                         // エラービット数

        for (int trial = 0; trial < N_Tri; trial++) {
            // 送信データ生成
            txData = unitIntUniformRand_();
            for (int l = 0; l < L; l++) {
                x(l) = symbol_(txData);
                h(l) = unitCNormalRand_();      // 伝送路
                n(l) = unitCNormalRand_();      // 雑音

                // 受信
                y(l) = x(l) + noiseSD_ * sqrt(L) * n(l) / h(l);
            }

            // 最尤推定
            for (int i = 0; i < numberOfSymbols_; i++) {
                obj(i) = 0.0;
                for(int l = 0; l < L; l++) {
                    obj(i) += std::norm(h(l)) * std::norm((y(l) - symbol_(i)));
                }
            }
            Eigen::VectorXd::Index col;
            obj.minCoeff(&col);
            rxData = col;

            // ビットエラー数をカウント
            berCnt += hammingDistance(grayNum_[rxData], grayNum_[txData]);
        }

        // ビット誤り率を計算
        return (double)berCnt / (double)N_Tri / (double)NUMBER_OF_BIT;
    }

    // 理論値
    // 誤りビット数ごとの平均BERの上界
    Eigen::VectorXd getBerUpperBoundVec(double theta_rad) {
        Eigen::VectorXd ber;
        double per;
        int hamming;        // ハミング距離

        set_alphaOmegaKappaGammaZeta(theta_rad);
        ber.resize(NUMBER_OF_BIT);
        ber.setZero();
        
        for(int i = 0; i < numberOfSymbols_; i++) {
            for(int i_hat = 0; i_hat < numberOfSymbols_; i_hat++) {
                hamming = hammingDistance(grayNum_[i], grayNum_[i_hat]);

                if(hamming > 0 && hamming <= NUMBER_OF_BIT) {
                    per = (zeta0_(i, i_hat) * (1 - sqrt(gamma0_(i, i_hat) / (2 + gamma0_(i, i_hat)))) + zeta1_(i, i_hat) * (1 - sqrt(gamma1_(i, i_hat) / (2 + gamma1_(i, i_hat))))) / 2;
                    ber(hamming - 1) += kappa_(i, i_hat) * per;
                }
            }
        }
        return ber / numberOfSymbols_;
    }

    // 平均BERの上界の近似
    Eigen::VectorXd get_nearlyBerUpperBoundVec(double theta) {
        Eigen::VectorXd lambda;
        int hamming;        // iとi_hatのハミング距離

        set_alphaOmegaKappaGammaZeta(theta);
        lambda.resize(NUMBER_OF_BIT);
        lambda.setZero();

        for(int i = 0; i < numberOfSymbols_; i++){
            for(int i_hat = 0; i_hat < numberOfSymbols_; i_hat++) {
                hamming = hammingDistance(grayNum_[i], grayNum_[i_hat]);

                //std::cout << "hamming(" << i << ", " << i_hat << ") : " << hamming << std::endl;
                // 誤りかつ誤りビット数が総ビット数を超えていないとき
                if(hamming > 0 && hamming <= NUMBER_OF_BIT) {
                    lambda(hamming - 1) += kappa_(i, i_hat) / omega0_(i, i_hat) / omega1_(i, i_hat);
                }
            }
        }

        return lambda * 3 * pow(noiseSD_, 4) / numberOfSymbols_;
    }

    // ニュートン法
    // func: 最適化対象の関数, initial_guess: 最適化の初期推定値, max_iter: 最大反復回数（デフォルトは100）, tol: 収束判定のための勾配の閾値（デフォルトは1e-9）
    double get_optimizeTheta_deg_byNewton(std::function<double(double)> func, double initial_deg, int max_iter = 100, double tol = 1e-12) {
        double x_rad = initial_deg * M_PI / 180.0;
        for (int i = 0; i < max_iter; ++i) {
            // 更新
            x_rad = x_rad - numerical_derivative(func, x_rad) / numerical_second_derivative(func, x_rad);
            
            // 収束判定
            if (std::abs(numerical_derivative(func, x_rad)) < tol) {
                std::cout << "------------" << i << " times search------------" << std::endl;
                break;
            }
        }
        return x_rad * 180.0 / M_PI;
    }

    // 解を求める（気合）
    // resolution：分解能
    double get_optimizeTheta_deg(std::function<double(double)> func, double start_deg, double end_deg, double resolution_deg = 1e-04) {
        //std::cout << "initial params" << std::endl;
        //std::cout << start_deg << ", " << end_deg << ", " << resolution_deg << std::endl;

        // rad に変換
        double start_rad = start_deg * M_PI / 180.0;
        double end_rad = end_deg * M_PI / 180.0;
        double resolution_rad = resolution_deg * M_PI / 180.0;

        // x, y の定義
        double y;               // y = func(x)
        double x_rad_min;           // 最小解
        double y_min = 100.0;     // BER（初期値は適当な値）

        for(double x = start_rad; x <= end_rad; x += resolution_rad) {
            y = func(x);
            if(y_min > y) {
                x_rad_min = x;
                y_min = y;
            }
            // std::cout << x * 180.0 / M_PI << ", " << y << std::endl;
        }
        return x_rad_min * 180.0 / M_PI;
    }

    // 4QAM：1次ダイバーシチ[ディジタルコミュニケーションp.898:式(14-3-7)]
    double get_4QAMTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);
        double ber = (1 - sqrt(gamma_b / (1.0 + gamma_b))) / 2;
        return ber;
    }

    // 4QAM：2次のダイバーシチ[ディジタルコミュニケーションp.906:式(14-4-15)]
    double get_4QAMTheory_2diversity(double EbN0dB) {
        double EbN0 = 0.5 * pow(10.0, 0.1 * EbN0dB);
        double mu = sqrt(EbN0 / (1.0 + EbN0));
        double ber = pow((1.0 - mu) / 2.0, 2) * (2.0 + mu);
        return ber;
    }

    // 4QAM：L次のダイバーシチ[ディジタルコミュニケーションp.906:式(14-4-15)]
    double get_4QAMTheory_Ldiversity(double EbN0dB, int L) {
        double gamma_c = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double mu = sqrt(gamma_c / (1.0 + gamma_c));
        double ber = 0;
        for(int l = 0; l < L; l++) {
            ber += combination(L - 1 + l, l) * pow((1.0 + mu) / 2.0, l);
        }
        return ber * pow((1.0 - mu) / 2.0, L);
    }

    // 4QAM：L次ダイバーシチ（積分）
    double get_4QAMTheory_Ldiversity_int(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_b;                         // 瞬時EbN0
        double gamma_b_min = 0.0;
        double gamma_b_max = 100.0;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0
        double ber = 0.0;
        int N_div = 100000;                     // 数値積分の分割数
        double integration = 0.0;               // 積分結果
        double P_AWGN;

        // -------- 積分（台形則）--------
        double dx = (gamma_b_max - gamma_b_min) / (double)N_div;        // 微小区間
        for (int i = 0; i <= N_div; ++i) {
            gamma_b = gamma_b_min + i * dx;
            P_AWGN = Q(sqrt(2 * gamma_b));
            if (i == 0 || i == N_div) {
                integration += P_AWGN
                        * (pow(gamma_b, L - 1) * exp(- gamma_b / gamma_c));
            } else {
                integration += 2 * P_AWGN
                        * (pow(gamma_b, L - 1) * exp(- gamma_b / gamma_c));
            }
        }
        integration = integration * dx / 2.0;
        // ------------------------------

        ber = integration / ((double)factorial(L - 1) * pow(gamma_c, L));
        return ber;
    }

    // 4QAM：L次ダイバーシチ（超幾何関数）
    double get_4QAMTheory_Ldiversity_hyp(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0

        return 1 / (2 * factorial(L - 1) * pow(gamma_c, L)) 
                * boost::math::tgamma(L + 0.5) / (sqrt(M_PI) * L * pow((1 / gamma_c + 1), L))
                * hypergeometric_2F1(0.5, L, L + 1, 1 / (1 + gamma_c));
    }

    // 16QAM理論値
    // 16QAM：1次のダイバーシチ
    double get_16QAMTheory(double EbN0dB) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB);

        return 3.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 / 2.0 / EbN0))
                + 1.0 / 4.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / EbN0 / 9.0))
                - 1.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / EbN0 / 25.0));
    }

    // 16QAM：L次ダイバーシチ（積分）
    double get_16QAMTheory_Ldiversity_int(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_b;                         // 瞬時EbN0
        double gamma_b_min = 0.0;
        double gamma_b_max = 500.0;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0
        double ber = 0.0;
        int N_div = 1000000;                     // 数値積分の分割数
        double integration = 0.0;               // 積分結果
        double P_AWGN;

        // -------- 積分（台形則）--------
        double dx = (gamma_b_max - gamma_b_min) / (double)N_div;        // 微小区間
        for (int i = 0; i <= N_div; ++i) {
            gamma_b = gamma_b_min + i * dx;
            P_AWGN = (3.0 * Q(sqrt(4.0 / 5.0 * gamma_b)) 
                    + 2.0 * Q(3.0 * sqrt(4.0 / 5.0 * gamma_b)) 
                    - Q(5.0 * sqrt(4.0 / 5.0 * gamma_b))) / 4.0;
            if (i == 0 || i == N_div) {
                integration += P_AWGN
                        * (pow(gamma_b, L - 1) * exp(- gamma_b / gamma_c));
            } else {
                integration += 2 * P_AWGN
                        * (pow(gamma_b, L - 1) * exp(- gamma_b / gamma_c));
            }
        }
        integration = integration * dx / 2.0;
        // ------------------------------

        ber = integration / ((double)factorial(L - 1) * pow(gamma_c, L));
        return ber;
    }

    // 16QAM：L次ダイバーシチ（超幾何関数）
    double get_16QAMTheory_Ldiversity_hyp(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0

        return 1 / (4 * factorial(L - 1) * pow(gamma_c, L)) 
                * ((3 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 2.0 / 5.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 2.0 / 5.0))
                + (2 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 18.0 / 5.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 18.0 / 5.0))
                - (boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 50.0 / 5.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 50.0 / 5.0)));
    }

    // 64QAM理論値
    // 64QAM：1次ダイバーシチ
    double get_64QAMTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);

        return 7.0 / 24.0 * (1.0 - 1.0 / sqrt(1.0 + 7.0 / gamma_b))
                + 1.0 / 4.0 * (1.0 - 1.0 / sqrt(1.0 + 7.0 / 9.0 / gamma_b))
                - 1.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 7.0 / 25.0 / gamma_b))
                + 1.0 / 24.0 * (1.0 - 1.0 / sqrt(1.0 + 7.0 / 81.0 / gamma_b))
                - 1.0 / 24.0 * (1.0 - 1.0 / sqrt(1.0 + 7.0 / 169.0 / gamma_b));
    }

    // 64QAM：L次ダイバーシチ（超幾何関数）
    double get_64QAMTheory_Ldiversity_hyp(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0

        // 以下に64QAMのL次ダイバーシチのBERを返す
        return 1 / (12 * factorial(L - 1) * pow(gamma_c, L)) 
                * ((7 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 1.0 / 7.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 1.0 / 7.0))
                + (6 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 9.0 / 7.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 9.0 / 7.0))
                - (boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 25.0 / 7.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 25.0 / 7.0))
                + (6 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 81.0 / 7.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 81.0 / 7.0))
                - (boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 169.0 / 7.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 169.0 / 7.0)));
    }

    // 256QAM理論値
    // 256QAM：1次ダイバーシチ
    double get_256QAMTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);

        return 15.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / gamma_b))
                + 7.0 / 32.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 9.0 / gamma_b))
                - 1.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 25.0 / gamma_b))
                + 5.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 81.0 / gamma_b))
                + 3.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 121.0 / gamma_b))
                - 1.0 / 16.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 169.0 / gamma_b))
                - 1.0 / 16.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 225.0 / gamma_b))
                + 1.0 / 16.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 289.0 / gamma_b))
                + 1.0 / 16.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 361.0 / gamma_b))
                - 1.0 / 32.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 441.0 / gamma_b))
                - 1.0 / 32.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 529.0 / gamma_b))
                + 1.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 625.0 / gamma_b))
                - 1.0 / 64.0 * (1.0 - 1.0 / sqrt(1.0 + 85.0 / 4.0 / 841.0 / gamma_b));
    }

    // 256QAM：L次ダイバーシチ（超幾何関数）
    double get_256QAMTheory_Ldiversity_hyp(double EbN0dB, int L = 2) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB) / (double)L;
        double gamma_c = EbN0;                  // 伝送路当たりの平均EbN0

        // 以下に256QAMのL次ダイバーシチのBERを返す
        return 1 / (factorial(L - 1) * pow(gamma_c, L)) 
                * ((15.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 4.0 / 85.0))
                + (7.0 / 16.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 9.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 9.0 * 4.0 / 85.0))
                - (1.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 25.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 25.0 * 4.0 / 85.0))
                + (5.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 81.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 81.0 * 4.0 / 85.0))
                + (3.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 121.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 121.0 * 4.0 / 85.0))
                - (1.0 / 8.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 169.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 169.0 * 4.0 / 85.0))
                - (1.0 / 8.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 225.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 225.0 * 4.0 / 85.0))
                + (1.0 / 8.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 289.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 289.0 * 4.0 / 85.0))
                + (1.0 / 8.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 361.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 361.0 * 4.0 / 85.0))
                - (1.0 / 16.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 441.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 441.0 * 4.0 / 85.0))
                - (1.0 / 16.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 529.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 529.0 * 4.0 / 85.0))
                + (1.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 625.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 625.0 * 4.0 / 85.0))
                - (1.0 / 32.0 * boost::math::tgamma(L + 0.5)) / (2 * sqrt(M_PI) * L * pow((1 / gamma_c + 841.0 * 4.0 / 85.0), L))
                * hypergeometric_2F1(0.5, L, L + 1, (1 / gamma_c) / (1 / gamma_c + 841.0 * 4.0 / 85.0)));
    }


    protected:
    int numberOfSymbols_;       // シンボル数
    double P;       // ビット当たりの送信電力

    // データ
    std::vector<int> num_;      // 数値データベクトル
    std::vector<int> grayNum_;      // グレイ符号ベクトル
    Eigen::VectorXd phi_;       // 位相ベクトル
    Eigen::Matrix2d rotationMatrix_;      // 回転行列
    Eigen::VectorXcd symbol_;       // シンボルベクトル
    Eigen::VectorXcd rotatedSymbol_;       // 回転シンボルベクトル

    // BERの上界
    Eigen::MatrixXd alpha0_;      // α_(i; i_hat, 0) = X_(i_hat, 0) - X_(i, 0)
    Eigen::MatrixXd alpha1_;      // α_(i; i_hat, 1) = X_(i_hat, 1) - X_(i, 1)
    Eigen::MatrixXd omega0_;       // ω_(i; i_hat, 0) = {α_(i; i_hat, 0)cosθ - α_(i; i_hat, 1)sinθ}^2
    Eigen::MatrixXd omega1_;       // ω_(i; i_hat, 1) = {α_(i; i_hat, 1)sinθ + α_(i; i_hat, 1)cosθ}^2
    Eigen::MatrixXd kappa_;     // κ_(i; i_hat) = (iとi_hatのハミング距離) / K(ビット数)
    Eigen::MatrixXd gamma0_;        // γ_(i; i_hat, 0) = ω_(i; i_hat, 0) / 2σ^2
    Eigen::MatrixXd gamma1_;        // γ_(i; i_hat, 1) = ω_(i; i_hat, 1) / 2σ^2
    Eigen::MatrixXd zeta0_;        // ζ_(i; i_hat, 0) = γ_(i; i_hat, 0) / γ_(i; i_hat, 0) - γ_(i; i_hat, 1)
    Eigen::MatrixXd zeta1_;        // ζ_(i; i_hat, 1) = γ_(i; i_hat, 1) / γ_(i; i_hat, 1) - γ_(i; i_hat, 0)

    // 通信
    Eigen::Vector2i txData_;     // 送信データベクトル
    Eigen::Vector2i rxData_;     // 受信データベクトル
    double noiseSD_;        // 雑音の標準偏差
    int N_Tri;      // 試行回数

    Eigen::Vector2cd y_;     // 受信信号ベクトル
    Eigen::Vector2cd s_;     // 送信信号ベクトル
    Eigen::Vector2d x0_;      // 回転型QAMシンボルベクトルx[0]
    Eigen::Vector2d x1_;      // 回転型QAMシンボルベクトルx[1]
    Eigen::Matrix2cd h_;     // 伝送路応答行列
    Eigen::Vector2cd n_;     // 雑音ベクトル

    // 乱数用
    unsigned long int seed = 10;      // seed値
    /*Normal Distribution*/
    normal_distribution<> unitNormalRand_;      // 正規乱数
    uniform_int_distribution<> unitIntUniformRand_;     // int型一様乱数
    cnormal_distribution<> unitCNormalRand_;        // 複素正規乱数

    // 距離
    Eigen::VectorXd distance_;      // 送信シンボルと受信シンボルのユークリッド距離ベクトル


    // グレイ符号化
    int setGrayCode(int num) {
        num = num ^ (num >> 1);
        return num;
    }

    // 数値データとグレイ符号のデータのセット
    void setNum() {
        for (int i = 0; i < numberOfSymbols_; i++) {
            num_[i] = i;
            grayNum_[i] = setGrayCode(num_[i]);
        }
    }

    // 伝送路
    // 選択性フェージング伝送路の初期化
    void initSelective() {
        h_.setZero();
        h_(0, 0) = unitCNormalRand_();
        h_(1, 1) = unitCNormalRand_();
        //std::cout << h << std::endl;
    }

    // 送信信号
    void set_x() {
        // 送信データを設定
        txData_(0) = unitIntUniformRand_();
        txData_(1) = unitIntUniformRand_();

        // 送信シンボルを設定
        x0_(0) = rotatedSymbol_(txData_(0)).real();
        x0_(1) = rotatedSymbol_(txData_(0)).imag();
        x1_(0) = rotatedSymbol_(txData_(1)).real();
        x1_(1) = rotatedSymbol_(txData_(1)).imag();

        // 送信信号を設定
        s_(0).real(x0_(0));
        s_(0).imag(x1_(0));
        s_(1).real(x0_(1));
        s_(1).imag(x1_(1));
    }

    // 雑音
    void set_n() {
        n_(0) = unitCNormalRand_();
        n_(1) = unitCNormalRand_();
        //std::cout << noiseSD_ << "*" << n << "=" << noiseSD_ * n << std::endl;
    }

    // 受信
    void set_y() {
        y_ = h_ * s_ + noiseSD_ * n_;
    }

    // 最尤推定で復調
    void set_rxData() {
        Eigen::VectorXd objective_0(numberOfSymbols_);        // 目的関数ベクトル第0項
        Eigen::VectorXd objective_1(numberOfSymbols_);        // 目的関数ベクトル第1項

        Eigen::Vector2d x0;      // 回転型QAMシンボルベクトルx[0]の候補
        Eigen::Vector2d x1;      // 回転型QAMシンボルベクトルx[1]の候補

        for(int i = 0; i < numberOfSymbols_; i++) {
            x0(0) = rotatedSymbol_(i).real();
            x0(1) = rotatedSymbol_(i).imag();
            x1(0) = rotatedSymbol_(i).real();
            x1(1) = rotatedSymbol_(i).imag();

            objective_0(i) = (h_ * ((h_.inverse() * y_).real() - x0)).norm();
            objective_1(i) = (h_ * ((h_.inverse() * y_).imag() - x1)).norm();
        }
        
        Eigen::VectorXd::Index col0;
        Eigen::VectorXd::Index col1;
        objective_0.minCoeff(&col0);
        objective_1.minCoeff(&col1);
        rxData_(0) = col0;
        rxData_(1) = col1;
    }

    // ハミング距離計算
    int hammingDistance(int num1, int num2) {
        int ham = 0;
        int xorResult;
        int bitMask = 1;

        xorResult = num1 ^ num2;

        for(int i = 0; i < NUMBER_OF_BIT; i++) {
            ham += (xorResult & bitMask) >> i;
            bitMask <<= 1;
        }
        
        return ham;
    }

    // ビット誤り数をカウント
    int getBitErrorCount() {
        /* 伝送路の初期化
        *  AWGN 伝送路：initAWGN()
        *  フラットフェージング伝送路：initFlat()
        *  選択性フェージング伝送路：initSelective()
        */
        set_x();     // 送信信号ベクトル生成
        initSelective();        // 伝送路応答行列を用意
        set_n();        // 雑音生成
        set_y();        // 受信
        set_rxData();        //復調

        // グレイ符号化
        rxData_(0) = setGrayCode(rxData_(0));
        rxData_(1) = setGrayCode(rxData_(1));
        txData_(0) = setGrayCode(txData_(0));
        txData_(1) = setGrayCode(txData_(1));

        // ハミング距離計算
        return hammingDistance(rxData_(0), txData_(0)) + hammingDistance(rxData_(1), txData_(1));
    }

    // コンビネーションを計算する関数(nCk)
    double combination(int n, int k) {
        if (k > n) return 0;
        return (double)(std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1)));
    }

    // 階乗
    int factorial(int n) {
        if (n == 0) return 1;
        return n * factorial(n - 1);
    }

    // iとi_hatでループしてα，ω，γ，ζ，κの行列に代入
    void set_alphaOmegaKappaGammaZeta(double theta) {
        // シンボル数 × シンボル数
        alpha0_.resize(numberOfSymbols_, numberOfSymbols_);
        alpha1_.resize(numberOfSymbols_, numberOfSymbols_);
        omega0_.resize(numberOfSymbols_, numberOfSymbols_);
        omega1_.resize(numberOfSymbols_, numberOfSymbols_);
        gamma0_.resize(numberOfSymbols_, numberOfSymbols_);
        gamma1_.resize(numberOfSymbols_, numberOfSymbols_);
        zeta0_.resize(numberOfSymbols_, numberOfSymbols_);
        zeta1_.resize(numberOfSymbols_, numberOfSymbols_);       
        kappa_.resize(numberOfSymbols_, numberOfSymbols_);

        // iとi_hatでループしてα，ω，γ(平均)，ζ，κの行列に代入
        for(int i = 0; i < numberOfSymbols_; i++) {
            for(int i_hat = 0; i_hat < numberOfSymbols_; i_hat++) {
                alpha0_(i, i_hat) = symbol_(i_hat).real() - symbol_(i).real();
                alpha1_(i, i_hat) = symbol_(i_hat).imag() - symbol_(i).imag();
                omega0_(i, i_hat) = pow(alpha0_(i, i_hat) * cos(theta) - alpha1_(i, i_hat) * sin(theta), 2);
                omega1_(i, i_hat) = pow(alpha0_(i, i_hat) * sin(theta) + alpha1_(i, i_hat) * cos(theta), 2);
                gamma0_(i, i_hat) = omega0_(i, i_hat) / 2 / pow(noiseSD_, 2);
                gamma1_(i, i_hat) = omega1_(i, i_hat) / 2 / pow(noiseSD_, 2);

                if(i == i_hat) {
                    kappa_(i, i_hat) = 0;
                    zeta0_(i, i_hat) = 0;
                    zeta1_(i, i_hat) = 0;
                } else {
                    kappa_(i, i_hat) = (double)hammingDistance(grayNum_[i], grayNum_[i_hat]) / (double)NUMBER_OF_BIT;
                    zeta0_(i, i_hat) = gamma0_(i, i_hat) / (gamma0_(i, i_hat) - gamma1_(i, i_hat));
                    zeta1_(i, i_hat) = gamma1_(i, i_hat) / (gamma1_(i, i_hat) - gamma0_(i, i_hat));
                }
            }
        }
    }

    // 数値微分の実装（中央差分法）
    double numerical_derivative(std::function<double(double)> func, double x, double h = 1e-6) {
        return (func(x + h) - func(x - h)) / (2 * h);
    }

    // 二階微分
    double numerical_second_derivative(std::function<double(double)> func, double x, double h = 1e-6) {
        return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
    }

    // Q関数
    double Q(double x) {
        return 1.0 / 2.0 * std::erfc(x / sqrt(2.0));
    }

    // 数値積分（台形則）
    // [min, max]: 積分範囲，N_div: 分割数，func: 関数
    double integration(double min, double max, int N_div, double(*func)(double)) {
        double dx = (max - min) / N_div;  // 区間幅
        double sum = 0.0;

        // 積分
        for (int i = 0; i <= N_div; ++i) {
            double x = min + i * dx;
            if (i == 0 || i == N_div) {
                sum += func(x);
            } else {
                sum += 2 * func(x);
            }
        }

        return sum * dx / 2.0;
    }

    // ガウスの超幾何関数2F1(a, b; c; z)
    double hypergeometric_2F1(double a, double b, double c, double z) {
        std::initializer_list<double> aj = {a, b};
        std::initializer_list<double> bj = {c};

        try {
            return boost::math::hypergeometric_pFq(aj, bj, z);
        } catch (const std::exception& e) {
            std::cerr << "Error computing hypergeometric_2F1: " << e.what() << std::endl;
            throw;
        }
    }
};
#endif /* SIMULATOR_H */