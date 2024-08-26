/*
 * File:   Simulator.h
 * Author: Tsunakawa
 *
 * Created on 2023/04/29, 15:49
*/

#ifndef SIMULATOR_H
#define SIMULATOR_H
#define _USE_MATH_DEFINES
#include </Volumes/USB1/eigen-3.4.0/Eigen/Dense>
#include </Volumes/USB1/eigen-3.4.0/Eigen/Eigen>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <bitset>
#include <functional>
#include "random_collection.h"

class Simulator {
    public:
    int NUMBER_OF_BIT;      // ビット数(QAM方式)

    // コンストラクタ
    Simulator() {
        // 次元を設定
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "Number of Bit? (QPSK:2, 16QAM:4)" << std::endl;
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
                        case 0b00: 
                            symbol_(i).real(-3 * sqrt(P));
                        break;
                        case 0b01: 
                            symbol_(i).real(-sqrt(P));
                        break;
                        case 0b11: 
                            symbol_(i).real(sqrt(P));
                        break;
                        case 0b10: 
                            symbol_(i).real(3 * sqrt(P));
                        break;
                    }
                    // 虚部の設計(1, 2ビット目で決まる)
                    switch(grayNum_[i] & 0b11) {
                        case 0b00: 
                            symbol_(i).imag(-3 * sqrt(P));
                        break;
                        case 0b01: 
                            symbol_(i).imag(-sqrt(P));
                        break;
                        case 0b11: 
                            symbol_(i).imag(sqrt(P));
                        break;
                        case 0b10: 
                            symbol_(i).imag(3 * sqrt(P));
                        break;
                    }
                }
            break;
        }
        //std::cout << symbol_ << std::endl;
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

    // Walsh Hadamard行列後のシンボル
    /*
    void setSymbol_WalshHadamard(){
        for(int i = 0; i < numberOfSymbols_; i++) {
            symbol_(i).real(M_SQRT1_2 * (cos(phi_(i)) + sin(phi_(i))));
            symbol_(i).imag(M_SQRT1_2 * (cos(phi_(i)) - sin(phi_(i))));
        }
    }
    */

    // 加法性雑音の標準偏差を設定
    void set_QPSKNoiseSD(double EbN0dB) {
        noiseSD_ = sqrt(0.5 * pow(10.0, -0.1 * EbN0dB));        // Eb/N0 [dB] から変換
    }

    void set_16QAMNoiseSD(double EbN0dB) {
        noiseSD_ = sqrt(0.25 * pow(10.0, -0.1 * EbN0dB));        // Eb/N0 [dB] から変換
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


    // 理論値
    // QPSK平均BERの上界
    double getBerUpperBound(double theta) {
        double berUpperBound = 0;
        double per = 0;
        Eigen::MatrixXd gamma;
        double zata_0;
        double zeta_1;

        // γ_0, γ_1
        gamma.resize(4, 2);     // γは誤りパターンの4通り存在する
        gamma(0, 0) = pow(cos(theta), 2);        // -+|++など4個
        gamma(0, 1) = pow(sin(theta), 2);
        gamma(1, 0) = pow(cos(theta), 2);        // +-|++など4個
        gamma(1, 1) = pow(sin(theta), 2);
        gamma(2, 0) = (1 - sin(2 * theta));        // --|++など2個
        gamma(2, 1) = (1 + sin(2 * theta));
        gamma(3, 0) = (1 + sin(2 * theta));        // -+|+-など2個
        gamma(3, 1) = (1 - sin(2 * theta));
        gamma = 2 * P / pow(noiseSD_, 2) * gamma;

        for(int i = 0; i < 4; i++) {
            zata_0 = gamma(i, 0) / (gamma(i, 0) - gamma(i, 1));
            zeta_1 = gamma(i, 1) / (gamma(i, 1) - gamma(i, 0));
            per = ((zata_0 * (1 - sqrt(gamma(i, 0) / (2 + gamma(i, 0))))) + (zeta_1 * (1 - sqrt(gamma(i, 1) / (2 + gamma(i, 1)))))) / 2;

            // 誤り方によって場合分け
            if (i == 0 || i == 1) {
                per = per * 4 / 2;      // 4個 × 重み(1/2)
            } else {
                per = per * 2;      // 2個 × 重み(1)
            }
            berUpperBound += per;
        }

        return berUpperBound / 4;
    }

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

    // 平均BERの上界の近似(QPSK)
    double nearly_getQPSKBerUpperBound(double theta) {
        double lambda = 1 / (pow(cos(theta), 2) * pow(sin(theta), 2)) + 1 / pow(cos(2 * theta), 2);
        return 3 * pow(noiseSD_, 4) * lambda / 4;
    }

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
    double get_optimizeTheta_rad_byNewton(std::function<double(double)> func, double initial_guess, int max_iter = 100, double tol = 1e-12) {
        double x = initial_guess;
        for (int i = 0; i < max_iter; ++i) {
            // 更新
            x = x - numerical_derivative(func, x) / numerical_second_derivative(func, x);
            
            // 収束判定
            if (std::abs(numerical_derivative(func, x)) < tol) {
                std::cout << "------------" << i << " times search------------" << std::endl;
                break;
            }
        }
        return x;
    }

    // 解を求める（気合）
    // resolution：分解能
    double get_optimizeTheta_deg(std::function<double(double)> func, double start_deg, double end_deg, double resolution_deg) {
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

    // QPSK_1次ダイバーシチ[ディジタルコミュニケーションp.898:式(14-3-7)]
    double get_QPSKTheory(double EbN0dB) {
        double gamma_b = pow(10.0, 0.1 * EbN0dB);
        double ber = (1 - sqrt(gamma_b / (1.0 + gamma_b))) / 2;
        return ber;
    }

    // QPSK_2次のダイバーシチ[ディジタルコミュニケーションp.898:式(14-4-15)]
    double get_QPSKTheory_2diversity(double EbN0dB) {
        double EbN0 = 0.5 * pow(10.0, 0.1 * EbN0dB);
        double mu = sqrt(EbN0 / (1.0 + EbN0));
        double ber = pow((1.0 - mu) / 2.0, 2) * (2.0 + mu);
        return ber;
    }

    // 16QAM_1次のダイバーシチ[ディジタルコミュニケーションp.917:式(14-4-47)]
    // 16QAM理論値
    double get16QAMTheory(double EbN0dB) {
        double EbN0 = pow(10.0, 0.1 * EbN0dB);

        return 3.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 / 2.0 / EbN0))
                + 1.0 / 4.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / EbN0 / 9.0))
                - 1.0 / 8.0 * (1.0 - 1.0 / sqrt(1.0 + 5.0 * 2.0 / EbN0 / 25.0));
    }

    // 16QAM_L次のダイバーシチ[ディジタルコミュニケーションp.917:式(14-4-47)]
    // デフォルトでは2次ダイバーシチ
    double get_QAMTheory_Ldiversity(double EbN0dB, int L = 2) {
        double EbN0 = 0.25 * pow(10.0, 0.1 * EbN0dB);
        double mu = sqrt(EbN0 / (1.0 + EbN0));
        double ber = pow((1.0 - mu) / 2.0, 2) * (2.0 + mu);
        return ber;
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

    Eigen::Vector2cd y;     // 受信信号ベクトル
    Eigen::Vector2cd x;     // 送信信号ベクトル
    Eigen::Vector2d s0;      // シンボルベクトルs[0]
    Eigen::Vector2d s1;      // シンボルベクトルs[1]
    Eigen::Matrix2cd h;     // 伝送路応答行列
    Eigen::Vector2cd n;     // 雑音ベクトル

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
        h.setZero();
        h(0, 0) = unitCNormalRand_();
        h(1, 1) = unitCNormalRand_();
        //std::cout << h << std::endl;
    }

    // 送信信号
    void set_x() {
        // 送信データを設定
        txData_(0) = unitIntUniformRand_();
        txData_(1) = unitIntUniformRand_();

        // 送信シンボルを設定
        s0(0) = rotatedSymbol_(txData_(0)).real();
        s0(1) = rotatedSymbol_(txData_(0)).imag();
        s1(0) = rotatedSymbol_(txData_(1)).real();
        s1(1) = rotatedSymbol_(txData_(1)).imag();

        // 送信信号を設定
        x(0).real(s0(0));
        x(0).imag(s1(0));
        x(1).real(s0(1));
        x(1).imag(s1(1));
    }

    // 雑音
    void set_n() {
        n(0) = unitCNormalRand_();
        n(1) = unitCNormalRand_();
        //std::cout << noiseSD_ << "*" << n << "=" << noiseSD_ * n << std::endl;
    }

    // 受信
    void set_y() {
        y = h * x + noiseSD_ * n;
    }

    // 最尤推定で復調
    void set_rxData() {
        Eigen::VectorXd objective_0(numberOfSymbols_);        // 目的関数ベクトル第0項
        Eigen::VectorXd objective_1(numberOfSymbols_);        // 目的関数ベクトル第1項

        for(int i = 0; i < numberOfSymbols_; i++) {
            s0(0) = rotatedSymbol_(i).real();
            s0(1) = rotatedSymbol_(i).imag();
            s1(0) = rotatedSymbol_(i).real();
            s1(1) = rotatedSymbol_(i).imag();

            objective_0(i) = (h * ((h.inverse() * y).real() - s0)).norm();
            objective_1(i) = (h * ((h.inverse() * y).imag() - s1)).norm();
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

    // コンビネーションを計算する関数
    double combination(int n, int k) {
        if (k > n) return 0;
        return std::tgamma(n + 1) / (std::tgamma(k + 1) * std::tgamma(n - k + 1));
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
};
#endif /* SIMULATOR_H */