# rotation_K=2

## 概要
二次元平面上で回転させたときのダイバーシチについて検討するプログラム．
極座標での表現


## 進捗
- Rotated QPSK OK
- QPSKや16QAMでBER近似が算出
- 厳密な上界の算出
- ニュートン法により最適な回転角導出
- QPSKと16QAMの1次ダイバーシチ
- QPSKと16QAMの2次ダイバーシチ

## 今後の課題
- 一般化
    - [x] 横軸 $\theta$ の上界
    - [x] 最適な $\theta$ の算出
    - [x] 横軸 $E_b/N_0$ の上界
    - [x] L次ダイバーシチ（hyp）
    - [x] L次ダイバーシチ（sim）
    - [x] 横軸 $\theta$ のシミュレーション
    - [x] 横軸 $E_b/N_0$ のシミュレーション
    - [x] 1次ダイバーシチ（64QAM）
    - [x] 2次ダイバーシチ（64QAM）
- 結果の取得
  - [x] 4QAM
    - 理論値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$
        - 1次
        - 2次
    - シミュレーション値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$
  - [x] 16QAM
    - 理論値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$
        - 1次
        - 2次
    - シミュレーション値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$
  - [ ] 64QAM
    - [x] 理論値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$
        - 1次
        - 2次
    - [ ] シミュレーション値
      - 横軸 $\theta$
      - 横軸 $E_b/N_0$

## boostライブラリの使い方？
1. boostをインストール
2. includeする
3. c_cpp_properties.jsonにインクルードパスを書く
4. ```g++ -std=c++17 -IC:\boost_1_86_0 -o diversity_hyp diversity_hyp.cpp```でコンパイル</br>
もし，エラーが出たらおそらくmutexに関することなので`boost_1_86_0/boost/math/special_functions/detail/polygamma.hpp`の409~412行目をコメントアウトしておく（boost1_86_0の場合）
5. `./diversity_hyp`でdiversity_hyp.exeを実行

## MTG
- 10/2</br>
  - 既に出ている理論値，シミュレーション，数値積分，ガウスの超幾何関数を用いた式，合計4つがL=1, 2それぞれで一致しているかを確認
  - 理由付けは微妙
  - 4QAM~256QAMまで全部30dBを比較するように図示
- 10/18</br>
  - 64QAMまでで良さそう
  - 修論をジャーナルに出す