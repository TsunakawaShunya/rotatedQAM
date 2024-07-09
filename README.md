# rotation_K=2

## 概要
二次元平面上で回転させたときのダイバーシチについて検討するプログラム．
極座標での表現


## 進捗
- Rotated QPSK OK
- QPSKや16QAMでBER近似が算出
- 厳密な上界の算出
- ニュートン法により最適な回転角導出

## 今後の課題
- 1次ダイバーシチ，シミュレーション値，厳密な上界，2次ダイバーシチを引く
    - QPSKと4QAM1次ダイバーシチで合っていない（4QAMが間違えている）
        - 手計算でも合わないから先生に聞きに行く
    - 2次ダイバーシチの一般化