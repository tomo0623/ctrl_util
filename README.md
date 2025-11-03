# ctrl_util

制御工学ユーティリティライブラリ

## インストール

### 方法1: GitHubから直接インストール

```bash
pip install git+https://github.com/yourusername/ctrl_util.git
```

### 方法2: インストールせずに使う（ファイルコピー）

`z_filter` フォルダをプロジェクトにコピーするだけで使えます：

```
your_project/
├── z_filter/          # このフォルダをコピー
│   ├── __init__.py
│   └── z_filter.py
└── your_script.py
```

```python
from z_filter import z_filter as zf

lpf = zf.Z_Filter_LPF(tau_lpf=0.1, sampling_freq=100.0)
```

### 方法3: ローカル開発

```bash
# リポジトリをクローン
git clone https://github.com/yourusername/ctrl_util.git
cd ctrl_util

# 開発モードでインストール
pip install -e .

# (uvの場合)開発モードでインストール
uv pip install -e .
```

## フィルタ機能 (z_filter)

### 基本的な使用方法

```python
from z_filter import z_filter as zf

# 1次LPF（時定数0.1秒、サンプリング周波数100Hz）
lpf = zf.Z_Filter_LPF(tau_lpf=0.1, sampling_freq=100.0)

# 入力信号を処理
output = lpf.update(1.0)
```

### バターワースフィルタ

```python
# 2次バターワースLPF（カットオフ10Hz、サンプリング周波数100Hz）
butterworth = zf.Z_Filter_Butterworth(order=2, cutoff_freq=10.0, sampling_freq=100.0)

# フィルタ処理
for input_val in input_data:
    filtered_output = butterworth.update(input_val)
```

### 利用可能なフィルタクラス

汎用的なフィルタクラスに加え、よく使うフィルタクラス(LPF、微分フィルタなど)は時定数やカットオフ周波数を設定するだけで使えるように用意されています：

- `Z_Filter`: 汎用離散時間フィルタ（分子・分母係数を直接指定）
- `Z_Filter_LPF`: 1次ローパスフィルタ（テンプレート, 時定数・サンプリング周波数を指定）
- `Z_Filter_ADF`: 1次微分フィルタ（テンプレート, 時定数・サンプリング周波数を指定）
- `Z_Filter_ADF_with_LPF`: 微分フィルタとLPFの組み合わせ（テンプレート, 時定数・サンプリング周波数を指定）
- `Z_Filter_Butterworth`: N次バターワースローパスフィルタ（テンプレート, 次数・カットオフ周波数・サンプリング周波数を指定）

### ログ出力の制御

デフォルトでは、フィルタライブラリは何も出力しません。内部の計算詳細を確認したい場合は、以下のようにログ出力を有効にできます：

```python
import logging
from z_filter import z_filter as zf

# デバッグ情報を表示（可制御正準系の行列、離散化行列など）
zf.enable_verbose_logging(logging.DEBUG)

# または、より少ない情報のみ表示
zf.enable_verbose_logging(logging.INFO)

# フィルタを作成・使用（内部計算が表示される）
lpf = zf.Z_Filter_LPF(tau_lpf=0.1, sampling_freq=100.0)

# ログ出力を無効に戻す
zf.disable_verbose_logging()
```

### リセット機能

フィルタの内部状態をリセットして、特定の出力値から開始することができます。

#### 基本的なリセット（LPFなど）

```python
lpf = zf.Z_Filter_LPF(tau_lpf=0.1, sampling_freq=100.0)

# フィルタ出力を10.0にリセット
lpf.reset(10.0)

# 次のupdateで出力は10.0に近い値から始まる
output = lpf.update(10.0)  # 出力 ≈ 10.0
```

#### 微分フィルタのリセット

微分フィルタの場合、定常状態の入力値を `u_ss` パラメータで指定する必要があります。

```python
adf = zf.Z_Filter_ADF(tau_derivative=0.1, sampling_freq=100.0)

# 微分出力を0にリセット（入力が5.0で一定と仮定）
adf.reset(0, u_ss=5.0)

# 次のupdateで微分出力は0に近い値から始まる
output = adf.update(5.0)  # 出力 ≈ 0
```

#### リセットモード

- **`mode="output"`（デフォルト）**: フィルタ出力が指定値になるように内部状態を設定
  - `reset(value)`: 定常状態の入力を `value` と仮定（LPF向け）
  - `reset(value, u_ss=u_val)`: 定常状態の入力を明示的に指定（微分フィルタ向け）
- **`mode="state"`**: 内部状態ベクトルを直接設定（上級者向け）

```python
# 内部状態ベクトルを直接指定
lpf.reset([0.5], mode="state")
```

### 伝達関数の係数から直接フィルタを構成

連続時間系の伝達関数の分母・分子係数から離散フィルタを作成できます。係数を与えるだけで、内部で双一次変換により離散化されます。

```python
import numpy as np

# 例: LPF付き1次微分フィルタ
tau_lpf = 0.1  # LPFの時定数
tau_adf = 0.05  # 微分フィルタの時定数
Fs = 100.0  # サンプリング周波数

# 伝達関数: H(s) = C(s) / A(s)
# $$H(s) = \frac{s/\tau_{lpf}}{1 + \frac{\tau_{lpf}+\tau_{adf}}{\tau_{lpf}\tau_{adf}}s + \frac{1}{\tau_{lpf}\tau_{adf}}s^2}$$

# 分子係数 C(s) = c[0]*s + c[1]  (sの降べきの順)
c = np.array([1/(tau_lpf*tau_adf), 0])

# 分母係数 A(s) = a[0]*s^2 + a[1]*s + a[2]  (sの降べきの順)
a = np.array([1, (tau_lpf+tau_adf)/(tau_lpf*tau_adf), 1/(tau_lpf*tau_adf)])

# フィルタを構成
filter = zf.Z_Filter(a, c, sampling_freq=Fs)

# 使用
output = filter.update(input_val)
```

**係数の指定方法:**

伝達関数を以下の形式で表現します：

$$H(s) = \frac{C(s)}{A(s)} = \frac{c_m s^m + c_{m-1} s^{m-1} + \cdots + c_0}{a_n s^n + a_{n-1} s^{n-1} + \cdots + a_0}$$

- 分母係数 `a`: `[a_n, a_{n-1}, ..., a_0]` （$s$ の降べきの順）
- 分子係数 `c`: `[c_m, c_{m-1}, ..., c_0]` （$s$ の降べきの順）

**プリワーピング処理:**

双一次変換（タスティン変換）では、高周波側で周波数特性が歪みます。プリワーピング処理を使用すると、特定の周波数で連続時間系と離散系の周波数特性を一致させることができます。

```python
# プリワーピングなし（デフォルト）
filter_normal = zf.Z_Filter(a, c, sampling_freq=100.0)

# プリワーピングあり（10Hzで周波数特性を一致させる）
filter_prewarped = zf.Z_Filter(
    a, c,
    sampling_freq=100.0,
    is_prewarping=True,
    prewarping_freq=10.0  # 一致させたい周波数 [Hz]
)
```

プリワーピングは、カットオフ周波数付近の特性を正確に再現したい場合に有用です。

## 加速度座標変換機能 (accel_transform)

剛体内の加速度センサの測定値を、剛体上の別の点での加速度に座標変換する機能。
ロボット工学や慣性航法システムでよく使われる剛体運動学の基本式を実装してる。
計算過程で角加速度が必要になるため、角速度の近似微分を計算しているため、Class初期化時のフィルタ設定に留意すること。

### 理論的背景

剛体上の2点（観測点と目標点）の加速度の関係は、以下の式で表される：

$$\boldsymbol{\alpha}_t = \boldsymbol{\alpha}_o + \dot{\boldsymbol{\omega}} \times \boldsymbol{r} + \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \boldsymbol{r})$$

ここで：
- $\boldsymbol{\alpha}_t$: 目標点の加速度
- $\boldsymbol{\alpha}_o$: 観測点の加速度（センサ測定値）
- $\dot{\boldsymbol{\omega}}$: 角加速度ベクトル（角速度の時間微分）
- $\boldsymbol{\omega}$: 角速度ベクトル（ジャイロセンサ測定値）
- $\boldsymbol{r}$: 観測点から目標点へのベクトル
  - なお、`×`は外積を意味し、すべての変数は3次元ベクトルであることに注意する

**各項の物理的意味：**
1. $\boldsymbol{\alpha}_o$: 観測点での並進加速度
2. $\dot{\boldsymbol{\omega}} \times \boldsymbol{r}$: オイラー加速度（角加速度による接線加速度）
3. $\boldsymbol{\omega} \times (\boldsymbol{\omega} \times \boldsymbol{r})$: 向心加速度（角速度による遠心加速度）

### 基本的な使用方法

```python
from accel_transform import acc_trans as accT
import numpy as np

# 初期化
obs2tar_vec = [0.1, 0.2, 0.3]  # 観測点→目標点のベクトル [m]
sampling_freq = 100.0  # サンプリング周波数 [Hz]
tau_lpf = 0.03  # LPF時定数 [秒]
tau_adf = 0.02  # 微分フィルタ時定数 [秒]

transformer = accT.AccTransform(obs2tar_vec, sampling_freq, tau_lpf, tau_adf)

# リアルタイム処理のループ
for i in range(len(data)):
    # センサデータの取得
    obs_acc = [ax, ay, az]  # 加速度センサ値 [m/s^2]
    ang_vel = [wx, wy, wz]  # ジャイロセンサ値 [rad/s]

    # 座標変換の実行
    tar_acc, ang_acc = transformer.update(obs_acc, ang_vel)

    # tar_acc: 目標点での加速度 [m/s^2]
    # ang_acc: 角加速度（参考値）[rad/s^2]
```

### パラメータの説明

#### コンストラクタ `AccTransform(obs2tar_vec, sampling_freq, tau_lpf, tau_adf)`

- **`obs2tar_vec`**: 観測点から目標点へのベクトル `[x, y, z]` [m]
  - 加速度センサが取り付けられている位置から、加速度を知りたい目標位置までの相対位置ベクトル
  - 右手座標系で定義

- **`sampling_freq`**: サンプリング周波数 [Hz]
  - センサデータの更新レート

- **`tau_lpf`**: LPF（ローパスフィルタ）時定数 [秒]
  - 角加速度計算時の微分フィルタに使用
  - 小さいほど応答が速いが、ノイズの影響を受けやすい

- **`tau_adf`**: ADF（近似微分フィルタ）時定数 [秒]
  - 角加速度計算時のノイズ除去に使用
  - 小さいほど精度が高いが、ノイズの影響を受けやすい

**フィルタ時定数の設定ガイドライン:**
- 最小値: `サンプリング周期 × 2` 以上
- 推奨値: 0.01～0.1秒程度（100Hzサンプリングの場合）
- ノイズが大きい場合は時定数を大きく、応答性を優先する場合は小さく設定

#### `update(obs_acc, ang_vel)` メソッド

- **入力:**
  - `obs_acc`: 観測点（センサ位置）での加速度ベクトル `[ax, ay, az]` [m/s^2]
  - `ang_vel`: 剛体の角速度ベクトル `[wx, wy, wz]` [rad/s]

- **出力:**
  - `tar_acc`: 目標点での加速度ベクトル `[ax, ay, az]` [m/s^2]
  - `ang_acc`: 角加速度ベクトル `[alpha_x, alpha_y, alpha_z]` [rad/s^2]（参考値）

### 応用例

#### 例1: ロボットアームの先端加速度計算

```python
# ロボットアームの根元にセンサがあり、先端の加速度を知りたい場合
obs2tar_vec = [0.0, 0.0, 0.5]  # アーム長0.5m（Z方向）
transformer = accT.AccTransform(obs2tar_vec, 100.0, 0.03, 0.02)

# センサデータから先端加速度を計算
tip_acc, _ = transformer.update(base_acc, gyro_data)
```

#### 例2: 車両の重心加速度推定

```python
# センサが車両前方に設置されており、重心の加速度を推定
obs2tar_vec = [-0.5, 0.0, -0.1]  # 重心方向へのベクトル
transformer = accT.AccTransform(obs2tar_vec, 100.0, 0.05, 0.03)

# 車両の角速度とセンサ加速度から重心加速度を計算
cg_acc, _ = transformer.update(sensor_acc, gyro_data)
```

### ログ出力の制御

`z_filter` と同様に、デフォルトでは何も出力しない。内部計算の詳細を確認したい場合は：

```python
import logging
from accel_transform import acc_trans as accT

# デバッグ情報を表示
accT.enable_verbose_logging(logging.DEBUG)

# 変換処理を実行（内部計算が表示される）
transformer = accT.AccTransform([0.1, 0.2, 0.3], 100.0, 0.03, 0.02)

# ログ出力を無効に戻す
accT.disable_verbose_logging()
```

### 注意事項

1. **座標系の統一**: 加速度ベクトル、角速度ベクトル、変換ベクトルは全て同じ座標系（右手座標系）で定義する必要がある
2. **重力加速度の扱い**: 加速度センサの出力には重力加速度が含まれています。必要に応じて重力成分を除去してから座標変換を行うこと（仮想的に加速度センサ搭載個所を変更したときの加速度センサ値を得たい場合は、センサ値をそのまま変換してよい）
3. **角加速度の精度**: 角速度を数値微分して角加速度を計算するため、ノイズの影響を受けやすくなります。適切なフィルタ時定数の設定が重要となる（もしノイズがひどい場合は、ADFはサンプリング周期の倍の時定数で固定し、LPF時定数のみを弄るのがおすすめ）
4. **サンプリング周波数**: 高速な運動を扱う場合は、サンプリング周波数を十分に高く設定してください（サンプリング定理に従う）
