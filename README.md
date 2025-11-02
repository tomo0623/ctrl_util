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
