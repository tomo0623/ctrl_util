# ctrl_util

制御工学ユーティリティライブラリ

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

### 利用可能なフィルタクラス

- `Z_Filter`: 汎用離散時間フィルタ（分子・分母係数を直接指定）
- `Z_Filter_LPF`: 1次ローパスフィルタ
- `Z_Filter_ADF`: 1次微分フィルタ
- `Z_Filter_ADF_with_LPF`: 微分フィルタとLPFの組み合わせ
- `Z_Filter_Butterworth`: N次バターワースローパスフィルタ（SciPy不要）