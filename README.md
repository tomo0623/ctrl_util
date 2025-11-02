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

# ログ出力を無効にして静かに戻す
zf.disable_verbose_logging()
```

### 利用可能なフィルタクラス

- `Z_Filter`: 汎用離散時間フィルタ（分子・分母係数を直接指定）
- `Z_Filter_LPF`: 1次ローパスフィルタ
- `Z_Filter_ADF`: 1次微分フィルタ
- `Z_Filter_ADF_with_LPF`: 微分フィルタとLPFの組み合わせ
- `Z_Filter_Butterworth`: N次バターワースローパスフィルタ（SciPy不要）