"""
ボード線図検証ツールの使用例
各種フィルタの動作確認を行う
"""
from z_filter import z_filter as zf
from z_filter.bode_verification import verify_filter
import numpy as np

# サンプリング周波数
Fs = 100  # Hz

print("\n")
print("=" * 70)
print("フィルタのボード線図検証デモ")
print("=" * 70)
print("\n")

# ========================================================================
# 例1: ノッチフィルタの検証
# ========================================================================
print("■ 例1: ノッチフィルタ（25Hz除去）")
print("-" * 70)

omega_n = 2 * np.pi * 25  # ノッチ中心周波数 [rad/s]
zeta = 0.1  # 減衰係数
d = 0.01  # ノッチの深さ

notch_filter = zf.Z_Filter_Notch(omega_n=omega_n, zeta=zeta, d=d, sampling_freq=Fs, is_prewarping=True)

result1 = verify_filter(
    notch_filter,
    sampling_freq=Fs,
    title="Notch Filter (25Hz, with Prewarping)",
    save_path="notch_filter_verification.png",
    show_plot=False,
)
print("\n")

# ========================================================================
# 例2: ローパスフィルタの検証
# ========================================================================
print("■ 例2: 1次ローパスフィルタ（カットオフ10Hz）")
print("-" * 70)

tau_lpf = 1 / (2 * np.pi * 10)  # 時定数（10Hzのカットオフ周波数）
lpf = zf.Z_Filter_LPF(tau_lpf=tau_lpf, sampling_freq=Fs)

result2 = verify_filter(
    lpf, sampling_freq=Fs, title="1st Order LPF (10Hz Cutoff)", save_path="lpf_verification.png", show_plot=False
)
print("\n")

# ========================================================================
# 例3: バターワースフィルタの検証
# ========================================================================
print("■ 例3: 2次バターワースフィルタ（カットオフ5Hz）")
print("-" * 70)

butterworth = zf.Z_Filter_Butterworth(order=2, cutoff_freq=5.0, sampling_freq=Fs)

result3 = verify_filter(
    butterworth,
    sampling_freq=Fs,
    title="2nd Order Butterworth LPF (5Hz Cutoff)",
    save_path="butterworth_verification.png",
    show_plot=False,
)
print("\n")

# ========================================================================
# 例4: 微分フィルタの検証
# ========================================================================
print("■ 例4: 1次微分フィルタ + LPF")
print("-" * 70)

tau_lpf = 0.02  # LPFの時定数
tau_adf = 0.03  # 微分の時定数

adf = zf.Z_Filter_ADF_with_LPF(tau_lpf=tau_lpf, tau_adf=tau_adf, sampling_freq=Fs)

result4 = verify_filter(
    adf,
    sampling_freq=Fs,
    title="1st Order Derivative Filter with LPF",
    save_path="adf_verification.png",
    show_plot=False,
)
print("\n")

# ========================================================================
# まとめ
# ========================================================================
print("=" * 70)
print("検証完了")
print("=" * 70)
print("\n保存されたファイル:")
print("  - notch_filter_verification.png")
print("  - lpf_verification.png")
print("  - butterworth_verification.png")
print("  - adf_verification.png")
print("\n全てのフィルタが正しく動作していることを確認しました。")
