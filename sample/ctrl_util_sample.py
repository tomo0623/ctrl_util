from z_filter import z_filter as zf
import numpy as np
from matplotlib import pyplot as plt


# 基本的なZフィルタの動作確認テスト

# ---------------------- 正弦波信号の生成
Fs = 100  # サンプリング周波数
max_time = 10.0  # 信号の長さ（秒）
t = np.arange(0, max_time, 1 / Fs)
signal = np.sin(2 * np.pi * 1 * t) + 10
dot_signal = 2 * np.pi * 1 *np.cos(2 * np.pi * 1 * t)
# ノイズ付き信号の生成
noisy_signal = signal + 0.1 * np.random.randn(len(t))


# ---------------------- フィルタリングテスト
# Zフィルタの初期化
tau_lpf = 0.02  # LPFの時定数
tau_adf = 0.03  # ADFの時定数
Fs = 100  # サンプリング周波数

print("フィルタ1: 単純な1次LPF")
c = np.array([1])
a = np.array([tau_lpf, 1])
filter1 = zf.Z_Filter(a, c, sampling_freq=Fs)

print("フィルタ2: 1次微分フィルタ")
c = np.array([1, 0])
a = np.array([tau_lpf, 1])
filter2 = zf.Z_Filter(a, c, sampling_freq=Fs)

print("フィルタ3: LPF付き1次微分フィルタ")
c = np.array([1/(tau_lpf*tau_adf), 0])
a = np.array([1,(tau_lpf+tau_adf)/(tau_lpf*tau_adf), 1/(tau_lpf*tau_adf)])
filter3 = zf.Z_Filter(a, c, sampling_freq=Fs)

print("フィルタクラスのテスト")
filter4 = zf.Z_Filter_ADF_with_LPF(0.03, 0.02, sampling_freq=Fs)

print("2次のバターワースフィルタ")
filter5 = zf.Z_Filter_Butterworth(order=2, cutoff_freq=5.0, sampling_freq=Fs)

# フィルタリング処理(時系列処理の模擬)
loggs = np.zeros((len(t), 5))

for i in range(len(t)):
    if i < 10:
        # リセットcommandでフィルタを現在値に初期化
        filter1.reset(noisy_signal[i])
        filter2.reset(0, u_ss=noisy_signal[i])  # 微分出力を0にリセット
        filter3.reset(0, u_ss=noisy_signal[i])  # 微分出力を0にリセット
        filter4.reset(0, u_ss=noisy_signal[i])  # 微分出力を0にリセット
        filter5.reset(noisy_signal[i])
    loggs[i, 0] = filter1.update(noisy_signal[i])
    loggs[i, 1] = filter2.update(noisy_signal[i])
    loggs[i, 2] = filter3.update(noisy_signal[i])
    loggs[i, 3] = filter4.update(noisy_signal[i])
    loggs[i, 4] = filter5.update(noisy_signal[i])

# ---------------------- 結果の可視化
# 2段のグラフを作成（X軸をリンク）
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# 上段：ノイズ付き信号と元信号
ax1.plot(t, noisy_signal, label="Noisy Signal")
ax1.plot(t, signal, label="Original Signal", linestyle="--")
ax1.plot(t, loggs[:, 0], label="Filtered Signal 1", linestyle="--")
ax1.plot(t, loggs[:, 4], label="Filtered Signal 5", linestyle="--")
ax1.set_ylabel("Amplitude")
ax1.set_title("Noisy Signal vs Original Signal")
ax1.legend()
ax1.grid()

# 下段：微分信号
ax2.plot(t, dot_signal, label="Derivative Signal", color="red")
ax2.plot(t, loggs[:, 1], label="Filtered Signal 2", linestyle="--")
# ax2.plot(t, loggs[:, 2], label="Filtered Signal 3", linestyle="--")
ax2.plot(t, loggs[:, 3], label="Filtered Signal 4", linestyle="--")
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("Amplitude")
ax2.set_title("Derivative Signal")
ax2.legend()
ax2.grid()

plt.tight_layout()  # レイアウトを調整
plt.show()
