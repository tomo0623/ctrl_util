from z_filter import z_filter as zf
from accel_transform import acc_trans as accT
import numpy as np
from matplotlib import pyplot as plt



# ----------------------------------------------------------------------- 基本的なZフィルタの動作確認テスト

if 1:  
    # ---------------------- 正弦波信号の生成
    Fs = 100  # サンプリング周波数
    max_time = 10.0  # 信号の長さ（秒）
    t = np.arange(0, max_time, 1 / Fs)
    signal = np.sin(2 * np.pi * 1 * t) + 10
    dot_signal = 2 * np.pi * 1 * np.cos(2 * np.pi * 1 * t)
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
    c = np.array([1 / (tau_lpf * tau_adf), 0])
    a = np.array([1, (tau_lpf + tau_adf) / (tau_lpf * tau_adf), 1 / (tau_lpf * tau_adf)])
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



# ----------------------------------------------------------------------- 基本的な加速度座標変換の動作確認テスト
if 1:  
    # ---------------------- 正弦波信号の生成
    Fs = 100  # サンプリング周波数
    max_time = 10.0  # 信号の長さ（秒）
    t = np.arange(0, max_time, 1 / Fs)
    signal = np.sin(2 * np.pi * 1 * t)
    dot_signal = 2 * np.pi * 1 * np.cos(2 * np.pi * 1 * t)
    # ノイズ付き信号の生成(角速度センサ値模擬)
    noisy_signal1 = (signal + 0.1 * np.random.randn(len(t))) * 1
    noisy_signal2 = (2 * signal + 0.1 * np.random.randn(len(t))) * 1
    noisy_signal3 = (3 * signal + 0.1 * np.random.randn(len(t))) * 1
    # 加速度センサ信号ベース
    base_signal = 0.5 * np.sin(t) + 0.05 * np.random.randn(len(t))

    # ---------------------- 座標変換テスト
    # 座標変換の初期化
    obs2tar_vec = [0.5, 0.5, 0.5]  # 観測座標系から目標座標系への変換ベクトル
    tau_lpf = 0.03  # LPFの時定数
    tau_adf = 0.02  # ADFの時定数

    acc_transformer = accT.AccTransform(obs2tar_vec, Fs, tau_lpf, tau_adf)

    # フィルタリング処理(時系列処理の模擬)
    loggs = np.zeros((len(t), 9))

    for i in range(len(t)):
        # 角速度センサ信号計測値
        ang_vel_vec = [noisy_signal1[i], noisy_signal2[i], noisy_signal3[i]]
        # 加速度センサ計測値(変換元座標での計測値)
        obs_acc = np.array([-0.2, 1.0, 9.81]) + np.array([1, 1, 1]) * base_signal[i]
        # 加速度座標変換の更新
        tar_acc, ang_acc = acc_transformer.update(obs_acc, ang_vel_vec)

        loggs[i, 0:3] = tar_acc
        loggs[i, 3:6] = ang_acc
        loggs[i, 6:9] = obs_acc

    # ---------------------- 結果の可視化
    # 3×3のグラフマトリクスを作成（X軸をリンク）
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)

    # 各軸に名前を割り当て（アクセスしやすくするため）
    ax11, ax12, ax13 = axes[0, 0], axes[0, 1], axes[0, 2]
    ax21, ax22, ax23 = axes[1, 0], axes[1, 1], axes[1, 2]
    ax31, ax32, ax33 = axes[2, 0], axes[2, 1], axes[2, 2]

    # X列（左列）
    # 上段
    ax11.plot(t, loggs[:, 6], label="org_acc_x")
    ax11.plot(t, loggs[:, 0], label="converted_acc_x")
    ax11.set_ylabel("acc. [m/s^2]")
    ax11.set_title("X-axis")
    ax11.legend()
    ax11.grid()
    # 中段
    ax21.plot(t, noisy_signal1, label="ang_vel_x")
    ax21.set_ylabel("ang vel. [rad/s]")
    ax21.legend()
    ax21.grid()
    # 下段
    ax31.plot(t, loggs[:, 3], label="filter output")
    ax31.plot(t, dot_signal, label="true", linestyle="--")
    ax31.set_xlabel("Time [s]")
    ax31.set_ylabel("ang acc. [rad/s^2]")
    ax31.legend()
    ax31.grid()

    # Y列（中央列）
    # 上段
    ax12.plot(t, loggs[:, 7], label="org_acc_y")
    ax12.plot(t, loggs[:, 1], label="converted_acc_y")
    ax12.set_title("Y-axis")
    ax12.legend()
    ax12.grid()
    # 中段
    ax22.plot(t, noisy_signal2, label="ang_vel_y")
    ax22.legend()
    ax22.grid()
    # 下段
    ax32.plot(t, loggs[:, 4], label="filter output")
    ax32.plot(t, 2 * dot_signal, label="true", linestyle="--")
    ax32.set_xlabel("Time [s]")
    ax32.legend()
    ax32.grid()

    # Z列（右列）
    # 上段
    ax13.plot(t, loggs[:, 8], label="org_acc_z")
    ax13.plot(t, loggs[:, 2], label="converted_acc_z")
    ax13.set_title("Z-axis")
    ax13.legend()
    ax13.grid()
    # 中段
    ax23.plot(t, noisy_signal3, label="ang_vel_z")
    ax23.legend()
    ax23.grid()
    # 下段
    ax33.plot(t, loggs[:, 5], label="filter output")
    ax33.plot(t, 3 * dot_signal, label="true", linestyle="--")
    ax33.set_xlabel("Time [s]")
    ax33.legend()
    ax33.grid()

    plt.tight_layout()  # レイアウトを調整
    plt.show()

# ----------------------------------------------------------------------- ノッチフィルタの動作確認テスト
if 1:
    Fs = 100  # サンプリング周波数
    max_time = 10.0  # 信号の長さ（秒）
    t = np.arange(0, max_time, 1 / Fs)

    # 元信号（1Hzの正弦波）
    base_signal = np.sin(2 * np.pi * 1 * t) + 5

    # 25Hzのノイズ成分を重畳
    noise_freq = 25  # Hz
    noise_signal = 2.0 * np.sin(2 * np.pi * noise_freq * t)

    # ノイズを含んだ信号
    noisy_signal = base_signal + noise_signal

    # ノッチフィルタの初期化
    # omega_n: ノッチの中心周波数 [rad/s]
    # zeta: ノッチの幅（減衰係数）
    # d: ノッチの深さ（0に近いほど深い）
    omega_n = 2 * np.pi * noise_freq  # 25Hzに対応
    zeta = 0.1  # ノッチの幅
    d = 0.01  # ノッチの深さ（0.01 = ほぼ完全に除去）

    print("ノッチフィルタ: 25Hzを除去")
    notch_filter = zf.Z_Filter_Notch(omega_n=omega_n, zeta=zeta, d=d, sampling_freq=Fs)

    # フィルタリング処理（時系列処理の模擬）
    filtered_signal = np.zeros(len(t))

    for i in range(len(t)):
        if i < 10:
            # リセットcommandでフィルタを現在値に初期化
            notch_filter.reset(noisy_signal[i])
        filtered_signal[i] = notch_filter.update(noisy_signal[i])

    # ---------------------- 結果の可視化
    # 3段のグラフを作成（X軸をリンク）
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    # 上段：ノイズを含んだ信号と元信号
    ax1.plot(t, noisy_signal, label="Noisy Signal (Base + 25Hz)", alpha=0.7)
    ax1.plot(t, base_signal, label="Base Signal (1Hz)", linestyle="--", linewidth=2)
    ax1.set_ylabel("Amplitude")
    ax1.set_title("Input Signal: Base Signal with 25Hz Noise")
    ax1.legend()
    ax1.grid()

    # 中段：フィルタ後の信号と元信号の比較
    ax2.plot(t, filtered_signal, label="Filtered Signal (Notch Filter)", linewidth=2)
    ax2.plot(t, base_signal, label="Base Signal (1Hz)", linestyle="--", alpha=0.7)
    ax2.set_ylabel("Amplitude")
    ax2.set_title("Output Signal: After Notch Filter (25Hz Removed)")
    ax2.legend()
    ax2.grid()

    # 下段：除去されたノイズ成分（理論値との比較）
    removed_noise = noisy_signal - filtered_signal
    ax3.plot(t, removed_noise, label="Removed Component", alpha=0.7)
    ax3.plot(t, noise_signal, label="Original 25Hz Noise", linestyle="--", alpha=0.7)
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Amplitude")
    ax3.set_title("Removed Noise Component vs Original 25Hz Noise")
    ax3.legend()
    ax3.grid()

    plt.tight_layout()
    plt.show()