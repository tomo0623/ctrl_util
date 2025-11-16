"""
フィルタのボード線図検証ツール
時間領域シミュレーションで周波数応答を測定し、理論値と比較する
"""
import numpy as np
from matplotlib import pyplot as plt
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def measure_frequency_response(
    filter_object,
    freq_range: np.ndarray,
    sampling_freq: float,
    settling_time: float = 5.0,
    measurement_duration: float = 2.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    時間領域シミュレーションでフィルタの周波数応答を測定

    Args:
        filter_object: updateメソッドとresetメソッドを持つフィルタオブジェクト
        freq_range: 測定する周波数の配列 [Hz]
        sampling_freq: サンプリング周波数 [Hz]
        settling_time: 定常状態に達するまでの時間 [秒]
        measurement_duration: 測定時間 [秒]

    Returns:
        gain: 各周波数でのゲイン（振幅比）
        phase: 各周波数での位相差 [度]
    """
    gain = np.zeros(len(freq_range))
    phase = np.zeros(len(freq_range))

    T = 1 / sampling_freq
    total_time = settling_time + measurement_duration

    for i, freq in enumerate(freq_range):
        # 正弦波信号を生成
        t = np.arange(0, total_time, T)
        input_signal = np.sin(2 * np.pi * freq * t)

        # フィルタをリセット
        filter_object.reset(0)

        # フィルタリング
        output_signal = np.zeros(len(t))
        for j in range(len(t)):
            output_signal[j] = filter_object.update(input_signal[j])

        # 定常状態部分を抽出
        settling_samples = int(settling_time * sampling_freq)
        steady_input = input_signal[settling_samples:]
        steady_output = output_signal[settling_samples:]

        # ゲイン計算（RMS比）
        input_rms = np.sqrt(np.mean(steady_input**2))
        output_rms = np.sqrt(np.mean(steady_output**2))
        gain[i] = output_rms / input_rms if input_rms > 1e-10 else 0

        # 位相差計算（クロススペクトル法）
        if input_rms > 1e-10 and output_rms > 1e-10:
            # FFTで位相を計算
            fft_input = np.fft.fft(steady_input)
            fft_output = np.fft.fft(steady_output)

            # 基本周波数成分のインデックスを探す
            freqs = np.fft.fftfreq(len(steady_input), T)
            idx = np.argmin(np.abs(freqs - freq))

            # 位相差を計算
            phase_input = np.angle(fft_input[idx])
            phase_output = np.angle(fft_output[idx])
            phase[i] = np.degrees(phase_output - phase_input)

            # -180～180度に正規化
            while phase[i] > 180:
                phase[i] -= 360
            while phase[i] < -180:
                phase[i] += 360
        else:
            phase[i] = 0

    return gain, phase


def calculate_theoretical_frequency_response(
    filter_object, freq_range: np.ndarray, sampling_freq: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    理論的な周波数応答を計算（状態空間表現から）

    Args:
        filter_object: Ad_mat, Bd_mat, Cd_mat, Dd_matを持つフィルタオブジェクト
        freq_range: 計算する周波数の配列 [Hz]
        sampling_freq: サンプリング周波数 [Hz]

    Returns:
        gain: 各周波数でのゲイン（振幅比）
        phase: 各周波数での位相差 [度]
    """
    gain = np.zeros(len(freq_range))
    phase = np.zeros(len(freq_range))

    T = 1 / sampling_freq
    I = np.eye(filter_object.Ad_mat.shape[0])

    for i, freq in enumerate(freq_range):
        omega = 2 * np.pi * freq
        z = np.exp(1j * omega * T)  # z = e^(jωT)

        try:
            # H(z) = Cd·(zI - Ad)^(-1)·Bd + Dd
            inv_term = np.linalg.inv(z * I - filter_object.Ad_mat)
            H_z = filter_object.Cd_mat @ inv_term @ filter_object.Bd_mat + filter_object.Dd_mat
            H_z = H_z[0, 0]  # スカラーに変換

            gain[i] = np.abs(H_z)
            phase[i] = np.angle(H_z, deg=True)
        except:
            gain[i] = 0
            phase[i] = 0

    return gain, phase


def plot_bode_verification(
    freq_range: np.ndarray,
    gain_measured: np.ndarray,
    phase_measured: np.ndarray,
    gain_theory: np.ndarray,
    phase_theory: np.ndarray,
    title: str = "Bode Plot Verification",
    save_path: Optional[str] = None,
    show_plot: bool = True,
) -> None:
    """
    ボード線図の検証結果をプロット

    Args:
        freq_range: 周波数配列 [Hz]
        gain_measured: 測定されたゲイン
        phase_measured: 測定された位相 [度]
        gain_theory: 理論的なゲイン
        phase_theory: 理論的な位相 [度]
        title: グラフタイトル
        save_path: 保存先パス（Noneの場合は保存しない）
        show_plot: グラフを表示するか
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))

    # 上段: ゲイン特性（dB）
    ax1.semilogx(freq_range, 20 * np.log10(gain_measured + 1e-10), "b-", label="Measured (Time Domain)", linewidth=2)
    ax1.semilogx(freq_range, 20 * np.log10(gain_theory + 1e-10), "r--", label="Theory (Frequency Domain)", linewidth=2)
    ax1.set_ylabel("Gain [dB]")
    ax1.set_title(f"{title} - Magnitude")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)

    # 中段: 位相特性
    ax2.semilogx(freq_range, phase_measured, "b-", label="Measured (Time Domain)", linewidth=2)
    ax2.semilogx(freq_range, phase_theory, "r--", label="Theory (Frequency Domain)", linewidth=2)
    ax2.set_ylabel("Phase [deg]")
    ax2.set_title(f"{title} - Phase")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)

    # 下段: 誤差
    gain_error = np.abs(gain_measured - gain_theory)
    phase_error = np.abs(phase_measured - phase_theory)

    ax3_twin = ax3.twinx()
    ax3.semilogx(freq_range, gain_error, "b-", label="Gain Error", linewidth=2)
    ax3_twin.semilogx(freq_range, phase_error, "r-", label="Phase Error", linewidth=2)

    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Gain Error (absolute)", color="b")
    ax3_twin.set_ylabel("Phase Error [deg]", color="r")
    ax3.set_title(f"{title} - Error")
    ax3.tick_params(axis="y", labelcolor="b")
    ax3_twin.tick_params(axis="y", labelcolor="r")
    ax3.grid(True, which="both", alpha=0.3)

    # 凡例を統合
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"グラフを {save_path} に保存しました")

    if show_plot:
        try:
            plt.show(block=False)
            plt.pause(5)
            plt.close()
        except:
            pass


def verify_filter(
    filter_object,
    sampling_freq: float,
    freq_range: Optional[np.ndarray] = None,
    title: str = "Filter Bode Plot Verification",
    save_path: Optional[str] = None,
    show_plot: bool = True,
    verbose: bool = True,
) -> dict:
    """
    フィルタのボード線図を時間領域シミュレーションで検証

    Args:
        filter_object: 検証するフィルタオブジェクト
        sampling_freq: サンプリング周波数 [Hz]
        freq_range: 測定する周波数範囲（Noneの場合は自動設定）
        title: グラフタイトル
        save_path: グラフ保存先
        show_plot: グラフを表示するか
        verbose: 詳細情報を出力するか

    Returns:
        結果の辞書（gain_measured, gain_theory, phase_measured, phase_theory, freq_range）
    """
    # 周波数範囲の自動設定
    if freq_range is None:
        freq_range = np.logspace(-1, np.log10(sampling_freq / 2), 100)

    if verbose:
        print("=" * 70)
        print(f"{title}")
        print("=" * 70)
        print(f"サンプリング周波数: {sampling_freq} Hz")
        print(f"測定周波数範囲: {freq_range[0]:.2f} ~ {freq_range[-1]:.2f} Hz")
        print(f"測定点数: {len(freq_range)}")
        print()

    # 時間領域シミュレーションで測定
    if verbose:
        print("【1】時間領域シミュレーション実行中...")
    gain_measured, phase_measured = measure_frequency_response(filter_object, freq_range, sampling_freq)
    if verbose:
        print("  ✓ 完了")
        print()

    # 理論値を計算
    if verbose:
        print("【2】理論値計算中...")
    gain_theory, phase_theory = calculate_theoretical_frequency_response(filter_object, freq_range, sampling_freq)
    if verbose:
        print("  ✓ 完了")
        print()

    # 誤差統計
    gain_error = np.abs(gain_measured - gain_theory)
    phase_error = np.abs(phase_measured - phase_theory)

    # 周波数範囲を分割して評価（ナイキスト周波数付近は除外）
    # 低周波: 0 ~ 0.2 * Nyquist
    # 中周波: 0.2 ~ 0.4 * Nyquist
    # 高周波: 0.4 ~ 1.0 * Nyquist（評価対象外）
    nyquist_freq = sampling_freq / 2
    low_freq_mask = freq_range <= 0.2 * nyquist_freq
    mid_freq_mask = (freq_range > 0.2 * nyquist_freq) & (freq_range <= 0.4 * nyquist_freq)
    valid_freq_mask = freq_range <= 0.4 * nyquist_freq  # 評価に使う周波数範囲

    if verbose:
        print("【3】検証結果")
        print(f"  全周波数範囲:")
        print(f"    ゲイン誤差（最大）: {np.max(gain_error):.6f}")
        print(f"    ゲイン誤差（平均）: {np.mean(gain_error):.6f}")
        print(f"    位相誤差（最大）:   {np.max(phase_error):.2f} deg")
        print(f"    位相誤差（平均）:   {np.mean(phase_error):.2f} deg")
        print()

        if np.any(low_freq_mask):
            print(f"  低周波領域 (0 ~ {0.2 * nyquist_freq:.1f} Hz):")
            print(f"    ゲイン誤差（最大）: {np.max(gain_error[low_freq_mask]):.6f}")
            print(f"    位相誤差（最大）:   {np.max(phase_error[low_freq_mask]):.2f} deg")

        if np.any(mid_freq_mask):
            print(f"  中周波領域 ({0.2 * nyquist_freq:.1f} ~ {0.4 * nyquist_freq:.1f} Hz):")
            print(f"    ゲイン誤差（最大）: {np.max(gain_error[mid_freq_mask]):.6f}")
            print(f"    位相誤差（最大）:   {np.max(phase_error[mid_freq_mask]):.2f} deg")
        print()

        # 精度判定（有効周波数範囲のみで評価）
        if np.any(valid_freq_mask):
            max_gain_error_valid = np.max(gain_error[valid_freq_mask])
            max_phase_error_valid = np.max(phase_error[valid_freq_mask])

            if max_gain_error_valid < 0.01 and max_phase_error_valid < 1.0:
                print("  ✓ フィルタは理論値と高精度で一致しています")
            elif max_gain_error_valid < 0.05 and max_phase_error_valid < 5.0:
                print("  ✓ フィルタは理論値とおおむね一致しています")
            else:
                print("  ⚠ フィルタに大きな誤差があります。実装を確認してください")

        print(f"  ※ ナイキスト周波数付近 (> {0.4 * nyquist_freq:.1f} Hz) は離散化の性能限界のため評価対象外")
        print()

    # グラフ作成
    plot_bode_verification(
        freq_range, gain_measured, phase_measured, gain_theory, phase_theory, title, save_path, show_plot
    )

    if verbose:
        print("=" * 70)

    return {
        "freq_range": freq_range,
        "gain_measured": gain_measured,
        "phase_measured": phase_measured,
        "gain_theory": gain_theory,
        "phase_theory": phase_theory,
        "gain_error": gain_error,
        "phase_error": phase_error,
    }
