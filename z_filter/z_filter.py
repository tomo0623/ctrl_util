# Z領域における離散フィルタ計算を行うClass
# ctrl_util/z_filter/z_filter.py

import numpy as np
import logging
from typing import Union, List

# ロガーの設定 --- ライブラリはデフォルトで静かに振る舞う
# ユーザー側で明示的にハンドラ/レベルを設定するか、以下のヘルパを呼ぶ
logger = logging.getLogger(__name__)
# ライブラリはハンドラを追加せず、NullHandlerを設定しておくのが推奨パターン
if not logger.handlers:
    logger.addHandler(logging.NullHandler())


def enable_verbose_logging(level: int = logging.INFO):
    """コンソール出力を有効にするヘルパ

    Args:
        level: 任意の logging レベル (デフォルト: logging.INFO)
    """
    logger.handlers.clear()
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(level)


def disable_verbose_logging():
    """コンソール出力を無効にしてライブラリを静かにするヘルパ"""
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())


def set_log_level(level: int):
    """ロガーのログレベルを設定するユーティリティ

    Args:
        level: logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR など
    """
    logger.setLevel(level)


def set_canonical_form(a_correct: Union[List[float], np.ndarray], c_correct: Union[List[float], np.ndarray]) -> tuple:
    """
    連続時間系の伝達関数を可制御正準系に変換する

    Args:
        a_correct: 正規の順番の分母多項式の係数 (a0, a1, a2, ...)
        c_correct: 正規の順番の分子多項式の係数 (c0, c1, c2, ...)

    Returns:
        可制御正準系の行列 A_c と B_c と C_c と D_c のタプル
    """
    a_correct = np.array(a_correct, dtype=float)
    c_correct = np.array(c_correct, dtype=float)

    if a_correct[-1] == 0:
        raise ValueError("分母多項式の先頭の係数がゼロのためエラー, プロパーな伝達関数を設定すること")

    # 正規化（a[-1] = 1にする）
    a_correct = a_correct / a_correct[-1]
    c_correct = c_correct / a_correct[-1]

    # 分子と分母の次数が同じ場合、多項式の長除法で直達項を分離
    D_c = 0.0
    if len(c_correct) == len(a_correct):
        # 直達項 D = c[-1] / a[-1] (既に正規化済みなので a[-1] = 1)
        D_c = c_correct[-1]
        # 分子から直達項を除去: c'(s) = c(s) - D * a(s)
        c_correct = c_correct - D_c * a_correct

    # 可制御正準系の行列を計算
    n = len(a_correct) - 1
    m = len(c_correct)
    A_c = np.zeros((n, n))
    B_c = np.zeros((n, 1))
    B_c[n - 1, 0] = 1.0
    C_c = np.zeros((1, n))

    for i in range(n):
        A_c[-1, i] = -a_correct[i]

        if i != 0:
            A_c[i - 1, i] = 1.0

        if i < m:
            C_c[0, i] = c_correct[i]

    # 可制御正準行列の確認
    logger.debug("可制御正準系の行列 A_c:")
    logger.debug(f"{A_c}")
    logger.debug("可制御正準系の行列 B_c:")
    logger.debug(f"{B_c}")
    logger.debug("可制御正準系の行列 C_c:")
    logger.debug(f"{C_c}")
    logger.debug("直達項 D_c:")
    logger.debug(f"{D_c}")

    return A_c, B_c, C_c, D_c


def c2d(
    Ac_mat: np.array,
    Bc_mat: np.array,
    Cc_mat: np.array,
    Dc_mat: np.array,
    Ts: float,
    method: str = "tustin",
) -> tuple:
    """
    連続時間系の可制御正準系を離散時間系に変換する
    Args:
        Ac_mat: 連続時間系の行列 A_c
        Bc_mat: 連続時間系の行列 B_c
        Cc_mat: 連続時間系の行列 C_c
        Dc_mat: 連続時間系の行列 D_c
        Ts: サンプリング周期 [秒]
        method: 離散化手法 tustin / euler
    Returns:
        A_d, B_d, C_d, D_d: 離散時間系の行列 A_d, B_d, C_d, D_d のタプル
    """

    # 離散化された行列を初期化
    A_d = np.zeros_like(Ac_mat)
    B_d = np.zeros_like(Bc_mat)
    C_d = np.zeros_like(Cc_mat)
    D_d = np.zeros_like(Dc_mat)

    # 離散化処理
    if method == "tustin":
        # タスティン変換による離散化
        phi_mat = np.linalg.inv(np.eye(Ac_mat.shape[0]) - (Ts / 2) * Ac_mat)
        A_d = (np.eye(Ac_mat.shape[0]) + (Ts / 2) * Ac_mat) @ phi_mat
        B_d = phi_mat @ Bc_mat * Ts
        C_d = Cc_mat @ phi_mat
        D_d = Cc_mat @ phi_mat @ Bc_mat * (Ts / 2) + Dc_mat
    elif method == "euler":
        # オイラー法で離散化
        A_d = Ac_mat * Ts + np.eye(Ac_mat.shape[0])
        B_d = Bc_mat * Ts
        C_d = Cc_mat
        D_d = Dc_mat

    # 離散化行列の確認
    logger.debug("離散化行列 A_d:")
    logger.debug(f"{A_d}")
    logger.debug("離散化行列 B_d:")
    logger.debug(f"{B_d}")
    logger.debug("離散化行列 C_d:")
    logger.debug(f"{C_d}")
    logger.debug("離散化行列 D_d:")
    logger.debug(f"{D_d}")

    return A_d, B_d, C_d, D_d


class Z_Filter:
    """Z変換を用いた離散時間フィルタクラス"""

    def __init__(
        self,
        denominator: Union[List[float], np.ndarray],
        numerator: Union[List[float], np.ndarray],
        sampling_freq: float = 1.0,
        is_prewarping: bool = False,
        prewarping_freq: float = None,
    ):
        """
        Z変換フィルタを初期化
        - 連続時間系の伝達関数の分子/分母を入力として受け取る
        - サンプリング周期を指定して離散化
        (プリワーピング処理は初期化時のサンプリング周期を補正して対応する)

        Args:
            denominator: 連続時間系の分母の係数 (an, an-1, ..., a2, a1, a0)
            numerator: 連続時間系の分子の係数 (cn, cn-1, ..., c2, c1, c0)
            sampling_freq: サンプリング周波数 [Hz]
        """

        self.a = np.array(denominator, dtype=float)
        self.c = np.array(numerator, dtype=float)
        self.fs = sampling_freq

        # プロパーな伝達関数を確認（分子の次数 <= 分母の次数）
        if len(self.a) < len(self.c):
            raise ValueError("プロパーな伝達関数ではないのでエラー（分子の次数 > 分母の次数）")

        # 正規化（a[0] = 1にする）
        if self.a[0] != 1.0:
            self.c = self.c / self.a[0]
            self.a = self.a / self.a[0]

        # 状態ベクトルの初期化（2次元配列として初期化）
        self.xnew_vec = np.zeros((len(self.a) - 1, 1), dtype=float)
        self.xold_vec = np.zeros((len(self.a) - 1, 1), dtype=float)

        # 連続時間の可制御正準系の行列式を計算
        # NOTE: 係数行列は逆順で渡すことで正規の多項式の形に対応
        self.Ac_mat, self.Bc_mat, self.Cc_mat, self.Dc_mat = set_canonical_form(self.a[::-1], self.c[::-1])

        # タスティン変換に基づく離散化
        if is_prewarping and prewarping_freq is not None:
            # プリワーピング処理
            T = 1 / self.fs
            omega_0 = 2 * np.pi * prewarping_freq
            omega_c = omega_0 / np.tan(omega_0 * T / 2)
            self.fs_modified = omega_c / (2 * np.pi)

            # 離散化処理
            self.Ad_mat, self.Bd_mat, self.Cd_mat, self.Dd_mat = c2d(
                self.Ac_mat, self.Bc_mat, self.Cc_mat, self.Dc_mat, 1 / self.fs_modified
            )
        else:
            # プリワーピング処理なしで離散化
            self.Ad_mat, self.Bd_mat, self.Cd_mat, self.Dd_mat = c2d(
                self.Ac_mat, self.Bc_mat, self.Cc_mat, self.Dc_mat, 1 / self.fs
            )

        # 内部リセット処理実行
        self.reset(np.zeros((len(self.a) - 1, 1), dtype=float))

    def reset(self, xini_vec) -> None:
        """
        #     内部リセット（リアルタイム処理対応）
        #     Args:
        #         xini_vec: フィルタ初期値
        #     Returns:
        #         なし
        """
        self.xnew_vec = xini_vec

    def update(self, u) -> float:
        """
        #     フィルタ更新関数（リアルタイム処理対応）
        #     Args:
        #         u: 入力（スカラー値）
        #     Returns:
        #         y: 出力（スカラー値）
        #
        """

        # 離散状態/出力方程式更新
        self.xold_vec = self.xnew_vec.copy()
        self.xnew_vec = self.Ad_mat @ self.xold_vec + self.Bd_mat * u
        y = self.Cd_mat @ self.xold_vec + self.Dd_mat * u
        return y[0, 0]

# 一次LPFフィルタクラス 
class Z_Filter_LPF(Z_Filter):
    """
    一次LPFフィルタクラス
    """

    def __init__(self, tau_lpf: float, sampling_freq: float = 1.0):
        """
        一次LPFフィルタの初期化
        Args:
            tau_lpf: LPFの時定数 [秒]
            sampling_freq: サンプリング周波数 [Hz]
        """
        c = np.array([1])
        a = np.array([tau_lpf, 1])
        super().__init__(a, c, sampling_freq=sampling_freq)

# 一次微分フィルタクラス
class Z_Filter_ADF(Z_Filter):
    """
    一次微分フィルタクラス
    """

    def __init__(self, tau_derivative: float, sampling_freq: float = 1.0):
        """
        一次微分フィルタの初期化
        Args:
            tau_derivative: 微分の時定数 [秒]
            sampling_freq: サンプリング周波数 [Hz]
        """
        c = np.array([1, 0])
        a = np.array([tau_derivative, 1])
        super().__init__(a, c, sampling_freq=sampling_freq)

# 一次微分フィルタとLPFを組み合わせたフィルタクラス
class Z_Filter_ADF_with_LPF(Z_Filter):
    """
    一次微分フィルタとLPFを組み合わせたフィルタクラス
    """

    def __init__(self, tau_lpf: float, tau_adf: float, sampling_freq: float = 1.0):
        """
        一次微分フィルタとLPFを組み合わせたフィルタの初期化
        Args:
            tau_lpf: LPFの時定数 [秒]
            tau_adf: 微分の時定数 [秒]
            sampling_freq: サンプリング周波数 [Hz]
        """
        c = np.array([1 / (tau_lpf * tau_adf), 0])
        a = np.array([1, (tau_lpf + tau_adf) / (tau_lpf * tau_adf), 1 / (tau_lpf * tau_adf)])
        super().__init__(a, c, sampling_freq=sampling_freq)

# N次のバターワースフィルタクラス
class Z_Filter_Butterworth(Z_Filter):
    """
    N次のバターワースフィルタクラス
    """

    def __init__(self, order: int, cutoff_freq: float, sampling_freq: float = 1.0):
        """
        N次のバターワースフィルタの初期化
        Args:
            order: フィルタの次数
            cutoff_freq: カットオフ周波数 [Hz]
            sampling_freq: サンプリング周波数 [Hz]
        """
        # バターワースフィルタの連続時間系係数を計算
        a_cont, c_cont = self._calculate_butterworth_coefficients(order, cutoff_freq)
        super().__init__(a_cont, c_cont, sampling_freq=sampling_freq)

    def _calculate_butterworth_coefficients(self, order: int, cutoff_freq: float) -> tuple:
        """
        バターワースフィルタの連続時間系係数を計算
        
        Args:
            order: フィルタの次数
            cutoff_freq: カットオフ周波数 [Hz]
            
        Returns:
            (denominator, numerator): 分母係数と分子係数のタプル
        """
        # 正規化カットオフ周波数 (rad/s)
        omega_c = 2 * np.pi * cutoff_freq
        
        # バターワース極の計算
        # s_k = omega_c * exp(j * (pi/2 + (2k-1)*pi/(2*n))) for k = 1, 2, ..., n
        # または s_k = omega_c * exp(j * (pi/2 + (2k+1)*pi/(2*n))) for k = 0, 1, ..., n-1
        all_poles = []
        for k in range(order):
            angle = np.pi / 2 + (2 * k + 1) * np.pi / (2 * order)
            pole = omega_c * np.exp(1j * angle)
            all_poles.append(pole)
        
        logger.debug(f"全ての極 (order={order}):")
        for i, pole in enumerate(all_poles):
            logger.debug(f"  極{i+1}: {pole:.6f} (real={pole.real:.6f})")
        
        # 左半平面の極のみを選択（安定な極）
        stable_poles = [p for p in all_poles if p.real < -1e-10]  # より厳密な判定
        
        logger.debug("左半平面の極 (安定極):")
        for i, pole in enumerate(stable_poles):
            logger.debug(f"  安定極{i+1}: {pole:.6f}")
        
        if len(stable_poles) == 0:
            raise ValueError(f"安定な極が見つかりません (order={order})")
        
        # 複素共役極を実係数多項式として処理
        denominator = [1.0]  # 最高次の係数
        processed_poles = set()  # 処理済みの極のインデックス
        
        for i, pole in enumerate(stable_poles):
            if i in processed_poles:
                continue
                
            if abs(pole.imag) < 1e-10:  # 実極の場合
                # (s - pole) との積
                new_denominator = [0.0] * (len(denominator) + 1)
                for j in range(len(denominator)):
                    new_denominator[j] += denominator[j]
                    new_denominator[j + 1] += -pole.real * denominator[j]
                denominator = new_denominator
                processed_poles.add(i)
                
            else:  # 複素極の場合（共役ペア）
                # 共役ペアを探す
                conjugate_pole = np.conj(pole)
                conjugate_index = None
                
                for j, other_pole in enumerate(stable_poles):
                    if j != i and j not in processed_poles:
                        if abs(other_pole - conjugate_pole) < 1e-10:
                            conjugate_index = j
                            break
                
                if conjugate_index is not None:
                    # (s - pole)(s - conj(pole)) = s² - 2*Re(pole)*s + |pole|²
                    real_part = pole.real
                    magnitude_squared = pole.real**2 + pole.imag**2
                    
                    # s² - 2*real_part*s + magnitude_squared との積
                    new_denominator = [0.0] * (len(denominator) + 2)
                    for j in range(len(denominator)):
                        new_denominator[j] += denominator[j]  # s²項
                        new_denominator[j + 1] += -2 * real_part * denominator[j]  # s項
                        new_denominator[j + 2] += magnitude_squared * denominator[j]  # 定数項
                    
                    denominator = new_denominator
                    processed_poles.add(i)
                    processed_poles.add(conjugate_index)
                else:
                    # 共役ペアが見つからない場合（エラー）
                    raise ValueError(f"複素極 {pole} の共役ペアが見つかりません")
        
        # 実数部のみを取得（虚数部は数値誤差）
        denominator = [coeff.real if np.isreal(coeff) else coeff for coeff in denominator]
        
        # DC gain = 1 となるように分子を調整
        # H(0) = numerator(0) / denominator(0) = 1
        # 分子は定数項のみ（ローパスフィルタ）
        dc_gain = denominator[-1]  # s=0での分母の値
        numerator = [dc_gain]  # 分子は定数項のみ
        
        logger.debug("計算された係数:")
        logger.debug(f"  分母: {denominator}")
        logger.debug(f"  分子: {numerator}")
        
        return denominator, numerator