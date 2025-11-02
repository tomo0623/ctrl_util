# Z領域における離散フィルタ計算を行うClass
# ctrl_util/z_filter/z_filter.py

import numpy as np
from typing import Union, List


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
    print("可制御正準系の行列 A_c:")
    print(A_c)
    print("可制御正準系の行列 B_c:")
    print(B_c)
    print("可制御正準系の行列 C_c:")
    print(C_c)
    print("直達項 D_c:")
    print(D_c)

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
    print("離散化行列 A_d:")
    print(A_d)
    print("離散化行列 B_d:")
    print(B_d)
    print("離散化行列 C_d:")
    print(C_d)
    print("離散化行列 D_d:")
    print(D_d)

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

        # 連続時間の可制御正準系の行列式を計算(
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
