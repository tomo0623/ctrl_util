# 加速度情報を剛体内の座標変換する
# ctrl_util/accel_transform/accel_trans.py

import numpy as np
import logging
from typing import Union, List
from z_filter import z_filter as zf

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
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
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


class AccTransform:
    """剛体内の加速度信号の座標変換クラス"""

    def __init__(
        self, obs2tar_vec: Union[List[float], np.ndarray], sampling_freq: float, tau_lpf: float, tau_adf: float
    ):
        """
        剛体内の加速度信号を座標変換の初期化
        - 座標変換用のアーム長(変換元→変換先のベクトル)を設定する
        - 微分計算用の各種フィルタ設定をする(サンプリング周期, フィルタ定数)
        - 右手座標系で規定された空間であることを前提とする
        - 角加速度算出の微分は近似微分計算でノイズ対策してるためフィルタ時定数は適切に設定すること

        Args:
            obs2tar_vec: 変換元から変換先へのベクトル[x, y, z]
            sampling_freq: サンプリング周波数 [Hz]
            tau_lpf: LPF時定数(近似微分計算用) [秒]
            tau_adf: ADF時定数(近似微分計算用) [秒]
        """
        self.obs2tar_vec = np.array(obs2tar_vec)
        self.fs = sampling_freq
        if tau_lpf < 1 / sampling_freq * 2 or tau_adf < 1 / sampling_freq * 2:
            raise ValueError(f"フィルタ設定がサンプリング定義に反している。時定数期待値>={1 / sampling_freq * 2}")
        self.tau_lpf = tau_lpf
        self.tau_adf = tau_adf

        # 微分計算用フィルタの初期化
        self.ang_acc_filter = [None] * 3
        for i in range(3):
            self.ang_acc_filter[i] = zf.Z_Filter_ADF_with_LPF(self.tau_lpf, self.tau_adf, sampling_freq=self.fs)

        # 設定値の確認
        logger.debug("変換ベクトル(アーム長) obs2tar_vec:")
        logger.debug(f"{self.obs2tar_vec}")
        logger.debug("角加速度計算フィルタ設定")
        logger.debug(f"Fs:{self.fs}, lpf:{self.tau_lpf}, adf:{self.tau_adf}")

    def update(self, obs_acc: Union[List[float], np.ndarray], ang_vel_vec: Union[List[float], np.ndarray]) -> tuple:
        """
        加速度座標変換の更新関数（リアルタイム処理対応）
        Args:
            obs_acc: 変換元加速度ベクトル[m/s^2], [x, y, z]
            ang_vel_vec: 剛体の角速度ベクトル[rad/s], [x, y, z]
        Returns:
            tar_acc: 変換先加速度ベクトル[m/s^2], [x, y, z]
            ang_acc: 角加速度ベクトル(参考値)[rad/s^2], [x, y, z]
        """
        obs_acc = np.array(obs_acc)
        ang_vel_vec = np.array(ang_vel_vec)

        # 角加速度の近似微分計算
        ang_acc = np.zeros(3)
        for i in range(3):
            ang_acc[i] = self.ang_acc_filter[i].update(ang_vel_vec[i])

        # 座標変換計算
        # NOTE: αt = αo + ω_dot×r + ω×(ω×r)
        # オイラー加速度（角加速度による接線加速度）
        cross1 = np.cross(ang_acc, self.obs2tar_vec)
        # 向心加速度（角速度による遠心加速度）, ベクトル3重積
        # NOTE: ω×(ω×r) = (ω・r)ω - (ω・ω)r, 今回は内積に分解せず外積のまま計算する
        cross2 = np.cross(ang_vel_vec, np.cross(ang_vel_vec, self.obs2tar_vec))

        # 観測点加速度に剛体回転の効果を加算して目標点加速度を計算
        tar_acc = obs_acc + cross1 + cross2

        return tar_acc, ang_acc
