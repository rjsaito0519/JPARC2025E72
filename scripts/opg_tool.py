import numpy as np
import matplotlib.pyplot as plt

def check_condition(values, threshold = 0.3):
    """
    外れ値がないか確認する
    外れ値があったらTrueを返す
    """
    values = np.array(values)
    if len(values) < 2:
        return False
    
    mean = np.sum(values)/len(values)
    for it in values:
        if np.abs(it - mean) > threshold:
            return True
    return False

def weighted_mean_and_error(values, weights):
    """
    重み付き平均とそのエラーを計算する関数。
    np.nanが含まれる場合は、そのデータを除外して計算する。

    Parameters:
        values (list or np.ndarray): データ点
        weights (list or np.ndarray): 重み

    Returns:
        tuple: (重み付き平均, エラー)
    """
    # 配列に変換
    values = np.array(values)
    weights = np.array(weights)
    
    # np.nanを除外するマスクを作成
    mask = ~np.isnan(values) & ~np.isnan(weights)
    values = values[mask]
    weights = weights[mask]

    # マスク後にデータがなくなった場合のエラーチェック
    if len(values) == 0 or len(weights) == 0:
        raise ValueError("All values or weights are NaN.")

    # 重み付き平均の計算
    weighted_mean = np.sum(weights * values) / np.sum(weights)

    # 標準誤差の計算
    weighted_error = np.sqrt(1 / np.sum(weights))

    return weighted_mean, weighted_error

# -- data summarize -----
def data_summarize(tree, ch, mppc_map):
    data = []
    prev_mppc = 0
    val, err = 0, np.inf
    for i in range(len(tree["run_num"])):
        if tree["run_num"][i] in mppc_map.keys() and tree["ch"][i] == ch:
            if (mppc_map[tree["run_num"][i]] - prev_mppc) > 1:
                data.append([ np.nan, np.nan ])
            if prev_mppc != mppc_map[tree["run_num"][i]]:
                if tree["result_err"][i][3] < err:
                    val = tree["result_val"][i][3]
                    err = tree["result_err"][i][3]
                data.append([ val, err ])
                prev_mppc = mppc_map[tree["run_num"][i]]
                val, err = 0, np.inf
            else:
                if tree["result_err"][i][3] < err:
                    val = tree["result_val"][i][3]
                    err = tree["result_err"][i][3]
    return np.array(data)

# -- data summarize -----
def data_summarize_for_jparc(tree, mppc_map, is_up):
    
    data = []
    for mppc_ch in range(1, 17):
        prev_mppc = 0
        val, err = 0, np.inf
        flag = False
        for i in range(len(tree["run_num"])):
            if tree["run_num"][i] in mppc_map[mppc_ch]:
                if is_up:
                    flag = tree["ch"][i] < 4
                else:
                    flag = tree["ch"][i] >= 4
                
                if flag and tree["result_err"][i][3] < err:
                    val = tree["result_val"][i][3]
                    err = tree["result_err"][i][3]
        if val == 0:
            data.append([np.nan, np.nan])
        else:
            data.append([ val, err ])
    return np.array(data)