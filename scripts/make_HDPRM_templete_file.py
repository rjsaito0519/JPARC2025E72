import os
import conf

# --- ユーザーが設定する項目 ---
hdprm_dir = f"{conf.param_dir}/HDPRM"
OUTPUT_FILE = f"{hdprm_dir}/HodoParam_1"

# 生成したいCIdと、それぞれが持つSegIdの最大値+1
# (例: CId=1はSegId 0-62 (計63個))
DETECTOR_CONFIG = {
     1: 63, # CId 1 BHT
     3: 15, # CId 3 BH2
     4: 4,  # CId 4 BAC
     5: 34, # CId 5 HTOF
     6: 8,  # CId 6 KVC
     7: 1,  # CId 7 T1
     8: 8,  # CId 8 CVC
     9: 1,  # CId 9 SAC3
    10: 5,  # CId 10 SFV
    11: 8   # CId 11 COBO
}

# PlIdのリスト (今回は [0] で固定)
PLID_LIST = [0]

# ★ CIdごとに(AorT, UorD)の組み合わせリストを定義します。
AT_UD_COMBINATIONS_BY_CID = {
    "default": [(1, 0), (1, 1), (0, 0), (0, 1)],
    
    # CId=3: AorT(0,1,2)
     3: [(1, 0), (1, 1), (1, 2), (0, 0), (0, 1)],

    10: [(1, 0), (1, 1)],

    11: [(1, 0)],
}

def get_initial_params(CId, PlId, SegId, AorT, UorD):
    """
    パラメータIDに基づいてp0, p1の初期値を返す関数。
    ここをカスタマイズしてください。
    """
    # ★ CId=3 (BH2) の AorT=2 の場合の初期値
    if CId == 3 and UorD == 2:
        p0 = 0.0  # (適当な初期値)
        p1 = 1.0  # (適当な初期値)
        return p0, p1

    if AorT == 1:
        # AorT=1 (TDC) の場合のデフォルト値
        p0 = 1.0
        p1 = -9.765625e-04
    
    elif AorT == 0:
        # AorT=0 (ADC) の場合のデフォルト値
        p0 = 0.0
        p1 = 1.0
    
    else:
        # AorT=2 やその他の場合のデフォルト値
        p0 = 0.0
        p1 = 0.0

    # p0, p1 をタプルで返す
    return p0, p1

# --- スクリプト本体 (変更不要) ---

def format_line(CId, PlId, SegId, AorT, UorD, p0, p1):
    """
    ファイルに書き込む行をフォーマットする。
    区切り文字はタブ文字 (\t) を使用します。
    """
    # 指数表記 (e) で統一
    return f"{CId}\t{PlId}\t{SegId}\t{AorT}\t{UorD}\t{p0:e}\t{p1:e}"
    
    # Pythonに表記を選ばせる場合 (g)
    # return f"{CId}\t{PlId}\t{SegId}\t{AorT}\t{UorD}\t{p0:g}\t{p1:g}"


def main():
    print(f"パラメータファイル {OUTPUT_FILE} を生成します...")
    
    try:
        with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
            # 1. ヘッダーを書き込む
            f.write("#CId\tPlId\tSegId\tAorT\tUorD\tp0\tp1\n")
            
            # 2. 設定に基づいてループ処理
            # CIdでループ
            for CId, seg_count in DETECTOR_CONFIG.items():
                
                # PlIdでループ
                for PlId in PLID_LIST:
                    
                    # (AorT, UorD) の組み合わせリストを取得
                    combinations = AT_UD_COMBINATIONS_BY_CID.get(CId, AT_UD_COMBINATIONS_BY_CID["default"])
                    
                    # (AorT, UorD) の組み合わせでループ
                    for AorT, UorD in combinations:
                    
                        # SegIdでループ
                        for SegId in range(seg_count):
                            
                            # 3. 初期値を取得
                            p0, p1 = get_initial_params(CId, PlId, SegId, AorT, UorD)
                            
                            # 4. 行をフォーマット
                            line = format_line(CId, PlId, SegId, AorT, UorD, p0, p1)
                            
                            # 5. ファイルに書き込み
                            f.write(line + "\n")
                                
        print(f"ファイル {OUTPUT_FILE} の生成が完了しました。")
        print(f"パス: {os.path.abspath(OUTPUT_FILE)}")

    except IOError as e:
        print(f"ファイルの書き込み中にエラーが発生しました: {e}")
    except Exception as e:
        print(f"予期せぬエラーが発生しました: {e}")

# スクリプトとして実行された場合のみmain()を呼び出す
if __name__ == "__main__":
    main()