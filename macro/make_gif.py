#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
指定ディレクトリ内の画像ファイルをファイル名の自然順（"img2" < "img10" になる順）に
結合し、GIF アニメーションを作成する汎用スクリプト。

画像の中身には関知しない。tpc_stage_display.py の --save で .jpg 出力した
連番フレーム（白背景・不透明、容量抑制用）をまとめて GIF 化する用途を主に想定しているが、
それ以外のディレクトリ・拡張子（--pattern で指定）にも使える。
透過 PNG など、アルファチャンネル付き画像を渡した場合は白背景に合成してから結合する。

Usage:
  python3 make_gif.py <image_dir> -o out.gif
  python3 make_gif.py <image_dir> -o out.gif --pattern "*.png" --duration 300 --loop 0
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import List

from PIL import Image

DEFAULT_PATTERN = "*.jpg"
DEFAULT_DURATION_MS = 500
DEFAULT_LOOP = 0  # 0 = 無限ループ


def _natural_sort_key(path: str):
    """"frame2.jpg" が "frame10.jpg" より前に来るよう、数字部分を int として比較する。"""
    name = os.path.basename(path)
    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", name)]


def list_images(image_dir: str, pattern: str) -> List[str]:
    paths = glob.glob(os.path.join(image_dir, pattern))
    paths.sort(key=_natural_sort_key)
    return paths


def _load_frame_rgb(path: str) -> Image.Image:
    """アルファチャンネル付き画像（透過PNG等）は白背景に合成してから RGB にする。"""
    img = Image.open(path)
    has_alpha = img.mode in ("RGBA", "LA") or (img.mode == "P" and "transparency" in img.info)
    if not has_alpha:
        return img.convert("RGB")
    img = img.convert("RGBA")
    bg = Image.new("RGB", img.size, (255, 255, 255))
    bg.paste(img, mask=img.split()[-1])
    return bg


def make_gif(
    image_dir: str,
    output: str,
    pattern: str = DEFAULT_PATTERN,
    duration: int = DEFAULT_DURATION_MS,
    loop: int = DEFAULT_LOOP,
) -> str:
    paths = list_images(image_dir, pattern)
    if not paths:
        raise FileNotFoundError(f"'{pattern}' に一致する画像が {image_dir} に見つかりません")
    frames = [_load_frame_rgb(p) for p in paths]
    frames[0].save(
        output,
        save_all=True,
        append_images=frames[1:],
        duration=duration,
        loop=loop,
    )
    return output


def _parse_cli(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="ディレクトリ内の画像をファイル名の自然順に結合してGIFアニメーションを作成する"
    )
    p.add_argument("image_dir", help="画像ファイルが入っているディレクトリ")
    p.add_argument("-o", "--output", default="output.gif", help="出力GIFパス（既定: output.gif）")
    p.add_argument(
        "--pattern", default=DEFAULT_PATTERN,
        help=f"入力画像の glob パターン（既定: {DEFAULT_PATTERN}）",
    )
    p.add_argument(
        "--duration", type=int, default=DEFAULT_DURATION_MS,
        help=f"1フレームあたりの表示時間 [ms]（既定: {DEFAULT_DURATION_MS}）",
    )
    p.add_argument(
        "--loop", type=int, default=DEFAULT_LOOP,
        help=f"ループ回数（0=無限ループ、既定: {DEFAULT_LOOP}）",
    )
    return p.parse_args(argv)


def main(argv: List[str]) -> None:
    args = _parse_cli(argv)
    paths = list_images(args.image_dir, args.pattern)
    if not paths:
        print(f"Error: '{args.pattern}' に一致する画像が {args.image_dir} に見つかりません", file=sys.stderr)
        sys.exit(1)
    print(f"{len(paths)} 枚の画像を結合します:")
    for path in paths:
        print(f"  {path}")
    out = make_gif(
        args.image_dir, args.output,
        pattern=args.pattern, duration=args.duration, loop=args.loop,
    )
    print(f"Saved: {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
