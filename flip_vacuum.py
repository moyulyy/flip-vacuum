# -*- coding: utf-8 -*-
"""
flip_vacuum.py

把 POSCAR / CONTCAR / *.vasp / *.cif 模型沿 c 方向翻转，
也就是把分数坐标做 z -> 1 - z 变换（真空层从一侧换到另一侧）。

可选同时翻转 a、b 方向（x -> 1-x, y -> 1-y）。

图形界面用鼠标选文件即可，默认输出在输入文件同一目录下，
文件名加后缀 `_flip`（例如 POSCAR -> POSCAR_flip，slab.cif -> slab_flip.cif）。

也支持命令行直接调用：
    python flip_vacuum.py POSCAR slab.cif
    python flip_vacuum.py --flip-ab POSCAR
"""

import os
import sys
import argparse

# ---------------------------------------------------------------------------
# 基础数学（不依赖 numpy）
# ---------------------------------------------------------------------------

def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def _inv3(m):
    d = _det3(m)
    if abs(d) < 1e-12:
        raise ValueError("晶格矩阵奇异，无法求逆")
    inv = [[0.0] * 3 for _ in range(3)]
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d
    return inv


def _mat_vec(m, v):
    return [sum(m[i][j] * v[j] for j in range(3)) for i in range(3)]


def _flip_frac(frac, flip_a=False, flip_b=False, flip_c=True):
    """分数坐标翻转：x -> 1-x 等，并回卷到 [0,1)。"""
    out = []
    for f in frac:
        x, y, z = f[0], f[1], f[2]
        if flip_a:
            x = 1.0 - x
        if flip_b:
            y = 1.0 - y
        if flip_c:
            z = 1.0 - z
        xyz = [x, y, z]
        for k in range(3):
            xyz[k] = xyz[k] % 1.0
            if abs(xyz[k] - 1.0) < 1e-9 or abs(xyz[k]) < 1e-9:
                xyz[k] = 0.0
        out.append(xyz)
    return out


# ---------------------------------------------------------------------------
# POSCAR / CONTCAR 读写（保留原有格式信息）
# ---------------------------------------------------------------------------

def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class PosData(object):
    def __init__(self):
        self.comment = ""
        self.scale = 1.0
        self.lattice = None          # 3x3，已乘过 scale
        self.species = None          # list[str] 或 None
        self.counts = None           # list[int]
        self.selective = False
        self.direct = True
        self.frac = None             # list[[x,y,z]]
        self.flags = None            # list[[t1,t2,t3]] 或 None


def read_poscar(path):
    with open(path, "r", errors="ignore") as fh:
        lines = [ln.rstrip("\n").rstrip("\r") for ln in fh.readlines()]

    if len(lines) < 8:
        raise ValueError("文件行数太少，不像 POSCAR: %s" % path)

    data = PosData()
    data.comment = lines[0]

    # 第 2 行：缩放系数（可能是 1 个数，也可能是 3 个数）
    scale_tokens = lines[1].split()
    scale_nums = [float(t) for t in scale_tokens]
    if len(scale_nums) == 1:
        scale = scale_nums[0]
        lattice = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
        if scale < 0:
            # 负数表示目标体积
            vol = abs(_det3(lattice))
            if vol < 1e-12:
                raise ValueError("晶格体积为 0，无法处理负缩放系数")
            factor = (abs(scale) / vol) ** (1.0 / 3.0)
            lattice = [[v * factor for v in row] for row in lattice]
            data.scale = 1.0
        else:
            lattice = [[v * scale for v in row] for row in lattice]
            data.scale = 1.0
    else:
        # 三个数：分别作用于三条晶格矢量
        lattice = []
        for i in (0, 1, 2):
            row = [float(x) for x in lines[2 + i].split()[:3]]
            lattice.append([v * scale_nums[i] for v in row])
        data.scale = 1.0
    data.lattice = lattice

    idx = 5
    if not _is_number(lines[idx].split()[0]):
        data.species = lines[idx].split()
        idx += 1
    data.counts = [int(t) for t in lines[idx].split()]
    idx += 1

    if lines[idx].strip() and lines[idx].strip()[0] in ("s", "S"):
        data.selective = True
        idx += 1

    mode = lines[idx].strip()
    data.direct = mode[:1].lower() == "d"
    idx += 1

    natoms = sum(data.counts)
    if data.species is None:
        # 没有元素名，就从注释行猜（常见写法：Fe Ni ...）
        tokens = [t for t in data.comment.replace(",", " ").split() if t]
        if len(tokens) >= len(data.counts):
            data.species = tokens[-len(data.counts):]
        else:
            data.species = ["X%d" % (i + 1) for i in range(len(data.counts))]

    inv = _inv3(lattice)
    frac = []
    flags = []
    for i in range(natoms):
        if idx + i >= len(lines):
            raise ValueError("原子坐标行数不足")
        toks = lines[idx + i].split()
        coords = [float(toks[0]), float(toks[1]), float(toks[2])]
        if not data.direct:
            coords = _mat_vec(inv, coords)
        frac.append([c % 1.0 for c in coords])
        if data.selective and len(toks) >= 6:
            flags.append(toks[3:6])
        else:
            flags.append(None)
    data.frac = frac
    data.flags = flags
    return data


def write_poscar(data, path, note=None):
    lat = data.lattice
    comment = data.comment
    if note:
        comment = comment + "  " + note

    out = [comment, "1.0"]
    for row in lat:
        out.append("  %20.12f %20.12f %20.12f" % (row[0], row[1], row[2]))
    if data.species:
        out.append("  " + "  ".join(str(s) for s in data.species))
    out.append("  " + "  ".join(str(c) for c in data.counts))
    if data.selective:
        out.append("Selective dynamics")
    out.append("Direct")

    for i, f in enumerate(data.frac):
        line = "  %19.12f  %19.12f  %19.12f" % (f[0], f[1], f[2])
        if data.selective:
            fl = data.flags[i] or ["T", "T", "T"]
            line += "  " + " ".join(str(x) for x in fl)
        out.append(line)

    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# CIF 读写（用 ASE 展开对称性后按 P1 写出）
# ---------------------------------------------------------------------------

def process_cif(in_path, out_path, flip_a=False, flip_b=False, flip_c=True):
    from ase.io import read, write

    frames = read(in_path, index=":")
    if not isinstance(frames, list):
        frames = [frames]

    new_frames = []
    for atoms in frames:
        a = atoms.copy()
        frac = _flip_frac(a.get_scaled_positions(wrap=False),
                          flip_a=flip_a, flip_b=flip_b, flip_c=flip_c)
        a.set_scaled_positions(frac)
        new_frames.append(a)

    write(out_path, new_frames, format="cif")


# ---------------------------------------------------------------------------
# 统一入口
# ---------------------------------------------------------------------------

def is_cif(path):
    return os.path.splitext(path)[1].lower() == ".cif"


def default_output_path(in_path):
    base, ext = os.path.splitext(in_path)
    if ext:
        return base + "_flip" + ext
    return in_path + "_flip"


def process_file(in_path, out_path=None, flip_a=False, flip_b=False,
                 flip_c=True, log=print):
    in_path = os.path.abspath(in_path)
    if out_path is None:
        out_path = default_output_path(in_path)

    if flip_c:
        axes = []
        if flip_a:
            axes.append("a")
        if flip_b:
            axes.append("b")
        axes.append("c")
        note = "(%s-flipped)" % "".join(axes)
    else:
        note = "(a/b-flipped)"

    if is_cif(in_path):
        process_cif(in_path, out_path, flip_a=flip_a,
                    flip_b=flip_b, flip_c=flip_c)
    else:
        data = read_poscar(in_path)
        data.frac = _flip_frac(data.frac, flip_a=flip_a,
                               flip_b=flip_b, flip_c=flip_c)
        write_poscar(data, out_path, note=note)

    log("已处理: %s" % in_path)
    log("输出到: %s" % out_path)
    return out_path


# ---------------------------------------------------------------------------
# 图形界面（iOS 弹窗风格）
# ---------------------------------------------------------------------------

def run_gui():
    import tkinter as tk
    from tkinter import filedialog, messagebox

    W, H = 560, 620          # 弹窗尺寸
    R = 20                   # 主卡片圆角
    M = 20                   # 左右留白
    TRANS = "#FF00FF"        # 透明色（用于挖圆角）

    BG, CARD, SEP = "#F2F2F7", "#FFFFFF", "#E4E4E9"
    TEXT, SUB, FAINT = "#1C1C1E", "#8A8A8E", "#C7C7CC"
    BLUE, BLUE_H, BLUE_P = "#007AFF", "#1F8CFF", "#0067D6"
    RED = "#FF3B30"
    HAIR = "#D9D9DE"

    # 中文优先使用雅黑 UI
    try:
        _probe = tk.Tk()
        _families = set(_probe.tk.call("font", "families"))
        _probe.destroy()
    except tk.TclError:
        _families = set()
    FONT = "Segoe UI"
    for _f in ("Microsoft YaHei UI", "Microsoft YaHei", "PingFang SC",
               "Noto Sans SC", "Segoe UI"):
        if _f in _families:
            FONT = _f
            break

    def rr(cv, x1, y1, x2, y2, r, **kw):
        """在 Canvas 上画圆角矩形（平滑多边形）。"""
        pts = [x1 + r, y1, x2 - r, y1, x2, y1, x2, y1 + r,
               x2, y2 - r, x2, y2, x2 - r, y2, x1 + r, y2,
               x1, y2, x1, y2 - r, x1, y1 + r, x1, y1]
        return cv.create_polygon(pts, smooth=True, **kw)

    def mix(c1, c2, t):
        a = tuple(int(c1[i:i + 2], 16) for i in (1, 3, 5))
        b = tuple(int(c2[i:i + 2], 16) for i in (1, 3, 5))
        return "#%02X%02X%02X" % tuple(
            int(a[i] + (b[i] - a[i]) * t) for i in range(3))

    class RoundButton(tk.Canvas):
        """圆角按钮（胶囊或圆角矩形）。"""
        def __init__(self, master, text, command=None, width=120, height=32,
                     radius=None, fill=BLUE, hover=None, press=None,
                     fg="#FFFFFF", bg=CARD, font=None):
            tk.Canvas.__init__(self, master, width=width, height=height,
                               highlightthickness=0, bd=0, bg=bg)
            r = height / 2.0 if radius is None else radius
            self._cmd = command
            self._colors = (fill,
                            hover or mix(fill, "#FFFFFF", 0.12),
                            press or mix(fill, "#000000", 0.15))
            self._shape = rr(self, 0.5, 0.5, width - 0.5, height - 0.5, r,
                             fill=fill, outline="")
            self._text = self.create_text(width / 2.0, height / 2.0,
                                          text=text, fill=fg,
                                          font=font or (FONT, 10, "bold"))
            self.configure(cursor="hand2")
            self.bind("<Enter>", lambda e: self._paint(1))
            self.bind("<Leave>", lambda e: self._paint(0))
            self.bind("<ButtonPress-1>", lambda e: self._paint(2))
            self.bind("<ButtonRelease-1>", self._release)

        def _paint(self, i):
            self.itemconfigure(self._shape, fill=self._colors[i])

        def set_text(self, t):
            self.itemconfigure(self._text, text=t)

        def _release(self, event):
            self._paint(0)
            if (0 <= event.x <= self.winfo_width()
                    and 0 <= event.y <= self.winfo_height() and self._cmd):
                self._cmd()

    class TextButton(tk.Label):
        """纯文字按钮（次要操作）。"""
        def __init__(self, master, text, command=None, fg=BLUE, bg=CARD,
                     font=None):
            tk.Label.__init__(self, master, text=text, fg=fg, bg=bg,
                              font=font or (FONT, 10), cursor="hand2")
            self._fg = fg
            self.bind("<Button-1>", lambda e: command() if command else None)
            self.bind("<Enter>", lambda e: self.configure(
                fg=mix(self._fg, "#FFFFFF", 0.35)))
            self.bind("<Leave>", lambda e: self.configure(fg=self._fg))

    class Toggle(tk.Canvas):
        """iOS 开关。"""
        TW, TH = 50, 30

        def __init__(self, master, value=False, command=None, bg=CARD):
            tk.Canvas.__init__(self, master, width=self.TW, height=self.TH,
                               highlightthickness=0, bd=0, bg=bg)
            self.value = bool(value)
            self._cmd = command
            self._pos = 1.0 if self.value else 0.0
            self._job = None
            self._track = rr(self, 0.5, 0.5, self.TW - 0.5, self.TH - 0.5,
                             (self.TH - 1) / 2.0,
                             fill=self._color(self._pos), outline="")
            self._knob = self.create_oval(0, 0, 0, 0, fill="#FFFFFF",
                                          outline="")
            self._render()
            self.configure(cursor="hand2")
            self.bind("<Button-1>", self.toggle)

        @staticmethod
        def _color(p):
            a, b = (229, 229, 234), (52, 199, 89)
            return "#%02X%02X%02X" % tuple(
                int(a[i] + (b[i] - a[i]) * p) for i in range(3))

        def _render(self):
            d = self.TH - 4
            x1 = 2 + self._pos * (self.TW - d - 4)
            self.coords(self._knob, x1, 2, x1 + d, 2 + d)
            self.itemconfigure(self._track, fill=self._color(self._pos))

        def toggle(self, event=None):
            self.value = not self.value
            target = 1.0 if self.value else 0.0
            start = self._pos
            if self._job:
                self.after_cancel(self._job)
                self._job = None

            def step(i, total=8):
                t = i / float(total)
                e = 1 - (1 - t) ** 3
                self._pos = start + (target - start) * e
                self._render()
                if i < total:
                    self._job = self.after(10, step, i + 1)
                else:
                    self._job = None

            step(1)
            if self._cmd:
                self._cmd(self.value)

    # ----- 窗口 -------------------------------------------------------
    root = tk.Tk()
    root.title("翻转真空层方向")
    root.overrideredirect(True)          # 去掉系统标题栏 -> 弹窗感
    root.resizable(False, False)

    transparent = True
    try:
        root.attributes("-transparentcolor", TRANS)
    except tk.TclError:
        transparent = False
    base = TRANS if transparent else BG
    root.configure(bg=base)

    px = int((root.winfo_screenwidth() - W) / 2)
    py = int((root.winfo_screenheight() - H) / 3)
    root.geometry("%dx%d+%d+%d" % (W, H, px, py))

    bg = tk.Canvas(root, width=W, height=H, bd=0, highlightthickness=0,
                   bg=base)
    bg.place(x=0, y=0)

    # 主卡片：底部一层淡淡的投影，避免出现"灰边圈"
    rr(bg, 2, 5, W - 2, H - 2, R, fill="#E2E2E7", outline="")
    rr(bg, 1, 1, W - 2, H - 4, R, fill=BG, outline=HAIR)

    # ----- 顶部标题 ---------------------------------------------------
    bg.create_text(W / 2, 36, text="翻转真空层方向", fill=TEXT,
                   font=(FONT, 15, "bold"))
    bg.create_text(W / 2, 58, text="POSCAR / CIF　·　c 轴真空层换边",
                   fill=SUB, font=(FONT, 9))
    bg.create_line(M + 8, 78, W - M - 8, 78, fill=SEP)

    close = bg.create_oval(W - 44, 24, W - 22, 46, fill="", outline="")
    close_x = bg.create_text(W - 33, 35, text="✕", fill="#9A9AA0",
                             font=(FONT, 10))
    for item in (close, close_x):
        bg.tag_bind(item, "<Button-1>", lambda e: root.destroy())
        bg.tag_bind(item, "<Enter>", lambda e: (
            bg.itemconfigure(close, fill="#E9E9EB"),
            bg.itemconfigure(close_x, fill="#6C6C70")))
        bg.tag_bind(item, "<Leave>", lambda e: (
            bg.itemconfigure(close, fill=""),
            bg.itemconfigure(close_x, fill="#9A9AA0")))

    # ----- 结构文件 ---------------------------------------------------
    bg.create_text(M + 8, 100, text="结构文件", anchor="w", fill=SUB,
                   font=(FONT, 10))
    cnt_item = bg.create_text(W - M - 8, 100, text="0 个文件", anchor="e",
                              fill=SUB, font=(FONT, 10))
    rr(bg, M, 112, W - M, 306, 14, fill=CARD, outline=SEP)

    listbox = tk.Listbox(root, selectmode="extended", bd=0,
                         highlightthickness=0, activestyle="none",
                         bg=CARD, fg=TEXT, font=(FONT, 10),
                         selectbackground=BLUE, selectforeground="#FFFFFF")
    listbox.place(x=M + 12, y=122, width=W - 2 * M - 24, height=144)

    placeholder = tk.Label(root, text="尚未添加文件\n点击「添加文件」选择 POSCAR / CIF",
                           bg=CARD, fg=FAINT, font=(FONT, 10), justify="center")
    ph = dict(x=M + 12, y=170, width=W - 2 * M - 24, height=48)
    placeholder.place(**ph)

    def show_placeholder():
        placeholder.place(**ph)
        placeholder.lift()

    def hide_placeholder():
        placeholder.place_forget()

    files = []

    def log(msg):
        log_box.configure(state="normal")
        log_box.insert("end", msg + "\n")
        log_box.see("end")
        log_box.configure(state="disabled")
        root.update_idletasks()

    def add_files():
        root.attributes("-topmost", True)
        paths = filedialog.askopenfilenames(
            title="选择结构文件",
            filetypes=[("结构文件",
                        "*.cif *.vasp *.poscar *.contcar POSCAR CONTCAR"),
                       ("CIF", "*.cif"),
                       ("VASP", "*.vasp POSCAR CONTCAR"),
                       ("所有文件", "*.*")])
        for p in paths:
            if p not in files:
                files.append(p)
                listbox.insert("end", "  " + os.path.basename(p)
                               + "    —    " + os.path.dirname(p))
        if files:
            hide_placeholder()
        bg.itemconfigure(cnt_item, text="%d 个文件" % len(files))
        if paths:
            log("已添加 %d 个文件，共 %d 个。" % (len(paths), len(files)))

    def clear_files():
        del files[:]
        listbox.delete(0, "end")
        show_placeholder()
        bg.itemconfigure(cnt_item, text="0 个文件")
        log("已清空文件列表。")

    RoundButton(root, "＋ 添加文件", command=add_files, width=118,
                height=32, radius=9, bg=CARD,
                font=(FONT, 10)).place(x=M + 12, y=270)
    TextButton(root, "清空列表", command=clear_files, fg=RED,
               font=(FONT, 10)).place(x=M + 142, y=277)
    bg.create_text(W - M - 12, 286, text="支持 .cif / POSCAR / CONTCAR",
                   anchor="e", fill=FAINT, font=(FONT, 8))

    # ----- 选项 -------------------------------------------------------
    bg.create_text(M + 8, 328, text="选项", anchor="w", fill=SUB,
                   font=(FONT, 10))
    rr(bg, M, 340, W - M, 454, 14, fill=CARD, outline=SEP)

    toggle_c = Toggle(root, value=True)
    toggle_c.place(x=W - M - 12 - Toggle.TW, y=354)
    bg.create_text(M + 20, 368, text="翻转 c 方向（z → 1-z）", anchor="w",
                   fill=TEXT, font=(FONT, 11))
    bg.create_text(M + 20, 385, text="把真空层从一侧换到另一侧", anchor="w",
                   fill=SUB, font=(FONT, 8))
    bg.create_line(M + 20, 397, W - M - 20, 397, fill=SEP)

    toggle_ab = Toggle(root, value=False)
    toggle_ab.place(x=W - M - 12 - Toggle.TW, y=411)
    bg.create_text(M + 20, 425, text="同时翻转 a、b 方向", anchor="w",
                   fill=TEXT, font=(FONT, 11))
    bg.create_text(M + 20, 442, text="x → 1-x，y → 1-y（一般不需要）",
                   anchor="w", fill=SUB, font=(FONT, 8))

    # ----- 运行日志 ---------------------------------------------------
    bg.create_text(M + 8, 476, text="运行日志", anchor="w", fill=SUB,
                   font=(FONT, 10))
    rr(bg, M, 488, W - M, 542, 14, fill=CARD, outline=SEP)

    log_box = tk.Text(root, bd=0, highlightthickness=0, bg=CARD,
                      fg="#48484A", font=("Consolas", 8), wrap="none",
                      state="disabled")
    log_box.place(x=M + 12, y=496, width=W - 2 * M - 24, height=38)
    log("就绪：点击「添加文件」选择 POSCAR / CIF。")

    # ----- 主按钮 -----------------------------------------------------
    def run():
        if not files:
            messagebox.showwarning("提示", "请先选择文件！", parent=root)
            return
        if not toggle_c.value and not toggle_ab.value:
            messagebox.showwarning("提示", "至少选择一个翻转方向！", parent=root)
            return
        ok = 0
        for p in list(files):
            try:
                process_file(p, flip_a=toggle_ab.value,
                             flip_b=toggle_ab.value, flip_c=toggle_c.value,
                             log=log)
                ok += 1
            except Exception as exc:  # noqa: BLE001
                log("!!! 处理失败: %s  ->  %s" % (os.path.basename(p), exc))
        log("完成：成功 %d / %d。" % (ok, len(files)))

    RoundButton(root, "开始处理", command=run, width=W - 2 * M, height=44,
                radius=14, font=(FONT, 11, "bold")).place(x=M, y=558)

    # ----- 拖动 / 关闭 -------------------------------------------------
    def start_move(event):
        root._off = (event.x_root - root.winfo_x(),
                     event.y_root - root.winfo_y())

    def do_move(event):
        root.geometry("+%d+%d" % (event.x_root - root._off[0],
                                  event.y_root - root._off[1]))

    bg.bind("<ButtonPress-1>", start_move)
    bg.bind("<B1-Motion>", do_move)
    root.bind("<Escape>", lambda e: root.destroy())
    root.attributes("-topmost", True)
    root.after(500, lambda: root.attributes("-topmost", False))
    root.focus_force()

    root.mainloop()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description="翻转 POSCAR/CIF 的 c 方向（真空层换边）。")
    parser.add_argument("files", nargs="*", help="输入文件")
    parser.add_argument("--flip-ab", action="store_true",
                        help="同时翻转 a、b 方向")
    parser.add_argument("-o", "--output", default=None,
                        help="输出文件（仅单个输入时有效）")
    args = parser.parse_args(argv)

    if not args.files:
        run_gui()
        return 0

    for p in args.files:
        out = args.output if (args.output and len(args.files) == 1) else None
        process_file(p, out, flip_a=args.flip_ab, flip_b=args.flip_ab,
                     flip_c=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
