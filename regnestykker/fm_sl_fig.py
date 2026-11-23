
# ---- fm_sl_fig.py: picture of one inline, near | mid | far, 1240-1620 ms, each panel scaled to the 1316/1332 event ----
# Run:  cat sl_common.py fm_sl_fig.py | ssh linux-pc '<geoscience python> -I -'
from PIL import Image, ImageDraw

S = stacks()
IL = 124
T_A, T_B = 1240, 1620
i0, i1 = T_A // 2, T_B // 2
ZOOM = 3
panels = []
for k in ("nea", "mid", "far"):
    v = np.asarray(S[k][IL - 2:IL + 3, :, i0:i1], np.float64).mean(0).T          # (time, crossline), 5 inlines averaged
    ref = np.sqrt(np.mean(np.asarray(S[k][IL - 2:IL + 3, :, 650:673], np.float64) ** 2))   # 1300-1346 ms
    x = np.clip(v / (3.0 * ref), -1, 1)
    rgb = np.stack([np.where(x > 0, 1, 1 + x), 1 - np.abs(x), np.where(x < 0, 1, 1 - x)], -1)
    img = Image.fromarray((rgb * 255).astype(np.uint8)).resize((v.shape[1], v.shape[0] * ZOOM), Image.BILINEAR)
    panels.append(img)
gap, top, left = 14, 26, 46
Wd = left + sum(p.width for p in panels) + gap * 2 + 8
Ht = top + panels[0].height + 8
canvas = Image.new("RGB", (Wd, Ht), (255, 255, 255))
d = ImageDraw.Draw(canvas)
x = left
for p, name in zip(panels, ("near (ca. 4-18 grader)", "mid (ca. 17-36 grader)", "far (ca. 35-48 grader)")):
    canvas.paste(p, (x, top))
    d.text((x + 6, 6), name, fill=(0, 0, 0))
    x += p.width + gap
for tm in range(1250, 1620, 50):
    y = top + (tm - T_A) // 2 * ZOOM
    d.line([(left - 6, y), (left - 1, y)], fill=(0, 0, 0))
    d.text((4, y - 5), f"{tm}", fill=(0, 0, 0))
for tm, col in ((1316, (0, 120, 0)), (1456, (200, 0, 160))):
    y = top + (tm - T_A) // 2 * ZOOM
    xx = left
    for p in panels:
        d.line([(xx, y), (xx + 10, y)], fill=col, width=2)
        d.line([(xx + p.width - 10, y), (xx + p.width, y)], fill=col, width=2)
        xx += p.width + gap
out = os.path.join(AVO, "sleipner_near_mid_far_inline.png")
canvas.save(out)
print("written", out, canvas.size)
