import cv2
import numpy as np

fps = 2
n_steps = 50

video = cv2.VideoWriter(
    "combined_video.mp4",
    cv2.VideoWriter_fourcc(*"mp4v"),
    fps,
    (1080, 1080)
)

def load(path, size):
    img = cv2.imread(path)
    if img is None:
        raise FileNotFoundError(path)
    return cv2.resize(img, size)

for i in range(1, n_steps + 1):

    # Left
    keff_burnup = load(f"keff_burnup_{i:03d}.png", (360, 540))
    keff_time   = load(f"keff_time_{i:03d}.png",   (360, 540))
    left = np.vstack((keff_burnup, keff_time))

    # Right
    flux = load(f"flux_xz_step{i}.png", (360, 1080))
    heat = load(f"heat_xz_step{i}.png", (360, 1080))

    # Combined 
    frame = np.hstack((left, flux, heat))

    video.write(frame)

video.release()
print("✅ Video generated successfully")

