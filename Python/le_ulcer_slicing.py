import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

# settings
bCutImage = 1
nImageHalfWidth = 350

aFilter = 'lowpass' # 'median' (default), 'lowpass'

dGB = 10 # 10 (default), green - blue difference | detect ulcer
d2RGB = 80 # 80 (default), 2 * red - green - blue | exclude dark background

nMedianFilterArea = [31, 31] # [32, 32] (default)
nLowpassFilterParameter = 0.02 # 0.02 (default), 0.05

bNeighbourBlobThreshold = 10

# smoothing filter
fb, fa = signal.butter(4, 0.01, 'low')

# LP filter
lfb, lfa = signal.butter(4, nLowpassFilterParameter, 'low')

# load image
aFilename = 'd:/data/leprosy/TABLE_Aranz_Image/S-01/S-1_20200921.png'
I = mpimg.imread(aFilename)

# cut image
nWidth = np.shape(I)[1]
nHeight = np.shape(I)[0]
if bCutImage == 1:
  d = nImageHalfWidth
  x = int(nWidth / 2)
  y = int(nHeight / 2)
  I = I[(y - d):(y + d), (x - d):(x + d), :]
  
# median filter
if aFilter == 'median':
  D = nMedianFilterArea
  J_R = signal.medfilt2d(I[:, :, 0], D)
  J_G = signal.medfilt2d(I[:, :, 1], D)
  J_B = signal.medfilt2d(I[:, :, 2], D)
  J = np.stack([J_R, J_G, J_B], axis=2)
  H = J.astype(float)
  S = J.astype(float)
elif aFilter == 'lowpass':
  J = I.astype(float)
  J_R = J[:, :, 0]
  J_G = J[:, :, 1]
  J_B = J[:, :, 2]
  # smooth rows
  J_R = signal.filtfilt(lfb, lfa, J_R)
  J_G = signal.filtfilt(lfb, lfa, J_G)
  J_B = signal.filtfilt(lfb, lfa, J_B)
  # smooth cols
  J_R = signal.filtfilt(lfb, lfa, J_R.T)
  J_G = signal.filtfilt(lfb, lfa, J_G.T)
  J_B = signal.filtfilt(lfb, lfa, J_B.T)
  # init
  J = np.stack([J_R.T, J_G.T, J_B.T], axis=2)
  H = J.astype(float)
  S = J.astype(float)

# colour difference
H = (np.abs(H[:, :, 1] - H[:, :, 2]) < dGB / 255) & (2 * H[:, :, 0] - H[:, :, 1] - H[:, :, 2] > d2RGB / 255)

# circle mask
x = np.sum(H, axis=0)
y = np.sum(H, axis=1)
x = signal.filtfilt(fb, fa, x)
y = signal.filtfilt(fb, fa, y)
cx = 0
cy = 0
# range circles
pR = np.arange(5, nImageHalfWidth, 5)
nR = len(pR)
S = np.zeros(nR)
W = np.arange(-nImageHalfWidth, nImageHalfWidth, 1)
W = np.tile(W, [len(W), 1])
for iR in range(0, nR):
  R = pR[iR]
  s = np.sqrt((W - cx) ** 2 + (W.T - cy) ** 2) < R
  s = H * s
  S[iR] = np.sum(s)
xR = bNeighbourBlobThreshold # threshold
dS = np.diff(S)
hS = (dS > xR).astype(float)
i = np.squeeze(np.where(hS > 0.0))
i = i[0] # TODO: if isempty()
if i < 20: # threshold
  iR = np.squeeze(np.where(hS[i:] < 1.0))
  iR = iR[0] + i
  # if not iR: # TODO: if isempty()
  #   iR = len(pR)
else:
  iR = 0
R = pR[iR]

# apply circle mask
M = np.sqrt((W - cx) ** 2 + (W.T - cy) ** 2) < R
H = M * H
  
# fill holes
# see, https://stackoverflow.com/questions/36294025/python-equivalent-to-matlab-funciton-imfill-for-grayscale

# show image
plt.subplot(1, 2, 1)
plt.imshow(I)
plt.axis('off')
plt.subplot(1, 2, 2)
plt.imshow(H, cmap='gray')
plt.axis('off')
