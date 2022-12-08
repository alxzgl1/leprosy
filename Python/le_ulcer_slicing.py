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
nLowpassFilterParameter = 0.01 # 0.01 (default), 0.05

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




#%%

# colour difference
H = (np.abs(H[:, :, 1] - H[:, :, 2]) < dGB / 255) & (2 * H[:, :, 0] - H[:, :, 1] - H[:, :, 2] > d2RGB / 255)

#%%

plt.imshow(H, cmap='gray')
plt.axis('off')
plt.show()